classdef sindy_fast
    methods(Static)
        
        function a_dot_dns = eval_adot_spline(a_dns)

            [dns_time, n_pod, n_wave] = size(a_dns);
            
            a_dot_dns = zeros(dns_time, n_pod, n_wave);            
            
            disp('Calculating adot from dns')

            for i_diff = 1:n_wave
                for i_pod = 1:n_pod
                    a_dns_re = spline(1:dns_time, real(a_dns(1:dns_time,i_pod,i_diff)));
                    a_dns_im = spline(1:dns_time, imag(a_dns(1:dns_time,i_pod,i_diff)));
                    
                    a_dot_dns_re = fnder(a_dns_re, 1);
                    a_dot_dns_imag = fnder(a_dns_im,1);
            
                    a_dot_dns(:,i_pod, i_diff) = ppval(a_dot_dns_re, 1:dns_time) + ...
                                                 1i*ppval(a_dot_dns_imag, 1:dns_time);
                end

                if mod(i_diff, 100) == 0
                    disp([num2str(i_diff),' / ', num2str(n_wave)])
                end
            end

        end
        
        function D = eval_D(ny, Lx, Lz, phi, n, m, nx, nz, pod_wave, dw, w)

            % This function evaluates the eddy viscousity model          

            kron = 0;
            if n == m
                kron = 1;
            end

            mean_mode = -((2*pi*nx/Lx)^2 + (2*pi*nz/Lz)^2)*kron;

            index = gp.find_wave_number_in_map(nx, nz, pod_wave);
    
            phiu_n = phi(1:3:end, n, index);
            phiv_n = phi(2:3:end, n, index);
            phiw_n = phi(3:3:end, n, index);
    
            phiu_m = phi(1:3:end, m, index);
            phiv_m = phi(2:3:end, m, index);
            phiw_m = phi(3:3:end, m, index); 

            dphiu_n = dw(:,:,1)*phiu_n;
            dphiv_n = dw(:,:,1)*phiv_n;
            dphiw_n = dw(:,:,1)*phiw_n;

            dphiu_m = dw(:,:,1)*phiu_m;
            dphiv_m = dw(:,:,1)*phiv_m;
            dphiw_m = dw(:,:,1)*phiw_m;

            visc_u_int = conj(dphiu_n).*dphiu_m;
            visc_v_int = conj(dphiv_n).*dphiv_m;
            visc_w_int = conj(dphiw_n).*dphiw_m;

            visc_u = 0;
            visc_v = 0;
            visc_w = 0;

            for i = 1:ny
                visc_u = visc_u + w(i)*visc_u_int(i);
                visc_v = visc_v + w(i)*visc_v_int(i);
                visc_w = visc_w + w(i)*visc_w_int(i);
            end
            visc = visc_u + visc_v + visc_w;

            D = -visc + mean_mode;

            % check for mean mode again
            if nx == 0 && nz == 0
                D = 0;
            end

        end 
        
        function [D, fullmode_pod] = get_D_struct(Np, phi, lin_nonlin_ind, fullmode, pod_wave)
            
            [~,~,~,~,ny,~,Lx,Lz] = modal_decomposition.read_geom;
            
            const_lin_nonlin_ind = lin_nonlin_ind; % This line saves the index map to be used in ODE45 (nothing mathematical)
            
            w = math.f_chebyshev_int_weight(ny);
            [~,dw] = math.f_chebdiff(ny, 1);
                            
            [D, ~, fullmode, fullmode_pod] = gp.prep_L_matrix(Np, fullmode);
            D = gp.find_POD_pair_conj_index_in_model(D, Np, const_lin_nonlin_ind,fullmode);
            
            for i = 1:size(D,2)
                size_n_seq = size(D(i).n_seq,2);
                nxnz = D(i).nxnz;
                for i_pod = 1:size_n_seq
                    nx = nxnz(1);
                    nz = nxnz(2);
                    D(i).coeff(i_pod) = sindy.eval_D(ny, Lx, Lz, phi, D(i).n, D(i).n_seq(i_pod), ...
                        nx, nz, pod_wave, dw, w);
                end
            end
        
        end
        
        function T = residual_galerkin_R(Np, L, N, a_dns, a_dns_dot, fullmode_pod,pod_wave)
            
            % size of input a 
            [a_t, ~, ~] = size(a_dns);

            % re-arrange a_dns
            a_dns_sindy = zeros(a_t, size(fullmode_pod,1));
            a_dns_dot_sindy = a_dns_sindy;
            for i = 1:size(fullmode_pod,1)
                nxnz = fullmode_pod(i,2:3);
                Np = fullmode_pod(i,1);
            
                I = gp.find_wave_number_in_map(nxnz(1), nxnz(2), pod_wave);
            
                a_dns_sindy(:,i) = a_dns(:, Np, I);
                a_dns_dot_sindy(:,i) = a_dns_dot(:,Np,I);
            
            end
            
            T = zeros(a_t, size(fullmode_pod,1));
            for i_t = 1:a_t
                    tmp_1 = gp.eval_adot_fast(i_t, a_dns_sindy(i_t,:)', L, N, Np);
                    T(i_t, :) =  a_dns_dot_sindy(i_t,:) - tmp_1';      
            end
        
        end

        function theta = library_galerkin_R(e, a_dns, D, pod_wave)
            
            [time, ~, ~] = size(a_dns);
            
            theta = ones(time, size(D, 2));
        
            % re-arrange a_dns into theta
            for i_t = 1:time
                for i = 1:size(D,2)
                    nxnz = D(i).nxnz;
                    Np = D(i).n;
            
                    I = gp.find_wave_number_in_map(nxnz(1),nxnz(2),pod_wave);
                            
                    theta(i_t, i) = a_dns(i_t, Np, I)*D(i).coeff(Np)*e(i_t);
                end
            end
            
         end        
        
        function Xi = sindy_solve(dof, lambda, k, Xdot, theta)

            % SINDY BY S.BRUNTON

                Xi = theta\Xdot;
            
                % dof : system degree of freedom
                for iter = 1:k
                    smallinds = (abs(Xi)<lambda);
                    Xi(smallinds) = 0;
                    for i = 1:dof
                        biginds = ~smallinds(:,i);
            %             biginds = ~smallinds(i)
                        Xi(biginds, i) = theta(:,biginds)\Xdot(:,i);
                    end
                end
            
            end

        function adot = eval_adot_galerkin_R(t, a, L, N, D, c, nw_pod)
            
            disp(['t = ', num2str(t)]);
    
            a_len = size(L,2);
    
            adot_lin = zeros(a_len,1);
            adot_nonlin = adot_lin;
            adot_eddy_visc = adot_lin;
    
            % total number of permutations of the PODs
%             perms_length = size(N(1).n_seq,1);        
            
            % linear evaluation
            conj_len = size(L(1).conj_ind_start,1);
            for i_conj = 1:conj_len
                st = L(1).conj_ind_start(i_conj);
                change_pod_ind = L(1).change_POD_ind(i_conj);
                for i_lin = change_pod_ind:st
        
                    for i_pod = 1:nw_pod
        
                        anxnz = a(L(i_lin).pod_pair(i_pod));
        
                        adot_lin(i_lin) = adot_lin(i_lin) + L(i_lin).coeff(i_pod)*anxnz;
                    end
                        
        
                    % nonlinear evaluation
        
                    i_nonlin = i_lin;
                    
                    kxkz_arr = N(i_nonlin).kxkz;
                    kxkz_len = size(kxkz_arr,1);
        
                    a_mxmz_loc = N(i_nonlin).a_mxmz_loc;
                    a_kxkz_loc = N(i_nonlin).a_kxkz_loc;
                    N_coeff_tmp = N(i_nonlin).coeff;
        
        
                    for i = 1:kxkz_len               
                        
                        m_loc = a_mxmz_loc(i, :);
                        k_loc = a_kxkz_loc(i, :);
        
                        tmp = a(m_loc).*a(k_loc);
        
                        adot_nonlin(i_nonlin) = adot_nonlin(i_nonlin) + ...
                                        N_coeff_tmp(i,:)*tmp;
                    end   
    
                    
                    % eddy viscousity evaluation (Galerkin-R)
                    i_visc = i_lin;
                    for i_pod_visc = 1:nw_pod
    
                        anxnz_visc = a(D(i_visc).pod_pair(i_pod_visc));
                        adot_eddy_visc(i_visc) = adot_eddy_visc(i_visc) + D(i_visc).coeff(i_pod_visc)*anxnz_visc;
    
                    end
    
                    % fills in e(t) and c
                    adot_eddy_visc(i_visc) = adot_eddy_visc(i_visc).*c(i_visc);
                    
                    e = abs(a(25));
    
                    adot_eddy_visc(i_visc) = e.*adot_eddy_visc(i_visc);
    
                    
                    
                end
            end


    
            adot = adot_lin + adot_nonlin + adot_eddy_visc;
    
            % force conjugation
            conj_len = size(L(1).conj_ind_start,1);
            for i_conj = 1:conj_len
                st = L(1).conj_ind_start(i_conj);
                ed = L(1).conj_ind_end(i_conj); % ed for end 
    
                change_pod_ind = L(1).change_POD_ind(i_conj);
    
                adot(st:ed) = flip(conj(adot(change_pod_ind:st-2)));
            end
    
        end
        
        function adot = eval_adot_galerkin_R_origin(t, a, L, N, D, c, nw_pod)
            
            disp(['t = ', num2str(t)]);
    
            a_len = size(L,2);
    
            adot_lin = zeros(a_len,1);
            adot_nonlin = adot_lin;
            adot_eddy_visc = adot_lin;
    
            % total number of permutations of the PODs
            perms_length = size(N(1).n_seq,1);        
            
            % linear evaluation
            for i_lin = 1:a_len
    
                for i_pod = 1:nw_pod
    
                    anxnz = a(L(i_lin).pod_pair(i_pod));
    
                    adot_lin(i_lin) = adot_lin(i_lin) + L(i_lin).coeff(i_pod)*anxnz;
                end
    
            
            
    
                % nonlinear evaluation
    
                i_nonlin = i_lin;
                
                kxkz_arr = N(i_nonlin).kxkz;
                kxkz_len = size(kxkz_arr,1);
    
                for i = 1:kxkz_len               
                    
                    for i_perms = 1:perms_length
    
                        m_loc = N(i_nonlin).a_mxmz_loc(i, i_perms);
                        k_loc = N(i_nonlin).a_kxkz_loc(i, i_perms);
    
                        adot_nonlin(i_nonlin) = adot_nonlin(i_nonlin) + ...
                            N(i_nonlin).coeff(i,i_perms)*a(m_loc)*a(k_loc);
                    end
    
                end

                
                % eddy viscousity evaluation (Galerkin-R)
                i_visc = i_lin;
                for i_pod_visc = 1:nw_pod

                    anxnz_visc = a(D(i_visc).pod_pair(i_pod_visc));
                    adot_eddy_visc(i_visc) = adot_eddy_visc(i_visc) + D(i_visc).coeff(i_pod_visc)*anxnz_visc;

                end

                % fills in e(t) and c
                adot_eddy_visc(i_visc) = adot_eddy_visc(i_visc).*c(i_visc);
                
                e = abs(a(25));

                adot_eddy_visc(i_visc) = e.*adot_eddy_visc(i_visc);

                
            end


    
            adot = adot_lin + adot_nonlin + adot_eddy_visc;
    
            % force conjugation
            conj_len = size(L(1).conj_ind_start,1);
            for i_conj = 1:conj_len
                st = L(1).conj_ind_start(i_conj);
                ed = L(1).conj_ind_end(i_conj); % ed for end 
    
                change_pod_ind = L(1).change_POD_ind(i_conj);
    
                adot(st:ed) = flip(conj(adot(change_pod_ind:st-2)));
            end
    
        end


    end
end