classdef gp
    methods(Static)
    
    function I = find_wave_number_in_map(nx, nz, pod_wave)
    
        r = size(pod_wave,1);
        
        for i = 1:r
        
            if ((pod_wave(i,1) == nx) && (pod_wave(i,2) == nz))
                break
            end
    
        end
    
        I = i;
    
    end

    function I = find_wave_number_pod_in_map(np, nx, nz, pod_wave)
    
        r = size(pod_wave,1);
        
        for i = 1:r
        
            if ((pod_wave(i,1) == np) && (pod_wave(i,2) == nx) && (pod_wave(i,3) == nz))
                break
            end
    
        end
    
        I = i;
    
    end    

    function L = find_POD_pair_conj_index_in_model(L, nw_pod, const_lin_nonlin_ind,fullmode)

        % Nothing mathematical. Just for the ease of calculation

        for i = 1:size(fullmode,1)
        c_nxnz = L(i).nxnz;
    
        counter = 0;
    
            for ii = 1:size(fullmode,1)
                match_nxnz = L(ii).nxnz;
        
                if ((match_nxnz == c_nxnz))
                    counter = counter + 1;
                    L(i).pod_pair(counter) = ii;
                end
            end
        end       

        conj_ind = zeros(nw_pod,1);
        conj_ind(1) = (size(const_lin_nonlin_ind,2) + 1)/2; % first conjugation index

        for i = 2:nw_pod
            conj_ind(i) = conj_ind(i-1) + size(const_lin_nonlin_ind,2);
        end

        L(1).conj_ind_start = conj_ind  + 1;
        L(1).conj_ind_end = conj_ind + (size(const_lin_nonlin_ind,2) + 1)/2 - 1;

        for i = 1:size(L(1).conj_ind_end,1)
            if L(1).conj_ind_end(i) > size(fullmode,1)
                L(1).conj_ind_end(i) = []; % Check if index out of range.
            end
        end

        L(1).change_POD_ind = L(1).conj_ind_start - conj_ind(1);
        L(1).change_POD_ind(1) = L(1).change_POD_ind(1) ;

    end

    function [out, loc] = filter_out_of_range(model, wave_map)
        % this function filters out the mx = nx - kx and mz that are no in
        % the model wavenumber range
    
        [r_wave, ~] = size(wave_map);
    
        for i = 1:r_wave
            a(i) = ismember(abs(wave_map(i,:)),model,'row');
            
            if a(i) == 1
                loc(i) = i;
            end
        end
    
        wave_map = wave_map(a,:);
        out = wave_map;
        
    end    
    
    function [lin_nonlin_ind, fullmode, model] = load_predefined_model(model_name)

        if model_name == 6
            model = [ 0, 0, 1; ...
                   0, 1, 1; 0,-1, 1; ...
                   0, 2, 1; 0,-2, 1; ...
                   1, 0, 1;-1, 0, 1; ...
                   1, 1, 1;-1,-1, 1; ...
                   1,-1, 1;-1, 1, 1];
        elseif model_name == 6.1
            model = [ 0, 0, 1; ...
                   0, 1, 1; 0,-1, 1; ...
                   0, 2, 1; 0,-2, 1; ...
                   1, 0, 1;-1, 0, 1; ...
                   1, 1, 1;-1,-1, 1; ...
                   1,-1, 1;-1, 1, 1; ...
                   1, 0, 2;-1, 0, 2; ...
                   1, 1, 2;-1,-1, 2; ...
                   1,-1, 2;-1, 1, 2];                               
        elseif model_name == 9
                model = [ 0, 0, 1; 0, 1, 1; ...
                   1, 0, 1; 1, 1, 1; ...
                   1,-1, 1; 0, 2, 1; ...
                   1, 2, 1; 1,-2, 1; ...
                   0, 3, 1];
        elseif model_name == 3
            model = [ 0, 0, 1; ...
                   0, 1, 1; 0,-1, 1];
        end

        % find max nx nz
        nx = model(:,1);
        nz = model(:,2);
        
        min_nx = min(nx);
        min_nz = min(nz);
        max_nx = max(nx);
        max_nz = max(nz);
        
        % map_kxkz = model(:,1:2);
        ind = 1;
        for i = -max_nx:max_nx
            for j = -max_nz:max_nz
                map_kxkz(ind,1) = i;
                map_kxkz(ind,2) = j;
        
                ind = ind + 1;
            end
        end
        map_kxkz = gp.filter_out_of_range(model(:,1:2), map_kxkz);
        fullmode = map_kxkz;

        % do the mx = nx - kx subtraction
        map_mxmz = zeros([size(map_kxkz),size(map_kxkz,1)]);
        for i = 1:size(map_kxkz,1)
            map_mxmz(:,:,i) = map_kxkz(i,:) - map_kxkz;
        end
        
        % linear and nonlinear operator indexing
        lin_nonlin_ind = struct();
        for i = 1:size(map_kxkz,1)
            lin_nonlin_ind(i).nxnz(1,:) = map_kxkz(i,:);
            
            % filter out wavenumber not in model
            a = gp.filter_out_of_range(map_kxkz, map_mxmz(:,:,i));
            
            lin_nonlin_ind(i).mxmz = a;

            % recover kxkz
            lin_nonlin_ind(i).kxkz = lin_nonlin_ind(i).nxnz(1,:) - ...
                                        lin_nonlin_ind(i).mxmz;
        end

    end
    
    function [lin_nonlin_ind, fullmode] = load_model(nx, nz)

        max_nx = nx;
        min_nx = -nx;

        max_nz = nz;
        min_nz = -nz;

        range_nx = -nx:nx;
        range_nz = -nz:nz;

        % map_kxkz = model(:,1:2);
        ind = 1;
        for i = min_nx:max_nx
            for j = min_nz:max_nz
                fullmode(ind,1) = i;
                fullmode(ind,2) = j;
        
                ind = ind + 1;
            end
        end
        map_kxkz = fullmode;

        map_mxmz = zeros([size(map_kxkz),size(map_kxkz,1)]);
        for i = 1:size(map_kxkz,1)
            map_mxmz(:,:,i) = map_kxkz(i,:) - map_kxkz;
        end
        
        % linear and nonlinear operator indexing
        lin_nonlin_ind = struct();
        for i = 1:size(map_kxkz,1)
            lin_nonlin_ind(i).nxnz(1,:) = map_kxkz(i,:);
            
            % filter out wavenumber not in model
            a = gp.filter_out_of_range(map_kxkz, map_mxmz(:,:,i));
            
            lin_nonlin_ind(i).mxmz = a;

            % recover kxkz
            lin_nonlin_ind(i).kxkz = lin_nonlin_ind(i).nxnz(1,:) - ...
                                        lin_nonlin_ind(i).mxmz;
        end
        
    end
    
    function [size_n_seq, n_seq] = gen_Np_perm(Np)

        % Generate permutation for POD summations

        n_seq = perms(1:Np);

        if Np ~= 1
        tmp = repmat(1:Np,Np,1)';
        n_seq = [n_seq; tmp];
        end

        size_n_seq = size(n_seq,1);
        if Np == 1
            n_seq = [1 1];
            size_n_seq = size(n_seq,1);
        end

    end

    function [L, fullpod, fullmode, fullmode_pod] = prep_L_matrix(nw_pod, fullmode)

        % Construct pod index before going into loops
        L = struct();
        fullmode = repmat(fullmode,nw_pod,1);
        for i = 1:size(fullmode, 1)
            L(i).nxnz = fullmode(i,1:2);
        end

        fullpod = 1:nw_pod;
        fullpod = repmat(fullpod, size(fullmode,1)/nw_pod,1);
        fullpod = reshape(fullpod, [size(fullpod,1)*size(fullpod,2) 1]);

        for i = 1:size(fullpod, 1)
            L(i).n = fullpod(i);
            L(i).n_seq = 1:nw_pod;
            L(i).coeff = 0; % And we also allocate memory for linear coefficients.
        end
        
        % Here we also inflate mode map to include POD
        fullmode_pod = [fullpod fullmode];

    end

    function N = eval_N(ny, Lx, Lz, phi, n, m, k, nx, nz, mx, mz, pod_wave,...
                        diff_weight, int_weight)
        
        kx = nx - mx;
        kz = nz - mz;
    
        In = galerkin_projection.find_wave_number_in_map(nx, nz, pod_wave);
        Im = galerkin_projection.find_wave_number_in_map(mx, mz, pod_wave);
        Ik = galerkin_projection.find_wave_number_in_map(kx, kz, pod_wave);
        
        phinU = phi(1:3:end, n, In);
        phinV = phi(2:3:end, n, In);
        phinW = phi(3:3:end, n, In);
    
        phikU = phi(1:3:end, k, Ik);
        phikV = phi(2:3:end, k, Ik);
        phikW = phi(3:3:end, k, Ik);
    
        phimU = phi(1:3:end, m, Im);
        phimV = phi(2:3:end, m, Im);
        phimW = phi(3:3:end, m, Im);
        
        dphimUdy = math.diff_phi(phimU, diff_weight);
        dphimVdy = math.diff_phi(phimV, diff_weight);
        dphimWdy = math.diff_phi(phimW, diff_weight);
        
        nonlinU = ((2*pi*1i*mx/Lx).*phikU.*phimU + phikV.*dphimUdy + ...
                   (2*pi*1i*mz/Lz).*phikW.*phimU).*conj(phinU);
        nonlinV = ((2*pi*1i*mx/Lx).*phikU.*phimV + phikV.*dphimVdy + ...
                   (2*pi*1i*mz/Lz).*phikW.*phimV).*conj(phinV);
        nonlinW = ((2*pi*1i*mx/Lx).*phikU.*phimW + phikV.*dphimWdy + ...
                   (2*pi*1i*mz/Lz).*phikW.*phimW).*conj(phinW);
        
        nU = 0;
        nV = 0;
        nW = 0;
        for i = 1:ny
            nU = nU + int_weight(i)*nonlinU(i);
            nV = nV + int_weight(i)*nonlinV(i);
            nW = nW + int_weight(i)*nonlinW(i);
        end
    
        N = nU + nV + nW;
        N = -N/sqrt(Lx*Lz);
    end    

    function N = eval_N_k(y, Lx, Lz, phi, n, m, k, nx, nz, kx, kz, mx, mz, pod_wave,...
                    diff_weight, int_weight)

    I_n = gp.find_wave_number_in_map(nx, nz, pod_wave);
    I_k = gp.find_wave_number_in_map(kx, kz, pod_wave);
    I_m = gp.find_wave_number_in_map(mx, mz, pod_wave);

    phiu_n = phi(1:3:end,n, I_n);
    phiv_n = phi(2:3:end,n, I_n);
    phiw_n = phi(3:3:end,n, I_n);

    phiu_k = phi(1:3:end,k, I_k);
    phiv_k = phi(2:3:end,k, I_k);
    phiw_k = phi(3:3:end,k, I_k);    

    phiu_m = phi(1:3:end,m, I_m);
    phiv_m = phi(2:3:end,m, I_m);
    phiw_m = phi(3:3:end,m, I_m);

    dphi_u = math.diff_phi(phiu_m, diff_weight);
    dphi_v = math.diff_phi(phiv_m, diff_weight);
    dphi_w = math.diff_phi(phiw_m, diff_weight);

    int_nonlin_u = ((2*pi*1i*mx/Lx)*phiu_k.*phiu_m + ...
                   phiv_k.*dphi_u + (2*pi*1i*mz/Lz)*phiw_k.*phiu_m).*conj(phiu_n);
    int_nonlin_v = ((2*pi*mx*1i/Lx)*phiu_k.*phiv_m + ...
                   phiv_k.*dphi_v + (2*pi*mz*1i/Lz)*phiw_k.*phiv_m).*conj(phiv_n);
    int_nonlin_w = ((2*pi*mx*1i/Lx)*phiu_k.*phiw_m + ...
                   phiv_k.*dphi_w + (2*pi*mz*1i/Lz)*phiw_k.*phiw_m).*conj(phiw_n);    

    nonlin_u = 0;
    nonlin_v = 0;
    nonlin_w = 0;
    for i = 1:length(y)
        nonlin_u = nonlin_u + int_weight(i)*int_nonlin_u(i);
        nonlin_v = nonlin_v + int_weight(i)*int_nonlin_v(i);
        nonlin_w = nonlin_w + int_weight(i)*int_nonlin_w(i);    
    end
    
    nonlin = nonlin_u + nonlin_v + nonlin_w;

    N = -nonlin/sqrt(Lx*Lz);

    end

    function L = eval_L(Re, ny, y, Lx, Lz, phi, n, m, nx, nz, pod_wave, w)

        % phi: wall normal mode from modal reduction     
        % ny: number of y division from cfd
        % nx, nz: stream/span-wise wave number
        % w: integration weights
    
        proj_wave_i = galerkin_projection.find_wave_number_in_map(nx, nz, pod_wave);
        
        invRe = 1/Re;
        
        kron = 0;
        if n == m
            kron = 1;
        end
        
        
        term1 = invRe*((2*pi*nx/Lx)^2 + (2*pi*nz/Lz)^2)*kron;
        
        [~, diff_weight] = math.f_chebdiff(ny, 1);
        
        dphiuMdy = math.diff_phi(phi(1:3:end,m,proj_wave_i), diff_weight);
        dphivMdy = math.diff_phi(phi(2:3:end,m,proj_wave_i), diff_weight);
        dphiwMdy = math.diff_phi(phi(3:3:end,m,proj_wave_i), diff_weight);
        
        
        dphiuMdy_innerprod = conj(dphiuMdy).*dphiuMdy;
        dphivMdy_innerprod = conj(dphivMdy).*dphivMdy;
        dphiwMdy_innerprod = conj(dphiwMdy).*dphiwMdy;
        
        dphiMdy_integral = 0;
        for i = 1:ny
            dphiMdy_integral = dphiMdy_integral + w(i)*dphiuMdy_innerprod(i) + ...
                                                  w(i)*dphivMdy_innerprod(i) + ...
                                                  w(i)*dphiwMdy_innerprod(i);
        end
        dphiMdy_integral = invRe*dphiMdy_integral;
        
        phiuMN = y.*conj(phi(1:3:end,m,proj_wave_i)).*(phi(1:3:end,n,proj_wave_i));
        phivMN = y.*conj(phi(2:3:end,m,proj_wave_i)).*(phi(2:3:end,n,proj_wave_i));
        phiwMN = y.*conj(phi(3:3:end,m,proj_wave_i)).*(phi(3:3:end,n,proj_wave_i));
        
        phiMN_integral = 0;
        for i = 1:ny
            phiMN_integral = phiMN_integral + w(i)*phiuMN(i) + w(i)*phivMN(i) ...
                                            + w(i)*phiwMN(i);
        end
        phiMN_integral = (2*pi*1i*nx/Lx)*phiMN_integral;
        
        phiuvMN_integrand = conj(phi(2:3:end,m,proj_wave_i)).*phi(1:3:end,n,proj_wave_i);
        
        phiuvMN_integral = 0;
        for i = 1:ny
            phiuvMN_integral = phiuvMN_integral + w(i)*phiuvMN_integrand(i);
        end
        
        L = - term1 - dphiMdy_integral - phiMN_integral - phiuvMN_integral;
     end
    
    function L = eval_L_k(Re, ny, y, Lx, Lz, phi, n, k, nx, nz, pod_wave, dw, w)
        invRe = 1/Re;
        
        term1 = -invRe*((2*pi*nx/Lx)^2 + (2*pi*nz/Lz)^2);
        
        I2n = gp.find_wave_number_in_map(nx, nz, pod_wave);
        I1k = I2n;
        
        phi2n = phi(2:3:end,n,I2n);
        phi1k_conj = conj(phi(1:3:end,k,I1k));
        
        term2_int = phi2n.*phi1k_conj;
        
        term2 = 0;
        for i = 1:ny
            term2 = term2 + term2_int(i)*w(i);
        end
        term2 = -term2;
        
        term3u_int = y.*phi(1:3:end,n, I2n).*conj(phi(1:3:end,k,I2n));
        term3v_int = y.*phi(2:3:end,n, I2n).*conj(phi(2:3:end,k,I2n));
        term3w_int = y.*phi(3:3:end,n, I2n).*conj(phi(3:3:end,k,I2n));
        
        term3u = 0;
        term3v = 0;
        term3w = 0;
        for i = 1:ny
            term3u = term3u + term3u_int(i)*w(i);
            term3v = term3v + term3v_int(i)*w(i);
            term3w = term3w + term3w_int(i)*w(i);
        end
        
        term3 = -2*pi*1i*nx*(term3u + term3w + term3v)/Lx;
        
        term4_u_n = phi(1:3:end,n,I2n);
        term4_v_n = phi(2:3:end,n,I2n);
        term4_w_n = phi(3:3:end,n,I2n);
        
        term4_u_k = phi(1:3:end,k,I2n);
        term4_v_k = phi(2:3:end,k,I2n);
        term4_w_k = phi(3:3:end,k,I2n);
        
        term4_u_diff = math.diff_phi(term4_u_n,dw);
        term4_v_diff = math.diff_phi(term4_v_n,dw);
        term4_w_diff = math.diff_phi(term4_w_n,dw);
        
        term4_u_k_diff = conj(math.diff_phi(term4_u_k,dw));
        term4_v_k_diff = conj(math.diff_phi(term4_v_k,dw));
        term4_w_k_diff = conj(math.diff_phi(term4_w_k,dw));
        
        term4_u_int = term4_u_diff.*term4_u_k_diff;
        term4_v_int = term4_v_diff.*term4_v_k_diff;
        term4_w_int = term4_w_diff.*term4_w_k_diff;
        
        term4_u = 0;
        term4_v = 0;
        term4_w = 0;
        for i = 1:ny
            term4_u = term4_u + w(i)*term4_u_int(i);
            term4_v = term4_v + w(i)*term4_v_int(i);
            term4_w = term4_w + w(i)*term4_w_int(i);
        end
        
        term4 = term4_u + term4_v + term4_w;
        term4 = -invRe*term4;
        
        L = term1 + term2 + term3 + term4;
    end
    
    function [L, N, fullmode_pod] = get_L_N_struct(Re, Np, phi, lin_nonlin_ind, fullmode, pod_wave)

        % This function encasuplate all processes to compute linear and
        % nonlinear coefficients

        [x,y,z,nx,ny,nz,Lx,Lz] = modal_decomposition.read_geom;
        
        const_lin_nonlin_ind = lin_nonlin_ind; % This line saves the index map to be used in ODE45 (nothing mathematical)
        
        w = math.f_chebyshev_int_weight(ny);
        [~,dw] = math.f_chebdiff(ny, 1);
                        
        [L, fullpod, fullmode, fullmode_pod] = gp.prep_L_matrix(Np, fullmode);
        L = gp.find_POD_pair_conj_index_in_model(L, Np, const_lin_nonlin_ind,fullmode);
        
        for i = 1:size(L,2)
            size_n_seq = size(L(i).n_seq,2);
            nxnz = L(i).nxnz;
            for i_pod = 1:size_n_seq
                L(i).coeff(i_pod) = gp.eval_L_k(Re, ny, y, Lx, Lz, phi, ...
                    L(i).n, L(i).n_seq(i_pod), nxnz(1), nxnz(2), pod_wave, dw, w);
            end
        end
        
        N = struct();
        for i = 1:size(fullpod, 1)
            N(i).n = fullpod(i); 
            N(i).coeff = 0;
        end
        
        % Inflate lin_nonlin_ind to match POD number of POD expansions
        totl_mode = size(lin_nonlin_ind,2);
        for i = 1:totl_mode*(Np-1)
            lin_nonlin_ind(i + totl_mode).nxnz = lin_nonlin_ind(i).nxnz;
            lin_nonlin_ind(i + totl_mode).mxmz = lin_nonlin_ind(i).mxmz;
            lin_nonlin_ind(i + totl_mode).kxkz = lin_nonlin_ind(i).kxkz;
        end
        
        [size_n_seq, n_seq] = gp.gen_Np_perm(Np);
        
        % Assign wavenumbers to struct. And search for mxmz, kxkz index in a
        for i = 1:size(N,2)
            N(i).n_seq = n_seq;
            N(i).nxnz = lin_nonlin_ind(i).nxnz;
            N(i).mxmz = lin_nonlin_ind(i).mxmz;
            N(i).kxkz = lin_nonlin_ind(i).kxkz;
        
            % search for corresponding mxmz, kxkz index in a (include Np)
            for i_mxmz = 1:size(N(i).mxmz,1)
        
                mxmz = N(i).mxmz(i_mxmz,:);
        
                for i_Np = 1:size_n_seq
                    I = gp.find_wave_number_pod_in_map(N(i).n_seq(i_Np,1), mxmz(1), mxmz(2), fullmode_pod);
                    N(i).a_mxmz_loc(i_mxmz,i_Np) = I; 
                end
        
            end
        
            for i_kxkz = 1:size(N(i).kxkz,1)
        
                kxkz = N(i).kxkz(i_kxkz,:);
        
                for i_Np = 1:size_n_seq 
                    I = gp.find_wave_number_pod_in_map(N(i).n_seq(i_Np,2), kxkz(1), kxkz(2), fullmode_pod);
                    N(i).a_kxkz_loc(i_kxkz,i_Np) = I; 
                end            
        
            end
        end
        
        % Fills in Nonlinear coefficients
        for i = 1:size(N,2)
            
            Np_c = N(i).n;
            nxnz = N(i).nxnz;
            kxkz_arr = N(i).kxkz;
            mxmz_arr = N(i).mxmz;
        
            N(i).coeff = zeros(size(kxkz_arr,1), size_n_seq);
        
            for i_mxmz = 1:size(mxmz_arr,1)
                kx = kxkz_arr(i_mxmz,1);
                kz = kxkz_arr(i_mxmz,2);
        
                mx = mxmz_arr(i_mxmz,1);
                mz = mxmz_arr(i_mxmz,2);
        
                for i_local_pod = 1:size_n_seq
        
                    m = N(i).n_seq(i_local_pod,1);
                    k = N(i).n_seq(i_local_pod,2);
                    N(i).coeff(i_mxmz, i_local_pod) = gp.eval_N_k(y, Lx, Lz, phi, Np_c, m, k, ...
                                                        nxnz(1),nxnz(2), kx, kz, mx, mz,...
                                                        pod_wave, dw, w);
                end
              
            end
        
        end        

    end

    function adot = eval_adot(t, a, L, N, nw_pod, a0)
        
        disp(['t = ', num2str(t)]);

        a_len = size(a0,1);

        adot_lin = zeros(a_len,1);
        adot_nonlin = adot_lin;
        
        % linear coefficient
        for i_lin = 1:a_len

            for i_pod = 1:nw_pod

                anxnz = a(L(i_lin).pod_pair(i_pod));

                adot_lin(i_lin) = adot_lin(i_lin) + L(i_lin).coeff(i_pod)*anxnz;
            end

        end
        
        % total number of permutations of the PODs
        perms_length = size(N(1).n_seq,1);

        % nonlinear coefficient
        for i_nonlin = 1:a_len
            
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
            
        end

        adot = adot_lin + adot_nonlin;

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

    % please stop. I have ran out of variable name......
end