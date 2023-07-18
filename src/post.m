classdef post
    methods(Static)

        function [u, v, w, x, y, z] = fluctuation_rms(data_set_string)

            % FROM Y. HWANG
        
            load(data_set_string);

            [x,y,z,~,ylen,~,Lx,Lz] = modal_decomposition.read_geom;            
        
            ind = and(fullmode(:,1) == 0, fullmode(:,2) == 0);
            a(:,ind) = [];
            
            fullmode_pod_no_mean = fullmode_pod(~ind,:);
            
            reduced_dof = size(fullmode_pod_no_mean, 1);
            
            u_fluc = zeros(ylen, 1);
            v_fluc = u_fluc;
            w_fluc = u_fluc;
            
            for i = 1:reduced_dof
            
                nxnz = fullmode_pod_no_mean(i, 2:3);
                Np = fullmode_pod_no_mean(i, 1);
            
                I = gp.find_wave_number_in_map(nxnz(1), nxnz(2), pod_wave);
            
                u_fluc = u_fluc + mean(real(a(:,i)).^2 + imag(a(:,i)).^2).*phi(1:3:end,Np,I).*conj(phi(1:3:end,Np,I))/(Lx*Lz);    
                v_fluc = v_fluc + mean(real(a(:,i)).^2 + imag(a(:,i)).^2).*phi(2:3:end,Np,I).*conj(phi(2:3:end,Np,I))/(Lx*Lz);   
                w_fluc = w_fluc + mean(real(a(:,i)).^2 + imag(a(:,i)).^2).*phi(3:3:end,Np,I).*conj(phi(3:3:end,Np,I))/(Lx*Lz);   

            end
            
            u = real(abs(u_fluc));      
            v = real(abs(v_fluc));
            w = real(abs(w_fluc));
        end        
        
        function [M,time] = modal_energy_from_file(data_set_string)

            load(data_set_string);
            
            [x,y,z,nx,ny,nz,Lx,Lz] = modal_decomposition.read_geom;
        
            % get size of a_gp
            [time, mode_count] = size(a_gp);
            
            % modal energy 
            M = zeros(time, mode_count/Np);
            
            n_seq_len = size(L(1).n_seq,2);
            
            
            for i_t = 1:time
            
                for i_pod = 1:mode_count/Np
            
                    for i_n_seq = 1:n_seq_len
                        M(i_t,i_pod) = M(i_t, i_pod) + abs(a_gp(i_t,L(i_pod).pod_pair(i_n_seq)));
                    end
            
                end
            
            
            end
            
            M = M/sqrt(Lx*Lz);

            time = 1:time;
    
            end
        
        function M = modal_energy(a, Np, L, Lx, Lz)
            
            % get size of a_gp
            [time, mode_count] = size(a);
            
            % modal energy 
            M = zeros(time, mode_count/Np);
            
            n_seq_len = size(L(1).n_seq,2);
            
            
            for i_t = 1:time
            
                for i_pod = 1:mode_count/Np
            
                    for i_n_seq = 1:n_seq_len
                        M(i_t,i_pod) = M(i_t, i_pod) + abs(a(i_t,L(i_pod).pod_pair(i_n_seq)));
                    end
            
                end
            
            
            end
            
            M = M/sqrt(Lx*Lz);


        end

        function [u, v, w, x, y, z] = flow_reconstruction(data_set_string)

            load(data_set_string);

            [x,y,z,nx,ny,nz,Lx,Lz] = modal_decomposition.read_geom;

            sim_time = size(a_gp,1);

            x_len = length(x);
            y_len = length(y);
            z_len = length(z);

            a = a_gp;

            % Search pod mode in full mode            
            mode_count = size(fullmode, 1);
            pod_mode_map = zeros(mode_count,1);
            i_ind = 1;
            for i_mc = 1:mode_count
                
                if (fullmode(i_mc,1) || fullmode(i_mc,2) ~= 0)
                    pod_mode_map(i_mc) = gp.find_wave_number_in_map(fullmode(i_mc,1), ...
                        fullmode(i_mc,2), pod_wave);
                else
                    mean_flow_index(i_ind) = i_mc;
                    i_ind  = i_ind + 1;
                end
            
            end

            % re-arrange phi indexes
            pod_mode_map_reconstruct = fullmode_pod;
            pod_mode_map_reconstruct(mean_flow_index,:) = [];
            pod_mode_map = nonzeros(pod_mode_map);
            phi = phi(:,:,pod_mode_map);
            
            u_fluc = zeros(y_len,1);
            v_fluc = u_fluc;
            w_fluc = u_fluc;
            
            for i = 1:size(phi,3)
                
                i_pod = pod_mode_map_reconstruct(i,1);
                
                u_fluc = u_fluc + mean(real(a(:,i)).^2 + imag(a(:,i)).^2)*phi(1:3:end,i_pod,i).*conj(phi(1:3:end,i_pod,i))/Lx/Lz;
                v_fluc = v_fluc + mean(real(a(:,i)).^2 + imag(a(:,i)).^2)*phi(2:3:end,i_pod,i).*conj(phi(2:3:end,i_pod,i))/Lx/Lz;
                w_fluc = w_fluc + mean(real(a(:,i)).^2 + imag(a(:,i)).^2)*phi(3:3:end,i_pod,i).*conj(phi(3:3:end,i_pod,i))/Lx/Lz;
            end
            
            u = real(u_fluc);
            v = real(v_fluc);
            w = real(w_fluc);
        end
        
        function c = xcorrelation(t, M_nxnz, M_mxmz, tau)
    
        total_time = size(t,1);
        
        M_nxnz_mean = mean(M_nxnz);
        M_mxmz_mean = mean(M_mxmz);
        
        M_tilde_nxnz = M_nxnz - M_nxnz_mean;
        M_tilde_mxmz = M_mxmz - M_mxmz_mean;
        
        % M_tilde_nxnz_rms = sqrt(mean(M_tilde_nxnz(tau + 1:total_time - tau).^2));
        % M_tilde_mxmz_rms = sqrt(mean(M_tilde_mxmz(tau + 1:total_time - tau).^2));
        if tau >= 0
            M_tilde_nxnz_rms = sqrt(mean(M_tilde_nxnz(1:total_time-tau).^2));
            M_tilde_mxmz_rms = sqrt(mean(M_tilde_mxmz(1:total_time-tau).^2));
        elseif tau < 0
            M_tilde_nxnz_rms = sqrt(mean(M_tilde_nxnz(1:total_time+tau).^2));
            M_tilde_mxmz_rms = sqrt(mean(M_tilde_mxmz(1:total_time+tau).^2));
        end
        
        i = 1;
        if tau >= 0
            for i_t = 1:total_time-tau
                c_temp(i) = M_tilde_nxnz(i_t + tau)*M_tilde_mxmz(i_t);
                i = i + 1;
            end
        elseif tau < 0
            for i_t = abs(tau)+1:total_time
                c_temp(i) = M_tilde_nxnz(i_t+tau)*M_tilde_mxmz(i_t);
                i = i + 1;
            end
        end
        
        c = mean(c_temp)/M_tilde_nxnz_rms/M_tilde_mxmz_rms;
    
    end

    end
end