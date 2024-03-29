classdef vis
    methods(Static)

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



    end
end