classdef sindy
    methods(Static)


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


    end
end