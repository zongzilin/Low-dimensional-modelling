classdef galerkin_projection
    methods(Static)
    
    function I = find_wave_number_in_map(nx, nz, pod_wave)
    
        r = size(pod_wave,1);
        nx = abs(nx);
        nz = abs(nz);
        
        for i = 1:r
        
            if ((pod_wave(i,1) == nx) && (pod_wave(i,2) == nz))
                break
            end
    
        end
    
        I = i;
    
    end

    function out = filter_out_of_range(model, wave_map)
        % this function filters out the mx = nx - kx and mz that are no in
        % the model wavenumber range
    
        [r_wave, ~] = size(wave_map);
    
        for i = 1:r_wave
            a(i) = ismember(abs(wave_map(i,:)),model,'row');
        end
    
        wave_map = wave_map(a,:);
        out = wave_map;
        
    end    
    
    function [lin_nonlin_ind, model] = load_predefined_model(model_name)

        if model_name == 6
            model = [ 0, 0, 1; ...
                   0, 1, 1; 0,-1, 1; ...
                   0, 2, 1; 0,-2, 1; ...
                   1, 0, 1;-1, 0, 1; ...
                   1, 1, 1;-1,-1, 1; ...
                   1,-1, 1;-1, 1, 1];
        elseif model_name == 9
                model = [ 0, 0, 1; 0, 1, 1; ...
                   1, 0, 1; 1, 1, 1; ...
                   1,-1, 1; 0, 2, 1; ...
                   1, 2, 1; 1,-2, 1; ...
                   0, 3, 1];
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
        map_kxkz = galerkin_projection.filter_out_of_range(model(:,1:2), map_kxkz);
        
        % do the mx = nx - kx subtraction
        map_mxmz = zeros([size(map_kxkz),size(map_kxkz,1)]);
        for i = 1:size(map_kxkz,1)
            map_mxmz(:,:,i) = map_kxkz(i,:) - map_kxkz;
        end
        
        % linear and nonlinear operator indexing
        lin_nonlin_ind = struct();
        for i = 1:size(map_kxkz,1)
            lin_nonlin_ind(i).kxkz(1,:) = map_kxkz(i,:);
            
            % filter out wavenumber not in model
            a = galerkin_projection.filter_out_of_range(map_kxkz, map_mxmz(:,:,i));
        
            lin_nonlin_ind(i).mxmz = a;
        end

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
    
    function L = eval_L(Re, ny, y, Lx, Lz, phi, n, m, nx, nz, pod_wave, w)
             
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

    end
end