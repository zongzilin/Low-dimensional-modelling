classdef galerkin_projection
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

    function N = eval_N(ny, Lx, Lz, phi, n, m, k, nx, nz, mx, mz, pod_wave,...
                        diff_weight, int_weight)
        
        kx = nx - mx;
        kz = nz - mz;
    
        In = find_wave_number_in_map(nx, nz, pod_wave);
        Im = find_wave_number_in_map(mx, mz, pod_wave);
        Ik = find_wave_number_in_map(kx, kz, pod_wave);
        
        phinU = phi(1:3:end, n, In);
        phinV = phi(2:3:end, n, In);
        phinW = phi(3:3:end, n, In);
    
        phikU = phi(1:3:end, k, Ik);
        phikV = phi(2:3:end, k, Ik);
        phikW = phi(3:3:end, k, Ik);
    
        phimU = phi(1:3:end, m, Im);
        phimV = phi(2:3:end, m, Im);
        phimW = phi(3:3:end, m, Im);
        
        dphimUdy = modal_decomposition.diff_phi(phimU, diff_weight);
        dphimVdy = modal_decomposition.diff_phi(phimV, diff_weight);
        dphimWdy = modal_decomposition.diff_phi(phimW, diff_weight);
        
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
        
        dphiuMdy = math.diff_phi(phi(1:3:end,1,proj_wave_i), diff_weight);
        dphivMdy = math.diff_phi(phi(2:3:end,1,proj_wave_i), diff_weight);
        dphiwMdy = math.diff_phi(phi(3:3:end,1,proj_wave_i), diff_weight);
        
        
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
        
        phiuMN = y.*conj(phi(1:3:end,1,proj_wave_i)).*(phi(1:3:end,1,proj_wave_i));
        phivMN = y.*conj(phi(2:3:end,1,proj_wave_i)).*(phi(2:3:end,1,proj_wave_i));
        phiwMN = y.*conj(phi(3:3:end,1,proj_wave_i)).*(phi(3:3:end,1,proj_wave_i));
        
        phiMN_integral = 0;
        for i = 1:ny
            phiMN_integral = phiMN_integral + w(i)*phiuMN(i) + w(i)*phivMN(i) ...
                                            + w(i)*phiwMN(i);
        end
        phiMN_integral = (2*pi*1i*nx/Lx)*phiMN_integral;
        
        phiuvMN_integrand = conj(phi(2:3:end,1,proj_wave_i)).*phi(1:3:end,1,proj_wave_i);
        
        phiuvMN_integral = 0;
        for i = 1:ny
            phiuvMN_integral = phiuvMN_integral + w(i)*phiuvMN_integrand(i);
        end
        
        L = - term1 - dphiMdy_integral - phiMN_integral - phiuvMN_integral;
     end

    end
end