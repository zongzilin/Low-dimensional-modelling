classdef gp
    
    % THIS FILE CONTAINS ALL THE RELEVANT FUNCTION TO COMPUTE GALERKIN
    % PROJECTION



    methods(Static)
    
    function I = find_wave_number_in_map(nx, nz, pod_wave)

        % This function finds nx (streamwise) and nz (spanwise) wavenumber
        % from a map matrix of size n-by-2. Column 1: all nx
        % Column 2: all nz
    
        r = size(pod_wave,1);
        
        for i = 1:r
        
            if ((pod_wave(i,1) == nx) && (pod_wave(i,2) == nz))
                break
            end
    
        end
    
        I = i;
    
    end

    function I = find_wave_number_pod_in_map(np, nx, nz, pod_wave)
        % This function finds nx (streamwise) and nz (spanwise) wavenumber
        % from a map matrix of size n-by-3       
        % pod_wave: column 1: all pod numbers
        %           column 2&3: all nx and nz wavenumbers
    
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
        % This function acts to find the the index to force conjugation in
        % model. 
        % Input: L                    ---- struct of initial L  
        %        nw_pod               ---- number of pod for the current model
        %        const_lin_nonlin_ind ---- this struct contains all pairs of linear
        %                                  and nonlinear wavenumber interactions.
        %                                  The prefix " const " signifies
        %                                  that this struct will not be
        %                                  modified throughout the code
        %        fullmode             ---- this is a n-by-2 matrix where
        %                                  first column: all nx wavenumber
        %                                  second column: all nz wavenumber
        % Output: L                   ---- the same input L but the indexes
        %                                  to force conjugation is written
        %                                  at the last field as
        %                                  "conj_ind_start" and
        %                                  "conj_ind_end"

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

        % This function filters out the nonlinear interaction of mx = nx - kx and mz that 
        % are not in the model wavenumber range
        % Inputs: model    ---- An n-by-2 matrix with first column: nx
        %                       second column: nz
        %         wave_map ---- A n-by-2 matrix from the output of function
        %                       load_model()
        % Outputs: out ---- A 2 column matrix of a filtered wavenumber map.
        %                   first column: nx
        %                   second column: nz
        %          loc ---- The indexes of filtered wavenumbers
    
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
        
        % This function hard codes the ROM from Smith(2005). This function
        % also contains a mini version of function filter_out_of_range()
        %
        % Inputs: model_name ---- A number that specify the model name as
        %                         in Smith(2005)
        % Outputs: lin_nonlin_ind ---- A full struct with all the possible
        %                              nx nz, kx kz, mx mz interaction
        %          fullmode       ---- All the wavenumber specified by
        %                              Smith(2005)
        %                              Column 1&2: nx&nz
        %          model          ---- All the wavenumber and pod specified
        %                              by Smith(2005)
        %                              Column 1&2: nx&nz
        %                              Column 3  : pods
        %
        % ============================ NOT USED ===========================
        % ============================ NOT USED ===========================

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
        elseif model_name == 1
            model = [1, 1, 1];
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

        % This function loads a generic model with nx from -nx to +nx and
        % -nz to +nz
        % Inputs: nx ---- streamwise wavenumber range
        %         nz ---- spanwise wavenumber range
        % Outputs: lin_nonlin_ind ---- As specified in
        %                              load_predefined_model()
        %          fullmode       ---- As specified in
        %                              load_predefined_model()

        switch nargin
            case 2
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
            case 1
                [lin_nonlin_ind, fullmode, ~] = gp.load_predefined_model(nx);                
        end
        
    end
    
    function [size_n_seq, n_seq] = gen_Np_perm(Np)

        % Generate permutation for POD summations
        % Inputs:  Np ---- No. of PODs in the current model
        % Outputs: n_seq      ---- a sequence of all possible 1:Np permutations
        %          size_n_seq ---- The number of permutations (i.e
        %                          length(n_seq))


        np_arr = 1:Np;
        np_mat = repmat(np_arr, Np);
        
        np_arr_1 = reshape(np_mat,Np^3,1); % ????
        np_arr_2 = reshape(np_mat',Np^3,1); % ????
        
        n_seq = [np_arr_1 np_arr_2];
        n_seq = n_seq(1:Np^2,:);

        size_n_seq = size(n_seq,1);

        % this function is a MESS
    end

    function [L, fullpod, fullmode, fullmode_pod] = prep_L_matrix(nw_pod, fullmode)
        
        % This function generates the complete linear coefficient matrix
        % before finally going into the ode timestepping
        % Inputs: nw_pod   ---- No. of PODs in the current model
        %         fullmode ---- An 2 column matrix with 
        %                       first column:  nx
        %                       second coulmn: nz
        % Outputs: L            ---- A complete linear coefficient struct with all required info to ode
        %                            timestep
        %          fullmode     ---- As defined above
        %          fullmode_pod ---- A 3 column matrix
        %                            Column 1&2: nx and nz
        %                            Column 3: all pods numbers

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

    function N = eval_N_k(y, Lx, Lz, phi, n, m, k, nx, nz, kx, kz, mx, mz, pod_wave,...
                    diff_weight, int_weight)

    % This funcion evaluates the nonlinear coefficient of the Galerkin
    % projection
    % Inputs: y (double array) --- An array of wall normal axis on
    %                              chebeshev points
    %         Lx (double)      --- Length of domain in streamwise direction
    %         Lz (double)      --- Length of domain in spanwise direction 
    %         phi (complex double array) --- Array of basis function
    %         n (integer)      --- nonlinear wavenumber coupling 
    %         m (integer)      --- nonlinear wavenumber coupling
    %         k (integer)      --- nonlinear wavenumber coupling
    %         nx, nz, kz, kz, mx, mz (integer) --- nonlinear wavenumber coupling
    %         pod_wave (integer array) --- Array of wavenumber map 
    %         diff_weight (double array) --- Differential weight of each
    %                                        chebeyshev points in y
    %         int_weight (double array) --- Integration weight of each
    %                                       chebeyshev points in y
    % Outputs: N (complex double) --- nonlinear coefficient of current
    %                                 nonlninear wavenumber coupling

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

    dphi_u = diff_weight(:,:,1)*phiu_m;
    dphi_v = diff_weight(:,:,1)*phiv_m;
    dphi_w = diff_weight(:,:,1)*phiw_m;

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
        
        %%%%%%%%%%%%%%%%%%%%%% DEPRECATED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    
    function [L, lam, diff, production] = eval_L_k(Re, ny, y, Lx, Lz, phi, n, k, nx, nz, pod_wave, dw, w)
    
        % This function evaluates the linear coefficients in Galerkin
        % projection. NOT USED in the final code


        invRe = 1/Re;
        
        kron = 0;
        if n == k
            kron = 1;
        end
        term1 = -invRe*((2*pi*nx/Lx)^2 + (2*pi*nz/Lz)^2)*kron;
        
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
        
        term3 = 2*pi*1i*nx*(term3u + term3w + term3v)/Lx;
        
        term4_u_n = phi(1:3:end,n,I2n);
        term4_v_n = phi(2:3:end,n,I2n);
        term4_w_n = phi(3:3:end,n,I2n);
        
        term4_u_k = phi(1:3:end,k,I2n);
        term4_v_k = phi(2:3:end,k,I2n);
        term4_w_k = phi(3:3:end,k,I2n);
        
        term4_u_diff = math.diff_phi(term4_u_n,dw(:,:,1));
        term4_v_diff = math.diff_phi(term4_v_n,dw(:,:,1));
        term4_w_diff = math.diff_phi(term4_w_n,dw(:,:,1));
        
        term4_u_k_diff = math.diff_phi(conj(term4_u_k),dw(:,:,1));
        term4_v_k_diff = math.diff_phi(conj(term4_v_k),dw(:,:,1));
        term4_w_k_diff = math.diff_phi(conj(term4_w_k),dw(:,:,1));
        
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
        
        diff = term1 + term4;
        production = term2;
        lam = term3;

        L = term1 + term2 + term3 + term4;
        
    end

    function [L, lam, diff, prod] = eval_L_m(Re, ny, y, Lx, Lz, phi, n, m, nx, nz, pod_wave, dw, w)

        % This function evaluates the linear coefficient in Galerkin
        % projection. 
        % Inputs: Re (double) --- Reynolds number
        %         ny (int) --- Number of grid points in wall normal direction
        %         y  (double array) --- An array of wall normal coordinates. len(y) = ny
        %         Lx (double) --- Length of streamwise axis 
        %         Lz (double) --- Length of spanwise axis
        %         phi (complex double array) --- wall normal POD basis function
        %         n  (integer) --- wave number coupling in linear Galerkin projection
        %         m  (integer) --- same as above
        %         pod_wave (complex double array) --- column 1: all pod numbers
        %                                             column 2&3: all nx and nz wavenumbers 
        %         dw (double array) --- differentiate weight for Cheybeshev points
        %         w (double array) --- integration weights for Cheybeshev
        %                              points
        % Outputs: L (complex double array) --- full linear matrix for
        %                                       Galerkin projection
        %          lam (complex double array) --- Linear matrix for
        %                                         galerkin projection of the laminar state
        %          diff (complex double array) --- Linear matrix for
        %                                          Galerkin projection of diffusion (i.e \nabla^2 u) 
        %          prod (complex double array) --- Linear matrix for
        %                                          Galerkin projection of production term


        invRe = 1/Re;
        
        lam = 0;
        diff = 0;
        prod = 0;

        index = gp.find_wave_number_in_map(nx, nz, pod_wave);

        phiu_n = phi(1:3:end, n, index);
        phiv_n = phi(2:3:end, n, index);
        phiw_n = phi(3:3:end, n, index);

        phiu_m = phi(1:3:end, m, index);
        phiv_m = phi(2:3:end, m, index);
        phiw_m = phi(3:3:end, m, index);

        prod_int = conj(phiv_m).*phiu_n;

        for i = 1:ny
            prod = prod + w(i)*prod_int(i);
        end
        prod = -prod;

        % laminar state !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        lam_u_int = y.*conj(phiu_m).*phiu_n;
        lam_v_int = y.*conj(phiv_m).*phiv_n;
        lam_w_int = y.*conj(phiw_m).*phiw_n;
        
        lam_u = 0;
        lam_v = 0;
        lam_w = 0;

        for i = 1:ny
            lam_u = lam_u + w(i)*lam_u_int(i);
            lam_v = lam_v + w(i)*lam_v_int(i);
            lam_w = lam_w + w(i)*lam_w_int(i);
        end
        % !!!!!!!!!!!!!!!!!!! (I MISSED A NEGATIVE SIGN HERE BUT IT
        %                       WORKS)
        lam = 2*pi*1i*nx*(lam_u + lam_v + lam_w)/Lx;
       
        
        % diffusion/viscous
        kron = 0;
        if n == m
            kron = 1;
        end
        diff_1 = -invRe*((2*pi*nx/Lx)^2 + (2*pi*nz/Lz)^2)*kron;

        diff_u = 0;
        diff_v = 0;
        diff_w = 0;

        diff_u_int = dw(:,:,1)*conj(phiu_m).*(dw(:,:,1)*phiu_n);
        diff_v_int = dw(:,:,1)*conj(phiv_m).*(dw(:,:,1)*phiv_n);
        diff_w_int = dw(:,:,1)*conj(phiw_m).*(dw(:,:,1)*phiw_n);

        for i = 1:ny
            diff_u = diff_u + w(i)*diff_u_int(i);
            diff_v = diff_v + w(i)*diff_v_int(i);
            diff_w = diff_w + w(i)*diff_w_int(i);
        end
        diff = -invRe*(diff_u + diff_v + diff_w) + diff_1;
        
        % Combine all linear matrix components
        L = lam + diff + prod;
        
    end
    
    function [L, N, fullmode_pod] = get_L_N_struct(Re, Np, phi, lin_nonlin_ind, fullmode, pod_wave)

        % This function encasuplate all processes to compute linear and
        % nonlinear coefficients
        % Inputs:   Re (double) --- Reynolds number
        %           Np (integer) --- Number of POD modes for each wavenumber
        %           phi (complex double array) --- Basis function array
        %           lin_nonlin_ind (integer array) --- Map of each
        %                                              wavenumber coupling
        %           fullmode (integer array) --- Array containing all POD
        %                                        modes and wavenumber
        %                                        indexes.
        %                                        (e.g: 1st col: POD modes,
        %                                        2nd col: wavenumber in x
        %                                        3rd col: wavenumber in z)
        %           pod_wave (integer array) --- An array of wavenumber map
        % Outputs: L (complex double array) --- An filled completed Linear
        %                                       matrix that is ready to go
        %                                       into timestepping
        %          N (complex double array) --- An filled completed
        %                                       nonlinear matrix that is
        %                                       ready to go into
        %                                       timestepping
        %          fullmode_pod (integer array) --- An finalised version of
        %                                           of POD modes to
        %                                           wavenumber map

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
                [L(i).coeff(i_pod), L(i).lam(i_pod), L(i).diff(i_pod), L(i).prod(i_pod)] = ...
                    gp.eval_L_m(Re, ny, y, Lx, Lz, phi, ...
                    L(i).n, L(i).n_seq(i_pod), nxnz(1), nxnz(2), pod_wave, dw, w);
            end
            disp(['linear coeff ',num2str(i),'/',num2str(size(L,2))])
        end

        % rearrange field (just for aesthetic)
        L = orderfields(L, [1:5,9,6:8,10:11]);
        L = orderfields(L, [1:6,10,7:9,11]);
        L = orderfields(L, [1:7,11,8:10]);

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
        
        % Assign wavenumbers to struct. And search for mxmz, kxkz index in
        % "a"
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

            if mod(i,100) == 0
            disp(['prep nonlinear coeff ',num2str(i),'/',num2str(size(N,2))])
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

             disp(['nonlinear coeff ',num2str(i),'/',num2str(size(N,2))])
        
        end        

    end
    
    function [L, fullmode_pod] = get_L_struct(Re, Np, phi, lin_nonlin_ind, fullmode, pod_wave)

    % This function encasuplate all processes to compute linear and
    % nonlinear coefficients
    % Inputs: SEE get_L_N_struct
    % Ouputs: SEE get_L_N_struct

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
                [L(i).coeff(i_pod), L(i).lam(i_pod), L(i).diff(i_pod), L(i).prod(i_pod)] = ...
                    gp.eval_L_m(Re, ny, y, Lx, Lz, phi, ...
                    L(i).n, L(i).n_seq(i_pod), nxnz(1), nxnz(2), pod_wave, dw, w);
            end
        end
    
        % rearrange field (just for aesthetic)
        L = orderfields(L, [1:5,9,6:8,10:11]);
        L = orderfields(L, [1:6,10,7:9,11]);
        L = orderfields(L, [1:7,11,8:10]);
    end
    
    function a_initial = set_initial_conditions(t_initial, a_dns, fullmode_pod, pod_wave)

    % This function sets the inital condition for ode timestepping. Initial
    % conditions are obtained from DNS data
    % Inputs: t_initial (integer) --- Timestep from DNS to set as initial
    %                                 conditions
    %         a_dns (complex double array) --- Modal time coefficients from
    %                                          DNS
    %         fullmode_pod (integer array) --- SEE get_L_N_struct
    %         pod_wave (integer array) --- SEE get_L_N_struct
    % Outputs: a_initial (complex double array) --- finial initial
    %                                               conditions to start timestepping.
            
            dof = size(fullmode_pod,1);
            a_initial = zeros(dof, 1);

            for i_init = 1:dof
                nxnz = fullmode_pod(i_init,2:3);
                Np = fullmode_pod(i_init,1);
            
                I = gp.find_wave_number_in_map(nxnz(1), nxnz(2), pod_wave);
            
                a_initial(i_init) = a_dns(t_initial, Np, I);
            
            end        

    end
    
    function [L_matrix,N_matrix] = built_matrix(L,N)

        % This function converts linear and nonlinear struct to matrix


        dof = size(L,2);

        L_matrix = zeros(dof,dof);
        
        for i = 1:dof
            disp(num2str(i))
            linear_loc = L(i).pod_pair;
            linear_coeff = L(i).coeff;
        
            L_matrix(i,linear_loc) = L_matrix(i,linear_loc) + linear_coeff;
        end


        tmp = 1:dof;
        nonlinear_index(:,2) = repmat(tmp',[dof,1]);
        nonlinear_index(:,1) = reshape(repmat(tmp,[dof,1]),[dof*dof,1]);

        NN = zeros(dof,dof*dof);

        for i = 1:dof

            % mxmz = N(i).mxmz;
            % kxkz = N(i).kxkz;
            coeff = N(i).coeff;
            mxmz_loc = N(i).a_mxmz_loc;
            kxkz_loc = N(i).a_kxkz_loc;
            
            for j = 1:size(mxmz_loc,1)
                for k = 1:size(mxmz_loc,2)
        
                    nonlinear_local_index = find(nonlinear_index(:,1) == mxmz_loc(j,k) & nonlinear_index(:,2) == kxkz_loc(j,k));
        
                    NN(i,nonlinear_local_index) = NN(i,nonlinear_local_index) + coeff(j,k);
        
        
                end
            end
        
        
        end

            N_matrix = NN;

        end






    function adot = eval_adot(t, a, L, N, nw_pod)
        
        disp(['t = ', num2str(t)]);

        a_len = size(L,2);

        adot_lin = zeros(a_len,1);
        adot_nonlin = adot_lin;

        % total number of permutations of the PODs
        perms_length = size(N(1).n_seq,1);        
        
        % linear coefficient
        for i_lin = 1:a_len

            for i_pod = 1:nw_pod

                anxnz = a(L(i_lin).pod_pair(i_pod));

                adot_lin(i_lin) = adot_lin(i_lin) + L(i_lin).coeff(i_pod)*anxnz;
            end

            % nonlinear coefficient

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

    function adot = eval_adot_fast(t, a, L, N, nw_pod)

        % This function evaluates the rate of change of modal time
        % coefficients. An improved faster version from eval_adot.
        % The differences enforced conjugate symmetry and loop unrolling
        % and some IO optimisation
        % Inputs: t (double) --- current time to evaluates rate of change (not
        %               computationally significant)
        %         a (complex double array) --- current modal time
        %                                      coefficients
        %         L (complex double array) --- Linear coefficient matrix
        %         N (complex double array) --- Nonlinear coefficient matrix
        %         nw_pod (double) --- number of POD basis for each wave
        %                             number
        % Outputs: adot (complex double array) --- current rate of change
        %                                          of modal time
        %                                          coefficiens

        disp(['t = ', num2str(t)]);

        a_len = size(L,2);

        adot_lin = zeros(a_len,1);
        adot_nonlin = adot_lin;

        conj_len = size(L(1).conj_ind_start,1);
        for i_conj = 1:conj_len
            st = L(1).conj_ind_start(i_conj);
            ed = L(1).conj_ind_end(i_conj); % ed for end 
            change_pod_ind = L(1).change_POD_ind(i_conj);

            for i_lin = change_pod_ind:st
                for i_pod = 1:nw_pod
    
                    anxnz = a(L(i_lin).pod_pair(i_pod));
    
                    adot_lin(i_lin) = adot_lin(i_lin) + L(i_lin).coeff(i_pod)*anxnz;
                end
    
                % nonlinear coefficient
    
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


    function adot = eval_adot_ultra_fast(t,a,L,N)

        disp(num2str(t))
        
        conj_len = size(L(1).conj_ind_start,1);

        for i_conj = 1:conj_len
            st = L(1).conj_ind_start(i_conj);
            ed = L(1).conj_ind_end(i_conj); % ed for end 
            change_pod_ind = L(1).change_POD_ind(i_conj);

            for i_dof = change_pod_ind:st

                tmp_pod_pair_lin = L(i_dof).pod_pair;
                t_adot_lin(i_dof,1) = L(i_dof).coeff(:)'*a(tmp_pod_pair_lin);
                
                a_mxmz_loc = N(i_dof).a_mxmz_loc(:);
                a_kxkz_loc = N(i_dof).a_kxkz_loc(:);
            
                AA = a(a_mxmz_loc).*a(a_kxkz_loc);
            
                coeff_nonlin = N(i_dof).coeff(:);
            
                t_adot_nonlin(i_dof,1) = coeff_nonlin'*AA;
            end
        
        end
        
        t_adot = t_adot_nonlin + t_adot_lin;
        
        conj_len = size(L(1).conj_ind_start,1);
        for i_conj = 1:conj_len
            st = L(1).conj_ind_start(i_conj);
            ed = L(1).conj_ind_end(i_conj); % ed for end 
        
            change_pod_ind = L(1).change_POD_ind(i_conj);
        
            t_adot(st:ed) = flip(conj(t_adot(change_pod_ind:st-2)));
        end
        
        adot = t_adot;

    end


    function [adot, adot_lin, adot_nonlin] = eval_adot_debug(t, a, L, N, nw_pod, a0)

    % Debug code for eval_adot.
    % Inputs: SEE eval_adot
    % Outputs: SEE eval_adot
    %          adot_lin (complex double array) --- linear contribution of
    %                                              Galerkin projection
    %          adot_nonlin (complex double array) --- nonlinear
    %          contribution of Galerkin projection
        
        disp(['t = ', num2str(t)]);

        a_len = size(L,2);

        adot_lin = zeros(a_len,1);
        adot_nonlin = adot_lin;

        conj_len = size(L(1).conj_ind_start,1);
        for i_conj = 1:conj_len
            st = L(1).conj_ind_start(i_conj);
            ed = L(1).conj_ind_end(i_conj); % ed for end 
            change_pod_ind = L(1).change_POD_ind(i_conj);

            for i_lin = change_pod_ind:st
                for i_pod = 1:nw_pod
    
                    anxnz = a(L(i_lin).pod_pair(i_pod));
    
                    adot_lin(i_lin) = adot_lin(i_lin) + L(i_lin).coeff(i_pod)*anxnz;
                end
    
                % nonlinear coefficient
    
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


    function [L, lam, diff, prod,Lv,Lu] = eval_L_m_verify(Re, ny, y, Lx, Lz, phi, n, m, nx, nz, pod_wave, dw, w)

        % This function evaluates the linear coefficient in Galerkin
        % projection. 
        % Inputs: Re (double) --- Reynolds number
        %         ny (int) --- Number of grid points in wall normal direction
        %         y  (double array) --- An array of wall normal coordinates. len(y) = ny
        %         Lx (double) --- Length of streamwise axis 
        %         Lz (double) --- Length of spanwise axis
        %         phi (complex double array) --- wall normal POD basis function
        %         n  (integer) --- wave number coupling in linear Galerkin projection
        %         m  (integer) --- same as above
        %         pod_wave (complex double array) --- column 1: all pod numbers
        %                                             column 2&3: all nx and nz wavenumbers 
        %         dw (double array) --- differentiate weight for Cheybeshev points
        %         w (double array) --- integration weights for Cheybeshev
        %                              points
        % Outputs: L (complex double array) --- full linear matrix for
        %                                       Galerkin projection
        %          lam (complex double array) --- Linear matrix for
        %                                         galerkin projection of the laminar state
        %          diff (complex double array) --- Linear matrix for
        %                                          Galerkin projection of diffusion (i.e \nabla^2 u) 
        %          prod (complex double array) --- Linear matrix for
        %                                          Galerkin projection of production term


        invRe = 1/Re;
        
        lam = 0;
        diff = 0;
        prod = 0;

        index = gp.find_wave_number_in_map(nx, nz, pod_wave);

        phiu_n = phi(1:3:end, n, index);
        phiv_n = phi(2:3:end, n, index);
        phiw_n = phi(3:3:end, n, index);

        phiu_m = phi(1:3:end, m, index);
        phiv_m = phi(2:3:end, m, index);
        phiw_m = phi(3:3:end, m, index);

        prod_int = conj(phiv_m).*phiu_n;

        for i = 1:ny
            prod = prod + w(i)*prod_int(i);
        end
        prod = -prod;

        % laminar state !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        lam_u_int = y.*conj(phiu_m).*phiu_n;
        lam_v_int = y.*conj(phiv_m).*phiv_n;
        lam_w_int = y.*conj(phiw_m).*phiw_n;
        
        lam_u = 0;
        lam_v = 0;
        lam_w = 0;

        for i = 1:ny
            lam_u = lam_u + w(i)*lam_u_int(i);
            lam_v = lam_v + w(i)*lam_v_int(i);
            lam_w = lam_w + w(i)*lam_w_int(i);
        end
        % !!!!!!!!!!!!!!!!!!! (I MISSED A NEGATIVE SIGN HERE BUT IT
        %                       WORKS)
        lam = 2*pi*1i*nx*(lam_u + lam_v + lam_w)/Lx;
        lam_V = 2*pi*1i*nx*lam_v/Lx;
        lam_U = 2*pi*1i*nx*lam_u/Lx;
        
        % diffusion/viscous
        kron = 0;
        if n == m
            kron = 1;
        end
        diff_1 = -invRe*((2*pi*nx/Lx)^2 + (2*pi*nz/Lz)^2)*kron;

        diff_u = 0;
        diff_v = 0;
        diff_w = 0;

        diff_u_int = dw(:,:,1)*conj(phiu_m).*(dw(:,:,1)*phiu_n);
        diff_v_int = dw(:,:,1)*conj(phiv_m).*(dw(:,:,1)*phiv_n);
        diff_w_int = dw(:,:,1)*conj(phiw_m).*(dw(:,:,1)*phiw_n);

        for i = 1:ny
            diff_u = diff_u + w(i)*diff_u_int(i);
            diff_v = diff_v + w(i)*diff_v_int(i);
            diff_w = diff_w + w(i)*diff_w_int(i);
        end
        diff = -invRe*(diff_u + diff_v + diff_w) + diff_1;
        diff_V = -invRe*(diff_v) + diff_1;
        diff_U = -invRe*(diff_u) + diff_1;
        
        % Combine all linear matrix components
        L = lam + diff + prod;
        Lv = lam_V + diff_V;
        Lu = lam_U + diff_U - prod;
        
        
    end

    function [L, fullmode_pod] = get_L_struct_verify(Re, Np, phi, lin_nonlin_ind, fullmode, pod_wave)

    % This function encasuplate all processes to compute linear and
    % nonlinear coefficients
    % Inputs: SEE get_L_N_struct
    % Ouputs: SEE get_L_N_struct

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
                [L(i).coeff(i_pod), L(i).lam(i_pod), L(i).diff(i_pod), L(i).prod(i_pod),L(i).Lv(i_pod),L(i).Lu(i_pod)] = ...
                    gp.eval_L_m_verify(Re, ny, y, Lx, Lz, phi, ...
                    L(i).n, L(i).n_seq(i_pod), nxnz(1), nxnz(2), pod_wave, dw, w);
            end
        end
    
        % rearrange field (just for aesthetic)
        % L = orderfields(L, [1:5,9,6:8,10:11]);
        % L = orderfields(L, [1:6,10,7:9,11]);
        % L = orderfields(L, [1:7,11,8:10]);
    end



    end

    % please stop. I have ran out of variable name......
end