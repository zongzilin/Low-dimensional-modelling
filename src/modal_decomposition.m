classdef modal_decomposition
    methods(Static)
    
        function [x,y,z,nx,ny,nz,Lx,Lz] = read_geom
        
            uperb = ncread('u0.nc','Velocity_X');
            
            nx = size(uperb,1);
            ny = size(uperb,2);
            nz = size(uperb,3);
            
            
            x = ncread('u0.nc','X');
            y = ncread('u0.nc','Y');
            z = ncread('u0.nc','Z');
            
            Lx = x(end) + x(2);
            Lz = z(end) + z(2);
        
        end        
        
        function [U, V, W] = read_snap(n)
            % READS U V W OF CURRENT T = n
            
            dset = ['u' num2str(n) '.nc'];
        
            U = ncread(dset,'Velocity_X');
            V = ncread(dset,'Velocity_Y');
            W = ncread(dset,'Velocity_Z');
        end    
    
        function [x,y] = applyPeriodicity(wave_x,wave_z,nx,nz)
    
            if wave_x >= 0 
                wave_x = wave_x + 1;
            else
                wave_x = nx + wave_x + 1;
            end
        
            if wave_z >= 0
                wave_z = wave_z + 1;
            else 
                wave_z = nz + wave_z + 1;
            end
        
            x = wave_x;
            y = wave_z;
        end
        
        function [Id,P,R,RP] = bult_sym(nx,ny,nz,Uhat_xz,Vhat_xz,What_xz)
        
            % P - flipped about (nx = 0, nz = 0)
            % R - flipped about (nz = 0)
            % RP - rotation about nz = 0 by \pi 
        
            % UVW-hats - (wave number in x, ny, wave number in z);
        
            Id = zeros([3,size(Uhat_xz)]);
            Id(1,:,:,:) = Uhat_xz;
            Id(2,:,:,:) = Vhat_xz;
            Id(3,:,:,:) = What_xz;
        
            ind_x = [1 nx:-1:2];
            ind_y = ny:-1:1;
            ind_z = [1 nz:-1:2];
        
            P(1,:,:,:) = -Id(1,ind_x,ind_y,ind_z);
            P(2,:,:,:) = -Id(2,ind_x,ind_y,ind_z);
            P(3,:,:,:) = -Id(3,ind_x,ind_y,ind_z);
        
            ind_x = 1:nx;
            ind_y = 1:ny;
            ind_z = [1 nz:-1:2];
        
            R(1,:,:,:) = Id(1,ind_x,ind_y,ind_z);
            R(2,:,:,:) = Id(2,ind_x,ind_y,ind_z);
            R(3,:,:,:) = -Id(3,ind_x,ind_y,ind_z);
        
            ind_x = [1 nx:-1:2];
            ind_y = ny:-1:1;
            ind_z = 1:nz;
        
            RP(1,:,:,:) = -Id(1,ind_x,ind_y,ind_z);
            RP(2,:,:,:) = -Id(2,ind_x,ind_y,ind_z);
            RP(3,:,:,:) = Id(3,ind_x,ind_y,ind_z);
        end    
        
        function [pod_wave,pod_wave_conj] = genwave(nx,nz)
                
            % FOR nz == 0
            pod_wave_nz0 = zeros(nz,2);
            pod_wave_nz0(:,2) = 0;
            pod_wave_nz0(:,1) = nx:-1:1;
        
            % FOR nx == 0
            pod_wave_nx0 = zeros(nx,2);
            pod_wave_nx0(:,1) = 0;
            pod_wave_nx0(:,2) = nz:-1:1;
        
            % FOR nx ~= 0 && nz ~= 0
            pod_wave = [];
        
            for i = 1 : nx
                for j = 1 : nz
                    
                    pod_wave = [pod_wave; i, j];
            
                end
            end
        
            pod_wave = [pod_wave; pod_wave(:,1), -pod_wave(:,2)];
        
            pod_wave = [0 0; pod_wave_nx0; pod_wave_nz0; pod_wave];
        
            pod_wave_conj = -pod_wave;
        
        
        end
        
        function [phi, pod_wave, proj_basis_amplitude] = get_dns_proj_basis(path2dns,no_snap, n_wave, n_pod)

            if isa(path2dns, string) == true % ?????
                error('path to dns data must be string');
                return
            end


            [~,~,~,nx,ny,nz,Lx,Lz] = modal_decomposition.read_geom;

            nt = no_snap;   % number of timestep

            [pod_wave,pod_wave_conj] = modal_decomposition.genwave(n_wave,n_wave);
            
            % INITIALISE AUTOCORRELATION TENSOR
            R_R = zeros(3,3,size(pod_wave,1),ny,ny);
            Id_R = R_R;
            RP_R = R_R;
            P_R = R_R;
            R = R_R;
            r = R_R; % TEMPORARY AUTOCORRELATION AT EACH DT
            
            % CHEBY WEIGHT FOR EIGENPROBLEM DISCRETISATION
            w = math.f_chebyshev_int_weight(ny);
            
            for i_snap = 1:no_snap
                
                [U, V, W] = modal_decomposition.read_snap(i_snap);
            
                % PERFORM FFT ON X AND Z
                U_fft = fft(fft(U,nx,1),nz,3);
                V_fft = fft(fft(V,nx,1),nz,3);
                W_fft = fft(fft(W,nx,1),nz,3);
            
                % SYMMETRY
                [Id,P,R_sym,RP] = modal_decomposition.bult_sym(nx,ny,nz,U_fft,V_fft,W_fft);
            
                % SNAPSHOT AUTOCORRELATION TENSOR R
                for i = 1:size(pod_wave,1)
                    wave_x = pod_wave(i,1);
                    wave_z = pod_wave(i,2);
            
                    [wave_x,wave_z] = modal_decomposition.applyPeriodicity(wave_x,wave_z,nx,nz);
            
                    for k = 1:ny
                        Id_Ui = Id(:,wave_x,k,wave_z);
                        R_Ui = R_sym(:,wave_x,k,wave_z);
                        RP_Ui = RP(:,wave_x,k,wave_z);
                        P_Ui = P(:,wave_x,k,wave_z);
            
                        for j = 1:ny
                            Id_Uj = Id(:,wave_x,j,wave_z);
                            R_Uj = R_sym(:,wave_x,j,wave_z);
                            RP_Uj = RP(:,wave_x,j,wave_z);
                            P_Uj = P(:,wave_x,j,wave_z);
            
                            Id_R(:,:,i,k,j) = Id_R(:,:,i,k,j) + Id_Ui*Id_Uj';
                            R_R(:,:,i,k,j) = R_R(:,:,i,k,j) + R_Ui*R_Uj';
                            RP_R(:,:,i,k,j) = RP_R(:,:,i,k,j) + RP_Ui*RP_Uj';
                            P_R(:,:,i,k,j) = P_R(:,:,i,k,j) + P_Ui*P_Uj';
                        end
                    end
                end
            
                r = P_R + RP_R + R_R + Id_R;
                R = R + r;
            
                % NORMALISE AUTOCORRELATION MATRIX
                for i = 1:size(pod_wave,1)
                    for j = 1:ny
                        for k = 1:ny
                            R(:,:,i,j,k) = w(j,1) * r(:,:,i,j,k)/(nx*nz)^2/4/nt;
                        end
                    end
                end
                
                disp(['TIMESTEP: ',num2str(i_snap),' DONE'])
            end
            
            
            % DISCRETISE EIGENVALUE PROBLEM
            A = zeros(3*ny, 3*ny, size(pod_wave,1));
            
            for i_wave = 1:size(pod_wave,1)
                for i = 1:ny
                    ind_i = 3*(i-1)+1:3*(i-1)+3;
                    for j = 1:ny
                        ind_j = 3*(j-1)+1:3*(j-1)+3;
                        A(ind_i,ind_j,i_wave) = R(:,:,i_wave,i,j);
                    end
                end
                A(:,:,i_wave) = Lx*Lz*A(:,:,i_wave);
            end
            
            n_waves = size(pod_wave,1);
            
            phi = zeros(3*ny,3*ny,n_waves);
            lambda = zeros(3*ny,3*ny,n_waves);
            
            w = cat(2,w,w,w);
            w = reshape(w',3*ny,1);
            
            PV = zeros(3*ny,1);
            PVi = PV;
            
            for i_wave = 1:n_waves
                [EV,lambda(:,:,i_wave)] = eig(A(:,:,i_wave));
            
                for n = 1:3*ny
                    V = EV(:,n);
                    Vi = 1i*V;
            
                    for i = 1:ny
                        ind_Y = 3*(i-1)+1:3*(i-1)+3;
                        ind_mY = 3*(ny-i)+1:3*(ny-i)+3;
            
                        PV(ind_Y) = V(ind_Y,1) - conj(V(ind_mY,1));
                        PVi(ind_Y) = Vi(ind_Y,1) - conj(Vi(ind_mY,1));
            
                    end
            
                    V_norm = norm(PV);
                    Vi_norm = norm(PVi);
            
                    if V_norm >= Vi_norm
                        V = PV/V_norm;
                    else
                        V = PVi/Vi_norm;
                    end
                    
                    phi(:,n,i_wave) = V / sqrt(V' * (w .* V));
            
                end
            end
            
            phi = cat(3,phi,conj(phi));
            pod_wave = [pod_wave; pod_wave_conj];

            proj_basis_amplitude = modal_decomposition.get_proj_basis_amplitude(phi, 1000, n_pod);

        end

        function proj_basis_amplitude = get_proj_basis_amplitude(phi, nt, n_pod)

                [~,~,~,nx,ny,nz,Lx,Lz] = modal_decomposition.read_geom;
                
                [wxz,~] = modal_decomposition.genwave(8, 8);
                
                w = math.f_chebyshev_int_weight(ny);
                
                A = zeros(nt,n_pod,size(wxz,1));
                
                for i = 1:nt
                    
                    [U,V,W] = modal_decomposition.read_snap(i);
                
                    U = fft(fft(U,nx,1),nz,3);
                    V = fft(fft(V,nx,1),nz,3);
                    W = fft(fft(W,nx,1),nz,3);
                
                    for j = 1:size(wxz,1)
                        [wave_x,wave_z] = modal_decomposition.applyPeriodicity(wxz(j,1),wxz(j,2),nx,nz);
                        for k = 1:n_pod
                            a = 0;
                            a_u = 0;
                            a_v = 0;
                            a_w = 0;
                            for l = 1:ny
                                L = 3*(l-1)+1:3*(l-1)+3;
                                a_u = a_u + sum(w(l)*U(wave_x,l,wave_z).*...
                                            conj(phi(L(1),k,j)));
                                a_v = a_v + sum(w(l)*V(wave_x,l,wave_z).*...
                                            conj(phi(L(2),k,j)));
                                a_w = a_w + sum(w(l)*W(wave_x,l,wave_z).*...
                                            conj(phi(L(3),k,j)));
                                a = a_u + a_v + a_w;
                            end
                            A(i,k,j) = a;
                        end
                    end
                    disp(['t = ',num2str(i)])
                end
                
                a_dns = sqrt(Lx*Lz)*A/(nx*nz);
                a_dns_conj = conj(a_dns);
                
                a_dns_out = cat(3,a_dns,a_dns_conj);
                a_dns = a_dns_out;
                
                proj_basis_amplitude = a_dns;
                
            
        end
    end
end
