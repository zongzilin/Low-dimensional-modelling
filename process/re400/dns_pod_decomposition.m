clear;
clc;

% Path to data and 'header'
addpath('../../src');
addpath('../../re400/size1');

clear;
clc;

close all

[x,y,z,nx,ny,nz,Lx,Lz] = modal_decomposition.read_geom;

dt = 0.01;
no_snap = 3;    % NUMBER OF SNAPSHOT TO LOAD
nt = no_snap;   % TAKE DT = 1

% No. of Fourier mode to compute in x and z
n_wave = 8;
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

a = lambda(:,:,1);
disp(real(a(1,1)));

phi = cat(3,phi,conj(phi));
pod_wave = [pod_wave; pod_wave_conj];
