clear;
clc;

% Path to data and 'header'
addpath('../../src');
addpath('../../re400/size1');

load('gp_9_mode_pod_3_time_coeff.mat');

% Number of POD modes we want to project into
n_pod = 1;

sim_time = size(a_gp,1);

% Initialise physical flow field
U = zeros(nx, ny, nz, sim_time);
V = U;
W = U;

x_len = length(x);
y_len = length(y);
z_len = length(z);

% Search pod mode in full mode
mode_count = size(fullmode, 1);
pod_mode_map = zeros(mode_count,1);
for i_mc = 1:mode_count
    pod_mode_map(i_mc) = gp.find_wave_number_in_map(fullmode(i_mc,1), ...
        fullmode(i_mc,2), pod_wave);
end
phi = phi(:,:,pod_mode_map);

for i_t = 1:sim_time
    disp(['Timestep: ', num2str(i_t)])
    
    for i_pod = 1:n_pod

        for i_m = 1:mode_count
            
            wnx = fullmode(i_m,1);
            wnz = fullmode(i_m,2);

            phiU = phi(1:3:end, n_pod, i_m);
            phiV = phi(2:3:end, n_pod, i_m);
            phiW = phi(3:3:end, n_pod, i_m);

            % Fills in physical flow field
            for i_x = 1:x_len
                for i_y = 1:y_len
                    for i_z = 1:z_len

                        exponen = exp(2*pi*1i*( wnx*x(i_x)/Lx + wnz*z(i_z)/Lz ));
                        sqrtlxlz = 1/sqrt(Lx*Lz);

                        U(i_x, i_y, i_z, i_t) = U(i_x, i_y, i_z, i_t) + sqrtlxlz*a_gp(i_t,i_m)*exponen*phiU(i_y);
                        V(i_x, i_y, i_z, i_t) = V(i_x, i_y, i_z, i_t) + sqrtlxlz*a_gp(i_t,i_m)*exponen*phiV(i_y);
                        W(i_x, i_y, i_z, i_t) = W(i_x, i_y, i_z, i_t) + sqrtlxlz*a_gp(i_t,i_m)*exponen*phiW(i_y);
                    end
                end
            end
    
        end
    end

end

U = real(U).^2;
V = real(V).^2;
W = real(W).^2;

U = mean(U,4);
U = mean(U,3);
U = mean(U,1);

V = mean(V,4);
V = mean(V,3);
V = mean(V,1);

W = mean(W,4);
W = mean(W,3);
W = mean(W,1);

U_rms = sqrt(U);
V_rms = sqrt(V);
W_rms = sqrt(W);


