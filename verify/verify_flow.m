clear
clc;

close all

% Path to data and 'header'
addpath('../src');
addpath('../process/re400');

load('gp_nx3_nz3_N3_T3000_dt1_time_coeff.mat');
% load('gp_trs6_N1_T3000_dt1_time_coeff.mat');

[~,~,~,~,~,~,Lx,Lz] = modal_decomposition.read_geom();

% Number of POD modes we want to project into
n_pod = 3;

sim_time = size(a_gp,1);

% Initialise physical flow field
U = zeros(nx, ny, nz);
V = U;
W = U;

x_len = length(x);
y_len = length(y);
z_len = length(z);

% Search unique pod mode in full mode
% mode_count = size(fullmode, 1);
% pod_mode_map = zeros(mode_count,1);
% for i_mc = 1:mode_count
%     pod_mode_map(i_mc) = gp.find_wave_number_in_map(fullmode(i_mc,1), ...
%         fullmode(i_mc,2), pod_wave);
% end
% pod_mode_map = unique(pod_mode_map);
% phi_fft = phi;
% phi = phi(:,:,pod_mode_map);


i_t = 1;
sqrtlxlz = 1/sqrt(Lx*Lz);
a_bench = a_gp(i_t,:);
a_verify = 1:size(a_gp,2);

for i_pod = 1:n_pod
    for i_m = 1:size(a_gp,2)
        
        wnx = fullmode(i_m,1);
        wnz = fullmode(i_m,2);

        wnxz_index = gp.find_wave_number_in_map(wnx, wnz, pod_wave);
    
        phiU = phi(1:3:end, i_pod, wnxz_index);
        phiV = phi(2:3:end, i_pod, wnxz_index);
        phiW = phi(3:3:end, i_pod, wnxz_index);
    
        % Fills in physical flow field
        for i_x = 1:x_len
            for i_y = 1:y_len
                for i_z = 1:z_len
    
                    exponen = exp(2*pi*1i*(wnx*x(i_x)/Lx + wnz*z(i_z)/Lz));
    
                    U(i_x, i_y, i_z) = U(i_x, i_y, i_z) + sqrtlxlz*a_verify(i_t,i_m)*exponen*phiU(i_y);
                    V(i_x, i_y, i_z) = V(i_x, i_y, i_z) + sqrtlxlz*a_verify(i_t,i_m)*exponen*phiV(i_y);
                    W(i_x, i_y, i_z) = W(i_x, i_y, i_z) + sqrtlxlz*a_verify(i_t,i_m)*exponen*phiW(i_y);
                end
            end
        end
    end
end

U = real(U);
V = real(V);
W = real(W);


% Nonlinear u \dot \grad u
% x derivatives (streamwise)
dudx = x_derivatives(U, x, y, z);
dvdx = x_derivatives(V, x, y, z);
dwdx = x_derivatives(W, x, y, z);

d2udx2 = x_derivatives2(U, x, y, z);
d2vdx2 = x_derivatives2(V, x, y, z);
d2wdx2 = x_derivatives2(W, x, y, z);

% y derivatives (wall normal)
[~, dw] = math.f_chebdiff(y_len, 1);
[~, dw2] = math.f_chebdiff(y_len, 2);
dudy = y_derivatives(U, x, y, z, dw);
dvdy = y_derivatives(V, x, y, z, dw);
dwdy = y_derivatives(W, x, y, z, dw);

d2udy2 = y_derivatives2(U, x, y, z, dw2);
d2vdy2 = y_derivatives2(V, x, y, z, dw2);
d2wdy2 = y_derivatives2(W, x, y, z, dw2);

% z derivatives (spanwise)
dudz = z_derivatives(U, x, y, z);
dvdz = z_derivatives(V, x, y, z);
dwdz = z_derivatives(W, x, y, z);

d2udz2 = z_derivatives2(U, x, y, z);
d2vdz2 = z_derivatives2(V, x, y, z);
d2wdz2 = z_derivatives2(W, x, y, z);

nonlinear_u = U.*dudx + V.*dudy + W.*dudz;
nonlinear_v = U.*dvdx + V.*dvdy + W.*dvdz;
nonlinear_w = U.*dwdx + V.*dwdy + W.*dwdz;

ns = nonlinear_u - (d2udx2 + d2udy2 + d2udz2)/400;

% fft of nonlinear velocities
u_fft = fft(fft(nonlinear_u, nx, 1), nz, 3);
v_fft = fft(fft(nonlinear_v, nx, 1), nz, 3);
w_fft = fft(fft(nonlinear_w, nx, 1), nz, 3);



for j = 1:size(fullmode_pod,1)
    [wave_x,wave_z] = modal_decomposition.applyPeriodicity(fullmode_pod(j,2),fullmode_pod(j,3),nx,nz);
    wnxz_phi_index = gp.find_wave_number_in_map(fullmode_pod(j,2), fullmode_pod(j,3), pod_wave);
        a = 0;
        a_u = 0;
        a_v = 0;
        a_w = 0;
        k = fullmode_pod(j,1);
        for l = 1:ny
            L = 3*(l-1)+1:3*(l-1)+3;
            a_u = a_u + sum(w(l)*u_fft(wave_x,l,wave_z).*...
                        phi(L(1),k,wnxz_phi_index));
            a_v = a_v + sum(w(l)*v_fft(wave_x,l,wave_z).*...
                        conj(phi(L(2),k,wnxz_phi_index)));
            a_w = a_w + sum(w(l)*w_fft(wave_x,l,wave_z).*...
                        conj(phi(L(3),k,wnxz_phi_index)));
            a = a_u + a_v + a_w;
        end
        A(j) = a;
end

a_dns_verify = [a_verify; sqrt(Lx*Lz)*A/(nx*nz)];





