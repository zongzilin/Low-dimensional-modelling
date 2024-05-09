clear
clc;

close all

% Path to data and 'header'
addpath('../src');
addpath('../process/re400_1');

load('gp_nx3_nz3_N20_T3000_dt1_time_coeff.mat');
% load('gp_trs6_N6_T3000_dt1_time_coeff.mat');
load('gp_nx4_nz4_N1_T3000_dt1_time_coeff.mat','pod_wave','Re');

[x,y,z,nx,ny,nz,Lx,Lz] = modal_decomposition.read_geom();

% Initialise physical flow field
U = zeros(nx, ny, nz);
V = U;
W = U;

x_len = length(x);
y_len = length(y);
z_len = length(z);

i_t = 1200;
sqrtlxlz = 1/sqrt(Lx*Lz);
a_verify = a_gp(i_t,:);


for i_m = 1:size(fullmode_pod,1)
    
    wnx = fullmode_pod(i_m,2);
    wnz = fullmode_pod(i_m,3);

    i_pod = fullmode_pod(i_m,1);

    wnxz_index = gp.find_wave_number_in_map(wnx, wnz, pod_wave);

    phiU = phi(1:3:end, i_pod, wnxz_index);
    phiV = phi(2:3:end, i_pod, wnxz_index);
    phiW = phi(3:3:end, i_pod, wnxz_index);

    % Fills in physical flow field
    for i_x = 1:x_len
        for i_y = 1:y_len
            for i_z = 1:z_len

                exponen = exp(2*pi*1i*(wnx*x(i_x)/Lx + wnz*z(i_z)/Lz));

                U(i_x, i_y, i_z) = U(i_x, i_y, i_z) + sqrtlxlz*a_verify(i_m)*exponen*phiU(i_y);
                V(i_x, i_y, i_z) = V(i_x, i_y, i_z) + sqrtlxlz*a_verify(i_m)*exponen*phiV(i_y);
                W(i_x, i_y, i_z) = W(i_x, i_y, i_z) + sqrtlxlz*a_verify(i_m)*exponen*phiW(i_y);
            end
        end
    end
end

U = real(U);
V = real(V);
W = real(W);

[~,dw] = math.f_chebdiff(y_len, 2);
w = math.f_chebyshev_int_weight(y_len);

[dudx, dvdx, dwdx] = x_derivatives_spectral(x, y, z, Lx, Lz, a_verify, fullmode_pod, pod_wave, phi);
[dudz, dvdz, dwdz] = z_derivatives_spectral(x, y, z, Lx, Lz, a_verify, fullmode_pod, pod_wave, phi);

dudy = y_derivatives(U, x, y, z, dw);
dvdy = y_derivatives(V, x, y, z, dw);
dwdy = y_derivatives(W, x, y, z, dw);

d2udy2 = y_derivatives2(U, x, y, z, dw);
d2vdy2 = y_derivatives2(V, x, y, z, dw);
d2wdy2 = y_derivatives2(W, x, y, z, dw);

[d2udx2, d2vdx2, d2wdx2] = x_derivatives_spectral2(x, y, z, Lx, Lz, a_verify, fullmode_pod, pod_wave, phi);
[d2udz2, d2vdz2, d2wdz2] = z_derivatives_spectral2(x, y, z, Lx, Lz, a_verify, fullmode_pod, pod_wave, phi);


ydudx = zeros(x_len, y_len, z_len);
ydvdx = ydudx;
ydwdx = ydudx;

% laminar state
for i_x = 1:x_len
    for i_z = 1:z_len
        for i_y = 1:y_len
            ydudx(i_x, i_y, i_z) = y(i_y)*dudx(i_x, i_y, i_z);
            ydvdx(i_x, i_y, i_z) = y(i_y)*dvdx(i_x, i_y, i_z);
            ydwdx(i_x, i_y, i_z) = y(i_y)*dwdx(i_x, i_y, i_z);
        end
    end
end

nonlinear_u = -(U.*dudx + V.*dudy + W.*dudz);
nonlinear_v = -(U.*dvdx + V.*dvdy + W.*dvdz);
nonlinear_w = -(U.*dwdx + V.*dwdy + W.*dwdz);

linear_u = (d2udx2 + d2udy2 + d2udz2)/Re - V - ydudx;
linear_v = (d2vdx2 + d2vdy2 + d2vdz2)/Re - ydvdx;
linear_w = (d2wdx2 + d2wdy2 + d2wdz2)/Re - ydwdx;

u_fft = fft(fft(nonlinear_u + linear_u, nx, 1), nz, 3);
v_fft = fft(fft(nonlinear_v + linear_v, nx, 1), nz, 3);
w_fft = fft(fft(nonlinear_w + linear_w, nx, 1), nz, 3);

a_out = zeros(1,size(fullmode_pod,1));
for j = 1:size(fullmode_pod,1)
    [wave_x,wave_z] = modal_decomposition.applyPeriodicity(fullmode_pod(j,2),fullmode_pod(j,3),nx,nz);
    wnxz_phi_index = gp.find_wave_number_in_map(fullmode_pod(j,2), fullmode_pod(j,3), pod_wave);
    phi_index_out(j) = wnxz_phi_index;
        a = 0;
        a_u = 0;
        a_v = 0;
        a_w = 0;
        k = fullmode_pod(j,1);
        for l = 1:ny
            L_fft = 3*(l-1)+1:3*(l-1)+3;
            a_u = a_u + sum(w(l)*u_fft(wave_x,l,wave_z).*...
                        conj(phi(L_fft(1),k,wnxz_phi_index)));
            a_v = a_v + sum(w(l)*v_fft(wave_x,l,wave_z).*...
                        conj(phi(L_fft(2),k,wnxz_phi_index)));
            a_w = a_w + sum(w(l)*w_fft(wave_x,l,wave_z).*...
                        conj(phi(L_fft(3),k,wnxz_phi_index)));
            a = a_u + a_v + a_w;
        end
        a_out(j) = a;
end


% computes nonlinear coefficients
[adot_vis] = gp.eval_adot_fast(1, a_verify', L, N, Np);

a_out = sqrt(Lx*Lz)*a_out/(nx*nz);
a_dns_verify = [adot_vis'; a_out]';

plot(a_dns_verify(:,1),'ro','markersize',15,'LineWidth', 2)
hold on
plot(a_dns_verify(:,2),'kx','markersize',15,'LineWidth', 2)

% legend('flow reconstruction','gp code')

grid on
grid minor

a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',15)

set(gca,'TickLabelInterpreter','latex')





