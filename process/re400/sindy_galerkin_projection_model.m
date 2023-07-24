cd '/home/zilin/Desktop/Project/process/re400'

clear;
clc;
close all

% Path to data and 'header'
addpath('../../src');
addpath('../../re400/size1');

data_name = 'box_data/gp_nx3_nz3_N1_T3000_dt1_time_coeff.mat';

load('../../low_dim_modelling/dns_inner_prod_hkw_6_T10000_dt1_n_pod10_wave8_8.mat','a_dns','phi','pod_wave');
load(data_name,'Np','L', 'N','fullmode_pod');

Re = 400;
% [lin_nonlin_ind, fullmode, model] = gp.load_predefined_model(9);
[lin_nonlin_ind, fullmode] = gp.load_model(3,3);

% computes the dns a_dot
ts = 1;
te = 10000;
const_a_dns = a_dns;
a_dns = a_dns(ts:te,:,:);

% use finite difference to construct dns a_dot
% use spline ppval to compute get adot?
a_dot_dns = sindy.eval_adot_spline(a_dns);

% residual/library
T = sindy.residual_galerkin_R(Np, L, N, a_dns, a_dot_dns, fullmode_pod,pod_wave);
[D, fullmode_pod] = sindy.get_D_struct(Np, phi, lin_nonlin_ind, fullmode, pod_wave);

for i = 1:te
%     e(i) = squeeze(conj(a_dns(i,1,:)))'*squeeze(a_dns(i,1,:));
    e(i) = abs(a_dns(i,1,1));
end
% e = sqrt(e);
theta = sindy.library_galerkin_R(e, a_dns, D, pod_wave);

% dof
dof = size(fullmode_pod,1);

% sindy solve
sparsity_knob = 0.0001;
iter = 10;
Xi = zeros(dof,1);
for i_s = 1:dof
    Xi(i_s) = sindy.sindy_solve(1, sparsity_knob, iter, T(:,i_s), theta(:,i_s));
end

c = Xi;

% prep initial condition from DNS
t_initial = 3000;
a_initial = zeros(size(N,2),1);

for i_init = 1:size(N,2)
    nxnz = fullmode_pod(i_init,2:3);
    Np = fullmode_pod(i_init,1);

    I = gp.find_wave_number_in_map(nxnz(1), nxnz(2), pod_wave);

    a_initial(i_init) = const_a_dns(t_initial, Np, I);

end

% set ODE
end_time = 5000;
tspan = 0:1:end_time;
opt = odeset();

[t_a, a] = ode45(@(t,a) sindy.eval_adot_galerkin_R(t, a, L, N, D, c, Np),...
                  tspan, a_initial, opt);

a_gp = a;

vis = gp.find_wave_number_pod_in_map(1, 0, 0, fullmode_pod);


dns_M = load('dns_M.mat','M');

M = post.modal_energy(a, Np, L, 1.75*pi, 1.2*pi);
M_1 = post.modal_energy_from_file(data_name);

subplot(2,1,1)
plot(dns_M.M(1:end_time,1,1),'linewidth',1.5)
hold on
plot(M(1:end-1,vis),'linewidth',1.5);

grid on
grid minor

legend('dns', 'sindy')

subplot(2,1,2)
plot(dns_M.M(1:end_time,1,1),'linewidth',1.5)
hold on
plot(M_1(1:end_time,vis),'linewidth',1.5);

set(gcf,'position',[500,1000,1200,600])

legend('dns', 'gp')

clear e theta a_dns a_dot_dns const_a_dns T