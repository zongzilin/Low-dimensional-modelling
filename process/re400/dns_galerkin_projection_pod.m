clear;
clc;
close all

% Path to data and 'header'
addpath('../../src');
addpath('../../re400/size1');

load('../../low_dim_modelling/dns_inner_prod_hkw_6_T10000_dt1_n_pod15_wave8_8.mat','phi','pod_wave','a_dns');

% [lin_nonlin_ind, fullmode, model] = gp.load_predefined_model(6);
[lin_nonlin_ind, fullmode] = gp.load_model(3, 3);

Re = 400;
Np = 1;

[L, N, fullmode_pod] = gp.get_L_N_struct(Re, Np, phi, lin_nonlin_ind, fullmode, pod_wave);

% prep initial condition from DNS
t_initial = 10500;
a_initial = zeros(size(N,2),1);


for i_init = 1:size(N,2)
    nxnz = fullmode_pod(i_init,2:3);
    Np = fullmode_pod(i_init,1);

    I = gp.find_wave_number_in_map(nxnz(1), nxnz(2), pod_wave);

    a_initial(i_init) = a_dns(t_initial, Np, I);

end

% set ODE
end_time = 1000;
tspan = 0:0.1:end_time;
opt = odeset();

[t_a, a] = ode45(@(t,a) gp.eval_adot(t, a, L, N, Np, a_initial),...
                  tspan, a_initial, opt);

a_gp = a;
dataTOplot = a(:,25);
plot(t_a,dataTOplot,'linewidth',1.5);
a_dns_to_plot = a_dns(t_initial:t_initial+end_time,1,1);
hold on
plot(0:1:end_time,a_dns_to_plot,'linewidth',1.5)






