clear
clc

close all

% Path to data and 'header'
addpath('../src');
addpath('../re400/size1');

load('../low_dim_modelling/dns_inner_prod_hkw_6_T10000_dt1_n_pod15_wave8_8.mat','pod_wave','phi');
load('../verify/a_dns_T10000_full.mat');
a_dns = squeeze(a_dns);

[lin_nonlin_ind, fullmode, model] = gp.load_predefined_model(9);
% [lin_nonlin_ind, fullmode] = gp.load_model(1, 2);

Re = 400;
Np = 4;

% get linear and nonlinear coefficients
[L, N, fullmode_pod] = gp.get_L_N_struct(Re, Np, phi, lin_nonlin_ind, fullmode, pod_wave);

a_verify = zeros(1, size(N,2));
for i = 1:size(N,2)
    nxnz = fullmode_pod(i,2:3);
    Np = fullmode_pod(i,1);

    I = gp.find_wave_number_in_map(nxnz(1), nxnz(2), pod_wave);

    a_verify(i) = a_dns(Np, I);

end

[adot, adot_lin, adot_nonlin] = gp.eval_adot_debug(1, a_verify', L, N, Np, a_verify');

adot_nonlin_verify = conj(a_verify)'.*adot_nonlin;

plot(adot_nonlin_verify, 'rs','markersize',15)

sum(adot_nonlin_verify)







