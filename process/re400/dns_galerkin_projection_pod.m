clear;
clc;

% Path to data and 'header'
addpath('../../src');
addpath('../../re400/size1');

load('dns_pod_modes_time_coeff.mat');

[x,y,z,nx,ny,nz,Lx,Lz] = modal_decomposition.read_geom;

% [lin_nonlin_ind, fullmode, model] = gp.load_predefined_model(6);
[lin_nonlin_ind, fullmode] = gp.load_model(4, 4);

w = math.f_chebyshev_int_weight(ny);
[~,dw] = math.f_chebdiff(ny, 1);

Re = 400;

% Number of pod modes that we want to project into
nw_pod = 1;
% Number of pod modes that we are using to project
n_pod = 1:nw_pod;
no_pod = length(n_pod);

% Construct pod index before going into loops
% L = struct();
% fullmode = repmat(fullmode,nw_pod,1);
% for i = 1:size(fullmode, 1)
%     L(i).nxnz = fullmode(i,1:2);
% end
% fullpod = 1:nw_pod;
% fullpod = repmat(fullpod, size(fullmode,1)/nw_pod,1);
% fullpod = reshape(fullpod, [size(fullpod,1)*size(fullpod,2) 1]);
% for i = 1:size(fullpod, 1)
%     L(i).n = fullpod(i); % This step we write "the number of pod modes that we want to project into" into convient location
%     L(i).coeff = 0; % And we also allocate memory for linear coefficients.
% end

[L, fullpod, fullmode] = gp.prep_L_matrix(nw_pod, fullmode);

for i = 1:size(lin_nonlin_ind,2)*nw_pod
        
    L(i).nxnz = fullmode(i,1:2);
    nxnz = L(i).nxnz;
    n = L(i).n;
    
    for i_pod = 1:no_pod
        L(i).coeff = L(i).coeff + gp.eval_L(Re, ny, y, Lx, Lz, phi, n, n_pod(i_pod), nxnz(1), nxnz(2), pod_wave, w);
    end
end

N = struct();
for i = 1:size(fullpod, 1)
    N(i).n = fullpod(i); 
    N(i).coeff = 0;
end

% Inflate lin_nonlin_ind to match POD number of POD expansions
totl_mode = size(lin_nonlin_ind,2);
for i = 1:totl_mode*(nw_pod-1)
    lin_nonlin_ind(i + totl_mode).nxnz = lin_nonlin_ind(i).nxnz;
    lin_nonlin_ind(i + totl_mode).mxmz = lin_nonlin_ind(i).mxmz;
    lin_nonlin_ind(i + totl_mode).kxkz = lin_nonlin_ind(i).kxkz;
end

% Nonlinear coefficient initialisation
for i = 1:size(lin_nonlin_ind,2)
    
    N(i).nxnz = lin_nonlin_ind(i).nxnz;
    N(i).kxkz = lin_nonlin_ind(i).kxkz;
    N(i).mxmz = lin_nonlin_ind(i).mxmz;
   
    nxnz = N(i).nxnz;

    mxmz_arr = N(i).mxmz;
    kxkz_arr = N(i).kxkz;

    N(i).coeff = zeros(size(mxmz_arr,1),1);

    for j = 1:size(mxmz_arr,1)
        kx = kxkz_arr(j,1);
        kz = kxkz_arr(j,2);

        mx = mxmz_arr(j,1);
        mz = mxmz_arr(j,2);

        for i_pod_m = 1:no_pod
            for i_pod_k = 1:no_pod
                N(i).coeff(j) = N(i).coeff(j) + gp.eval_N_k(y, Lx, Lz, phi, N(i).n , n_pod(i_pod_m), n_pod(i_pod_k), nxnz(1), nxnz(2), kx, kz, mx, mz, pod_wave,...
                                                dw, w);
            end
        end

    end
end


t_init_cond = 2000;
a0 = struct();
for i_pod = 1:size(fullpod,1)
    a0(i_pod).n = fullpod(i_pod);
end

for i = 1:size(fullpod,1)
    
    nxnz = fullmode(i,:);
    I = gp.find_wave_number_in_map(nxnz(1), nxnz(2), pod_wave);
    
    a0(i).nxnz = nxnz;
    a0(i).time_coeff = a_dns(t_init_cond, a0(i).n , I);

end



tspan = 0:0.5:3000;
opt = odeset('RelTol',1e-4,'AbsTol',1e-6);

a_in = vertcat(a0.time_coeff);
[t_a, a] = ode45(@(t,a) func_hdl(t, a, L, N, lin_nonlin_ind, fullmode),...
                  tspan, a_in, opt);
a_gp = -a;

clear a;

ddd = a_gp(:,6);
plot(t_a,ddd);

function adot = func_hdl(t, a, L, N, lin_nonlin_ind, fullmode)

    adot_lin = zeros(size(lin_nonlin_ind,2),1);
    adot_nonlin = adot_lin;
    
    disp(['t = ', num2str(t)])

    for i = 1:(size(lin_nonlin_ind,2) + 1)/2
        nxnz = L(i).nxnz;

        I = gp.find_wave_number_in_map(nxnz(1), nxnz(2), fullmode);
        
        anxnz = a(I);
        adot_lin(i) = L(i).coeff*anxnz;
    end

    for i = 1:(size(lin_nonlin_ind,2) + 1)/2
        kxkz_arr = N(i).kxkz;
        mxmz_arr = N(i).mxmz;

        for j = 1:size(kxkz_arr,1)
            kxkz = kxkz_arr(j,:);
            mxmz = mxmz_arr(j,:);

            Ik = gp.find_wave_number_in_map(kxkz(1),kxkz(2), fullmode);
            Im = gp.find_wave_number_in_map(mxmz(1),mxmz(2), fullmode);

            ak = a(Ik);
            am = a(Im);

            adot_nonlin(i) = adot_nonlin(i) + N(i).coeff(j)*ak*am;
        end
    end


    adot = adot_lin + adot_nonlin;
    
    conj_ind = (size(lin_nonlin_ind,2) + 1)/2;
    adot(conj_ind + 1:end) = flip(conj(adot(1:conj_ind - 1)));

end














