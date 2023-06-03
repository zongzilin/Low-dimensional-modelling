clear;
clc;

% Path to data and 'header'
addpath('../../src');
addpath('../../re400/size1');

% load('dns_pod_modes_time_coeff.mat');
load('../../low_dim_modelling/dns_inner_prod_hkw_6_T10000_dt1_n_pod15_wave8_8.mat');
% load('../../dim_reduc/pod_modes_hkw_6_T10000_dt1_wave8_8.mat');
box = load('box_coeff_debug.mat');
% box = load('box_debug_time_coeff_integrate.mat');
% boxini = load('../../low_dim_modelling/NL_coeff_6.mat');
% box = load('box_debug_time_coeff_integrate_1.mat');

[x,y,z,nx,ny,nz,Lx,Lz] = modal_decomposition.read_geom;

[lin_nonlin_ind, fullmode, model] = gp.load_predefined_model(6);
% [lin_nonlin_ind, fullmode] = gp.load_model(3, 3);

const_lin_nonlin_ind = lin_nonlin_ind; % This line saves the index map to be used in ODE45 (nothing mathematical)

w = math.f_chebyshev_int_weight(ny);
[~,dw] = math.f_chebdiff(ny, 1);

Re = 400;

Np = 3; % Number of POD projections

% need a summing squence(or map) for all linear/nonlinear, all wavenumber
[L, fullpod, fullmode, fullmode_pod] = gp.prep_L_matrix(Np, fullmode);
L = gp.find_POD_pair_conj_index_in_model(L, Np, const_lin_nonlin_ind,fullmode);

for i = 1:size(L,2)
    size_n_seq = size(L(i).n_seq,2);
    nxnz = L(i).nxnz;
    for i_pod = 1:size_n_seq
        L(i).coeff(i_pod) = gp.eval_L_k(Re, ny, y, Lx, Lz, phi, ...
            L(i).n, L(i).n_seq(i_pod), nxnz(1), nxnz(2), pod_wave, dw, w);
    end
end


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

% all possible POD number summation
n_seq = perms(1:Np);
if Np ~= 1
tmp = repmat(1:Np,Np,1)';
n_seq = [n_seq; tmp];
end
size_n_seq = size(n_seq,1);
if Np == 1
    n_seq = [1 1];
    size_n_seq = size(n_seq,1);
end
clear tmp;

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
end


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

end

% prep initial condition from DNS
t_initial = 10000;
a_initial = zeros(size(N,2),1);
state_selector = [ 0, 0, 1; ...
                   0, 1, 1; 0,-1, 1; ...
                   0, 2, 1; 0,-2, 1; ...
                   1, 0, 1;-1, 0, 1; ...
                   1, 1, 1;-1,-1, 1; ...
                   1,-1, 1;-1, 1, 1];

for i_init = 1:size(N,2)
    nxnz = fullmode_pod(i_init,2:3);
    Np = fullmode_pod(i_init,1);

    I = gp.find_wave_number_in_map(nxnz(1), nxnz(2), pod_wave);

    a_initial(i_init) = a_dns(t_initial, Np, I);

end

% set ODE
end_time = 2000;
tspan = 0:1:end_time;
opt = odeset();

[t_a, a] = ode45(@(t,a) eval_adot(t, a, L, N, Np, fullmode_pod, a_initial),...
                  tspan, a_initial, opt);

dataTOplot = a(:,6);
plot(t_a,dataTOplot);
a_dns_to_plot = a_dns(t_initial:t_initial+end_time,1,1);


% a_initial = 1:size(N,2);
% a_initial = a_initial';
% [adot, adot_lin, adot_nonlin] = eval_adot(1, a_initial, L, N, Np, fullmode_pod, a_initial);

function [adot] = eval_adot(t, a, L, N, nw_pod, fullmode_pod, a0)
        
        disp(['t = ', num2str(t)]);

        a_len = size(a0,1);

        adot_lin = zeros(a_len,1);
        adot_nonlin = adot_lin;
        
        % linear coefficient
        for i_lin = 1:a_len

            for i_pod = 1:nw_pod

                anxnz = a(L(i_lin).pod_pair(i_pod));

                adot_lin(i_lin) = adot_lin(i_lin) + L(i_lin).coeff(i_pod)*anxnz;
            end

        end
        
        % total number of permutations of the PODs
        perms_length = size(N(1).n_seq,1);

        % nonlinear coefficient
        for i_nonlin = 1:a_len
            
            nxnz = N(i_nonlin).nxnz;
            kxkz_arr = N(i_nonlin).kxkz;
            mxmz_arr = N(i_nonlin).mxmz;

            kxkz_len = size(kxkz_arr,1);
            mxmz_len = size(mxmz_arr,1);

            for i = 1:kxkz_len
                mxmz = N(i_nonlin).mxmz(i,:);
                kxkz = N(i_nonlin).kxkz(i,:);                
                
                for i_perms = 1:perms_length

%                     m_loc = N(i_nonlin).a_mxmz_loc(i, i_perms);
%                     k_loc = N(i_nonlin).a_kxkz_loc(i, i_perms);
                    m_loc = gp.find_wave_number_pod_in_map(N(i_nonlin).n_seq(i_perms,1), mxmz(1), mxmz(2), fullmode_pod);
                    k_loc = gp.find_wave_number_pod_in_map(N(i_nonlin).n_seq(i_perms,2), kxkz(1), kxkz(2), fullmode_pod);

                    mult = a(m_loc)*a(k_loc);

                    adot_nonlin(i_nonlin) = adot_nonlin(i_nonlin) + N(i_nonlin).coeff(i,i_perms)*mult;
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


function [adot] = eval_adot_debug(t, a, L, N, nw_pod,fullmode_pod, a0, y, Lx, Lz, phi, pod_wave, dw, w)
        
        disp(['t = ', num2str(t)]);

        a_len = size(a0,1);

        adot_lin = zeros(a_len,1);
        adot_nonlin = adot_lin;
        
        % linear coefficient
        for i_lin = 1:a_len

            for i_pod = 1:nw_pod

                anxnz = a(L(i_lin).pod_pair(i_pod));

                adot_lin(i_lin) = adot_lin(i_lin) + L(i_lin).coeff(i_pod)*anxnz;
            end

        end
        
        % total number of permutations of the PODs
        perms_length = size(N(1).n_seq,1);

        % nonlinear coefficient
        for i_nonlin = 1:a_len
            
            nxnz = N(i_nonlin).nxnz;
            nx = nxnz(1);
            nz = nxnz(2);
            kxkz_arr = N(i_nonlin).kxkz;
            mxmz_arr = N(i_nonlin).mxmz;

            kxkz_len = size(kxkz_arr,1);
            mxmz_len = size(mxmz_arr,1);

            Np_c = N(i_nonlin).n;

            for i = 1:kxkz_len
                kxkz = N(i_nonlin).kxkz(i,:);  
                
                kx = kxkz(1);
                kz = kxkz(2);
                for i_m = 1:mxmz_len
                    mxmz = N(i_nonlin).mxmz(i_m,:);

                    mx = mxmz(1);
                    mz = mxmz(2);
                    for i_perms = 1:perms_length
    
                        m_np_c = N(i_nonlin).n_seq(i_perms,1);
                        k_np_c = N(i_nonlin).n_seq(i_perms,2);
        
                        m_loc = N(i_nonlin).a_mxmz_loc(i,m_np_c); % ######## %
                        k_loc = N(i_nonlin).a_kxkz_loc(i,k_np_c);

                        nnn = gp.eval_N_k(y, Lx, Lz, phi, Np_c, m_np_c, k_np_c, ...
                                                nx, nz, kx, kz, mx, mz,...
                                                pod_wave, dw, w);
    
                        adot_nonlin(i_nonlin) = adot_nonlin(i_nonlin) + nnn*a(m_loc)*a(k_loc);
                    end
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


