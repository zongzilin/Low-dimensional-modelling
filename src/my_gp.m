clear;
clc;

load('dns_inner_prod_hkw_6_T10000_dt1_n_pod15_wave8_8.mat');

[x,y,z,nx,ny,nz,Lx,Lz] = modal_decomposition.read_geom;

[lin_nonlin_ind, fullmode, model] = gp.load_predefined_model(6);
% [lin_nonlin_ind, fullmode] = gp.load_model(3, 3);

w = math.f_chebyshev_int_weight(ny);
[~,dw] = math.f_chebdiff(ny, 1);

Re = 400;
n = 1;
m = 1;
k = 1;

L = struct();
for i = 1:size(lin_nonlin_ind,2)

    L(i).nxnz = fullmode(i,1:2);
    nxnz = L(i).nxnz;
    
    L(i).coeff = gp.eval_L(Re, ny, y, Lx, Lz, phi, n, m, nxnz(1), nxnz(2), pod_wave, w);

end

N = struct();
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

        N(i).coeff(j) = gp.eval_N_k(y, Lx, Lz, phi, n, m, k, nxnz(1), nxnz(2), kx, kz, mx, mz, pod_wave,...
                    dw, w);

    end

end

a0 = zeros(size(fullmode,1),1);
t_init_cond = 10000;

for i = 1:size(fullmode)
    
    nxnz = fullmode(i,:);
    I = gp.find_wave_number_in_map(nxnz(1), nxnz(2), pod_wave);

    a0(i) = a_dns(t_init_cond, 1, I);

end

tspan = 0:1:3000;
opt = odeset();

[t_a, a] = ode45(@(t,a) func_hdl(t, a, L, N, lin_nonlin_ind, fullmode),...
                  tspan, a0, opt);
ddd = a(:,9);
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
















