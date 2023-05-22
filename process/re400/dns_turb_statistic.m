clear;
clc;

% Path to data and 'header'
addpath('../../src');
addpath('../../re400/size1');

sim_time = 2000;
[x,y,z,x_len,y_len,z_len,Lx,Lz] = modal_decomposition.read_geom;

Ubase = load('Ubase.asc');

uperb = zeros(x_len, y_len, z_len, sim_time);
vperb = uperb;
wperb = uperb;

for i_t = 1:sim_time

    disp(['Time step:', num2str(i_t)])

    [u, v, w] = read_snap(i_t);

    uperb(:,:,:,i_t) = u.^2;
    vperb(:,:,:,i_t) = v.^2;
    wperb(:,:,:,i_t) = w.^2;


end

% Average in time
uperb = mean(uperb,4);
vperb = mean(vperb,4);
wperb = mean(wperb,4);

% Average in x
uperb = squeeze(mean(uperb,1));
vperb = squeeze(mean(vperb,1));
wperb = squeeze(mean(wperb,1));

% Average in z
uperb = squeeze(mean(uperb,2));
vperb = squeeze(mean(vperb,2));
wperb = squeeze(mean(wperb,2));

