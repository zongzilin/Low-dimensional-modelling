classdef modal_redProjlib
    methods(Static)

    function [x,y,z,nx,ny,nz,Lx,Lz] = read_geom
    
        uperb = ncread('u0.nc','Velocity_X');
        
        nx = size(uperb,1);
        ny = size(uperb,2);
        nz = size(uperb,3);
        
        
        x = ncread('u0.nc','X');
        y = ncread('u0.nc','Y');
        z = ncread('u0.nc','Z');
        
        Lx = x(end) + x(2);
        Lz = z(end) + z(2);
    
    end        
    
    function [U, V, W] = read_snap(n)
        % READS U V W OF CURRENT T = n
        
        dset = ['u' num2str(n) '.nc'];
    
        U = ncread(dset,'Velocity_X');
        V = ncread(dset,'Velocity_Y');
        W = ncread(dset,'Velocity_Z');
    end    

    function [x,y] = applyPeriodicity(wave_x,wave_z,nx,nz)

        if wave_x >= 0 
            wave_x = wave_x + 1;
        else
            wave_x = nx + wave_x + 1;
        end
    
        if wave_z >= 0
            wave_z = wave_z + 1;
        else 
            wave_z = nz + wave_z + 1;
        end
    
        x = wave_x;
        y = wave_z;
    end
    
    function [Id,P,R,RP] = bult_sym(nx,ny,nz,Uhat_xz,Vhat_xz,What_xz)
    
        % P - flipped about (nx = 0, nz = 0)
        % R - flipped about (nz = 0)
        % RP - rotation about nz = 0 by \pi 
    
        % UVW-hats - (wave number in x, ny, wave number in z);
    
        Id = zeros([3,size(Uhat_xz)]);
        Id(1,:,:,:) = Uhat_xz;
        Id(2,:,:,:) = Vhat_xz;
        Id(3,:,:,:) = What_xz;
    
        ind_x = [1 nx:-1:2];
        ind_y = ny:-1:1;
        ind_z = [1 nz:-1:2];
    
        P(1,:,:,:) = -Id(1,ind_x,ind_y,ind_z);
        P(2,:,:,:) = -Id(2,ind_x,ind_y,ind_z);
        P(3,:,:,:) = -Id(3,ind_x,ind_y,ind_z);
    
        ind_x = 1:nx;
        ind_y = 1:ny;
        ind_z = [1 nz:-1:2];
    
        R(1,:,:,:) = Id(1,ind_x,ind_y,ind_z);
        R(2,:,:,:) = Id(2,ind_x,ind_y,ind_z);
        R(3,:,:,:) = -Id(3,ind_x,ind_y,ind_z);
    
        ind_x = [1 nx:-1:2];
        ind_y = ny:-1:1;
        ind_z = 1:nz;
    
        RP(1,:,:,:) = -Id(1,ind_x,ind_y,ind_z);
        RP(2,:,:,:) = -Id(2,ind_x,ind_y,ind_z);
        RP(3,:,:,:) = Id(3,ind_x,ind_y,ind_z);
    end    
    
    function [pod_wave,pod_wave_conj] = genwave(nx,nz)
            
        % FOR nz == 0
        pod_wave_nz0 = zeros(nz,2);
        pod_wave_nz0(:,2) = 0;
        pod_wave_nz0(:,1) = nx:-1:1;
    
        % FOR nx == 0
        pod_wave_nx0 = zeros(nx,2);
        pod_wave_nx0(:,1) = 0;
        pod_wave_nx0(:,2) = nz:-1:1;
    
        % FOR nx ~= 0 && nz ~= 0
        pod_wave = [];
    
        for i = 1 : nx
            for j = 1 : nz
                
                pod_wave = [pod_wave; i, j];
        
            end
        end
    
        pod_wave = [pod_wave; pod_wave(:,1), -pod_wave(:,2)];
    
        pod_wave = [0 0; pod_wave_nx0; pod_wave_nz0; pod_wave];
    
        pod_wave_conj = -pod_wave;
    
    
    end
    end
end