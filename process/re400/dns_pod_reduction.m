clear;
clc

addpath('../../src');
addpath('../../re400/size1_1');

bnch = load('bench_dns_pod_modes.mat');

% t_start = 8000;
% t_end = 9000;
% 
% nt = t_end - t_start;
% 
% [pod_wave,pod_wave_conj] = modal_decomposition.genwave(8,8);
% 
% [~,geom.y,~,geom.nx,geom.ny,geom.nz,geom.Lx,geom.Lz] = modal_decomposition.read_geom;
% 
% w = math.f_chebyshev_int_weight(33);
% 
% U = zeros(geom.nx,geom.ny,geom.nz,nt);
% V = U;
% W = U;
% 
% % for i = 1:nt
% %     [U(:,:,:,i), V(:,:,:,i), W(:,:,:,i)] = modal_decomposition.read_snap(i);
% % 
% %     disp(num2str(i))
% % end
% % 
% % Uhat = fft(fft(U,geom.nx,1),geom.nz,3);
% % Vhat = fft(fft(V,geom.nx,1),geom.nz,3);
% % What = fft(fft(W,geom.nx,1),geom.nz,3);
% % 
% % for i = 1:nt
% %     [Id(:,:,:,:,i),P(:,:,:,:,i),R_sym(:,:,:,:,i),RP(:,:,:,:,i)] = modal_decomposition.bult_sym(geom.nx...
% %                                                             ,geom.ny,geom.nz,...
% %                                                             Uhat(:,:,:,i),Vhat(:,:,:,i),What(:,:,:,i));
% % end
% % 
% % % formulate autocorrelation matrix
% % total_wave = size(pod_wave,1);
% % R_tensor = zeros(3,3,total_wave,geom.ny,geom.ny);
% % 
% % for i_wave = 1:total_wave
% %     
% %     nxnz = pod_wave(i_wave,1:2);
% %     [i_nx,i_nz] = modal_decomposition.applyPeriodicity(nxnz(1),nxnz(2),16,16);
% % 
% %     for i = 1:geom.ny
% %         ui_Id = squeeze(Id(:,i_nx,i,i_nz,:));
% %         ui_P = squeeze(P(:,i_nx,i,i_nz,:));
% %         ui_R = squeeze(R_sym(:,i_nx,i,i_nz,:));
% %         ui_RP = squeeze(RP(:,i_nx,i,i_nz,:));
% %         for j = 1:geom.ny
% %             uj_Id = squeeze(Id(:,i_nx,j,i_nz,:));
% %             uj_P = squeeze(P(:,i_nx,j,i_nz,:));
% %             uj_R = squeeze(R_sym(:,i_nx,j,i_nz,:));
% %             uj_RP = squeeze(RP(:,i_nx,j,i_nz,:));
% %             R_tensor(:,:,i_wave,i,j) = R_tensor(:,:,i_wave,i,j) + ... 
% %                         ui_Id * uj_Id' + ui_P*uj_P' + ui_R*uj_R' + ui_RP*uj_RP';
% %         end
% %     end
% % 
% % end
% 
% total_wave = size(pod_wave,1);
% R = zeros(3,3,total_wave,geom.ny,geom.ny);
% 
% for i = t_start:t_end
% 
%     [U, V, W] = modal_decomposition.read_snap(i);
% 
%     disp(['summing snapshots ',num2str(i)])    
% 
%     Uhat = fft(fft(U,geom.nx,1),geom.nz,3);
%     Vhat = fft(fft(V,geom.nx,1),geom.nz,3);
%     What = fft(fft(W,geom.nx,1),geom.nz,3);    
% 
%     [Id,P,R_sym,RP] = modal_decomposition.bult_sym(geom.nx...
%                                                 ,geom.ny,geom.nz,...
%                                                 Uhat,Vhat,What);  
% 
%     for i_wave = 1:total_wave
% 
%         nxnz = pod_wave(i_wave,:);
%         [i_nx,i_nz] = modal_decomposition.applyPeriodicity(nxnz(1),nxnz(2),geom.nx,geom.nz);
% 
%         for i_y = 1:geom.ny
%             ui_Id = squeeze(Id(:,i_nx,i_y,i_nz));
%             ui_P = squeeze(P(:,i_nx,i_y,i_nz));
%             ui_R = squeeze(R_sym(:,i_nx,i_y,i_nz));
%             ui_RP = squeeze(RP(:,i_nx,i_y,i_nz));            
%             for j_y = 1:geom.ny
%                 uj_Id = squeeze(Id(:,i_nx,j_y,i_nz));
%                 uj_P = squeeze(P(:,i_nx,j_y,i_nz));
%                 uj_R = squeeze(R_sym(:,i_nx,j_y,i_nz));
%                 uj_RP = squeeze(RP(:,i_nx,j_y,i_nz));   
%                 R(:,:,i_wave,i_y,j_y) = R(:,:,i_wave,i_y,j_y) + ... 
%                             ui_Id * uj_Id' + ui_P*uj_P' + ui_R*uj_R' + ui_RP*uj_RP';                
%             end
%         end
% 
%     end
% end
% 
% 
% for i_wave = 1 : total_wave
%     for i = 1 : geom.ny
%         for j = 1 : geom.ny
%             
%             R(:,:,i_wave,i,j) = w(j,1)*R(:,:,i_wave,i,j)/(geom.nx*geom.nz)^2/4/nt;
%             
%         end
%     end
% end
% 
% A = zeros(3*geom.ny, 3*geom.ny, total_wave);
% for i_wave = 1 : total_wave
%     for i = 1 : geom.ny
%         ind_i = 3*(i-1) + 1 : 3*(i-1) + 3;
%         for j = 1 : geom.ny 
%             ind_j = 3*(j-1) + 1 : 3*(j-1) + 3;
%             A(ind_i,ind_j,i_wave) = R(:,:,i_wave,i,j);           
%         end
%     end    
%     A(:,:,i_wave) = geom.Lx *geom.Lz * A(:,:,i_wave);
% end

load('bench_dns_pod_modes_11aug.mat');

phi = zeros(3*geom.ny,3*geom.ny,total_wave);
lambda = zeros(3*geom.ny,3*geom.ny,total_wave);

% Integration Weights
% w = cat(2,w,w,w);
% w = reshape(w',3 * geom.ny ,1);

% Eigenvector under P Symmetry
PV = zeros(3*geom.ny,1);
PVi = zeros(3*geom.ny,1);

for i_wave = 1 : total_wave
    
    % Solve Eigenvalue Problem
    [EV,lambda(:,:,i_wave)] = eig(A(:,:,i_wave));
    
    % Enfore P Symmetry
    for n = 1 : 3 * geom.ny
        
        V = EV(:,n);
        Vi = 1i * V;
        
        % P Symmetry Summation : P(V) + V -> PV
        for i = 1 : geom.ny
            
            ind_Y = 3*(i-1) + 1 : 3*(i-1) + 3;
            ind_mY = 3*(geom.ny-i) + 1 : 3*(geom.ny-i) + 3;
            
            PV(ind_Y) = V(ind_Y,1) - conj(V(ind_mY,1));
            PVi(ind_Y) = Vi(ind_Y,1) - conj(Vi(ind_mY,1));

        end
        
        % V Norm
        V_norm = norm(PV);
        Vi_norm = norm(PVi);

        if V_norm >= Vi_norm
            V = PV / V_norm;
        else
            V = PVi / Vi_norm;
        end
        
        % Normalise with L2 Norm
        phi(:,n,i_wave) = V / sqrt(V' * (w .* V));
        
    end
    
    disp(['Eigenvalue problem of wavenumber (',num2str(pod_wave(i_wave,1)),...
          ',',num2str(pod_wave(i_wave,2)),') solved. Mode: ',num2str(i_wave),...
          '/',num2str(total_wave)]);
   
end

phi = cat(3,phi,conj(phi));
pod_wave = [pod_wave; pod_wave_conj];