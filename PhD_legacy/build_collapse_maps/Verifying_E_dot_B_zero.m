
close all

REINIT_ALL_TOKAMAK_DATA=1;

if REINIT_ALL_TOKAMAK_DATA==1
    clear all
    initialize_folder_names;
    initialize_collapse_map_calculation_context
    rescaling_to_XZsmall_maps
    DT_INTERPOLATION_METHOD='quadratic'     % by default
    CALCULATE_VD_DATA_FILE=0;
end
% Epot_evol=-Epot_evol;

load('../B_maps/B0301.mat');
load('../E_maps/E0301.mat');

frame_rank=31
phi_index=1

psi_star_omega_map_half(:,:)=psi_star_2D_evol_lin(frame_rank,:,:);
% Using symmetry to reconstruct a poloidal turn
psi_star_omega_map_rank=zeros(size_r,NB_THETA);
psi_star_omega_map_rank(:,1:round(0.5*NB_THETA))=psi_star_omega_map_half(:,:);
psi_star_omega_map_rank(:,round(0.5*NB_THETA):NB_THETA)=psi_star_omega_map_half(:,round(0.5*NB_THETA):-1:1);

psi_star_omega_map=zeros(size_r,NB_THETA);
psi_star_PR_map=psi_star_omega_map';

if PHI_OMEGA_RATIO==1
    phi_rank=phi_index
else
    phi_rank=PHI_OMEGA_RATIO*phi_index-1
end

% if phi_index==1
%     phi_rank=1
% end

psi_star_omega_map(:,:)=rotate_map_phi(psi_star_omega_map_rank,phi_rank);
psi_star_PR_map=psi_star_omega_map';

delta_phi_rank=phi_rank+1;
psi_star_omega_map_rank_next(:,:)=rotate_map_phi(psi_star_omega_map_rank,delta_phi_rank);
psi_star_PR_map_rank_next=psi_star_omega_map_rank_next';

delta_phi_rank=phi_rank-1;
psi_star_omega_map_rank_prev(:,:)=rotate_map_phi(psi_star_omega_map_rank,delta_phi_rank);
psi_star_PR_map_rank_prev=psi_star_omega_map_rank_prev';

calculate_Btot_phi_rank;


BHpolX_initial_XZsmall_map=zeros(sizeX,sizeZ);
BHpolX_initial_XZsmall_map(:,:)=BHpol_X_map(Xinf:Xsup,Zinf:Zsup);

BHpolZ_initial_XZsmall_map=zeros(sizeX,sizeZ);
BHpolZ_initial_XZsmall_map(:,:)=BHpol_Z_map(Xinf:Xsup,Zinf:Zsup);


BstarX_XZ_map=BpolX_XZ_map-BHpolX_initial_XZsmall_map;
BstarZ_XZ_map=BpolZ_XZ_map-BHpolZ_initial_XZsmall_map;

% calculate_Bstar_phi_rank;

%
% B_data=reshape(BstarX_PR_map(:,1:size_r),NP*size_r,1);
% BstarX_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,B_data(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
% BstarX_XZ_map=BstarX_XZ_map';
% BstarX_XZ_map(isnan(BstarX_XZ_map))=0;
%
% B_data=reshape(BstarZ_PR_map(:,1:size_r),NP*size_r,1);
% BstarZ_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,B_data(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
% BstarZ_XZ_map=BstarZ_XZ_map';
% BstarZ_XZ_map(isnan(BstarZ_XZ_map))=0;
%
% verif_dot_product=EX_XZ_map.*BstarX_XZ_map+EZ_XZ_map.*BstarZ_XZ_map+Efield_3_XZ_map.*Bfield_3_XZmap;
%verif_dot_product=abs(EX_XZ_map.*Bstar_X_XZ_map+EZ_XZ_map.*Bstar_Z_XZ_map+Efield_3_RZmap.*Bfield_3_RZmap);
%verif_dot_product=abs(ER_XZ_zoom_map.*Bpol_X_zoom_map+EZ_XZ_zoom_map.*Bpol_Z_zoom_map+Ephi_XZ_zoom_map_num.*Bphi_XZ_zoom_map);


Epot_omega_map(:,:)=Epot_evol(frame_rank,:,:);

E_potential_PR_map_phi=zeros(NB_THETA,NB_THETA,size_r);
for index=1:NB_THETA
    %         phi_rank=round(2*phi_index-1);
    E_potential_PR_map_phi(index,:,:)=rotate_map_phi(Epot_omega_map,index)';
end

psi_star_dot_omega_map(:,:)=psi_star_dot_evol(frame_rank,:,:);

psi_star_dot_PR_map_phi=zeros(NB_THETA,NB_THETA,size_r);
for index=1:NB_THETA
    %         phi_rank=round(2*phi_index-1);
    psi_star_dot_PR_map_phi(index,:,:)=rotate_map_phi(psi_star_dot_omega_map,index)';
end



psi_star_dot_PR_map(:,:)=psi_star_dot_PR_map_phi(phi_rank,:,:);
% load psi_star_dot_PR_map_ref3.mat
% psi_star_dot_PR_map(:,:)=psi_star_dot_PR_map_ref3(:,1:size_r);
% 
E_data=reshape(psi_star_dot_PR_map(:,1:size_r),NP*size_r,1);
% Efield_3_XZ_map=griddata(finesse_data_X(1:size_r*NP),finesse_data_Z(1:size_r*NP),E_data,XX_small,ZZ_small,'cubic');
Efield_3_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,E_data(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
Efield_3_XZ_map=Efield_3_XZ_map';
Efield_3_XZ_map(isnan(Efield_3_XZ_map))=0;
%
Bfield_3_XZmap=Bphi_XZsmall_map./Rpos_XZsmall_map;

E_potential_PR_map(:,:)=E_potential_PR_map_phi(phi_rank,:,:);
if phi_rank>1
    E_potential_PR_map_prev(:,:)=E_potential_PR_map_phi(phi_rank-1,:,:);
else
    E_potential_PR_map_prev(:,:)=E_potential_PR_map_phi(end-1,:,:);
end
E_potential_PR_map_next(:,:)=E_potential_PR_map_phi(phi_rank+1,:,:);

calculate_Efield_phi_rank;
% EX_PR_map=squeeze(Efield_X_map_phi(phi_index,:,:));
% EZ_PR_map=squeeze(Efield_Z_map_phi(phi_index,:,:));

        E_data=reshape(EZ_PR_map(:,1:size_r),NP*size_r,1);
        EZ_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,E_data(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
        EZ_XZ_map=EZ_XZ_map';
        
        E_data=reshape(EX_PR_map(:,1:size_r),NP*size_r,1);
        EX_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,E_data(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
        EX_XZ_map=EX_XZ_map';

        E_data=reshape(Ephi_PR_map(:,1:size_r),NP*size_r,1);
        Ephi_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,E_data(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
        Ephi_XZ_map=Ephi_XZ_map';

        E_data=reshape(grad_Phi_tor_PR_map(:,1:size_r),NP*size_r,1);
        gPhi_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,E_data(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
        gPhi_XZ_map=gPhi_XZ_map';

        %  
        
%         B_data=reshape(BpolZ_PR_map(:,1:size_r),NP*size_r,1);
%         BZ_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,B_data(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
%         BZ_XZ_map=BZ_XZ_map';
%         
%         B_data=reshape(BpolX_PR_map(:,1:size_r),NP*size_r,1);
%         BX_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,B_data(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
%         BX_XZ_map=BX_XZ_map';
% 
%         B_data=reshape(Btor_PR_map(:,1:size_r),NP*size_r,1);
%         Bphi_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,B_data(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
%         Bphi_XZ_map=Bphi_XZ_map';
%                
        
% E_data=reshape(EX_PR_map(:,1:size_r),NP*size_r,1);
% EX_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,E_data(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
% EX_XZ_map=EX_XZ_map';

% EX_map=EX_PR_map;
% EZ_map=EZ_PR_map;

EX_map(:,:)=Efield_X_map_phi(phi_index,:,:);
EZ_map(:,:)=Efield_Z_map_phi(phi_index,:,:);

% BX_map=BpolX_PR_map;
% BZ_map=BpolZ_PR_map;
% Btot_map=Btot_PR_map;
%
BX_map(:,:)=bX_map_phi(phi_index,:,:);
BZ_map(:,:)=bZ_map_phi(phi_index,:,:);

% Opposite ephi and Bphi
Bphi_map(:,:)=sqrt(1-(BX_map.^2+BZ_map.^2));

Btot_map(:,:)=Btot_map_phi(phi_index,:,:);
BX_map=BX_map.*Btot_map;
BZ_map=BZ_map.*Btot_map;
% Bphi_map=Btor_PR_map(:,1:size_r);
Bphi_map=Bphi_map.*Btot_map;


% verif_dot_product=(BpolX_PR_map.*EX_PR_map+BpolZ_PR_map.*EZ_PR_map+Ephi_PR_map.*Bphi_map);
verif_dot_product=(BX_map.*EX_map+BZ_map.*EZ_map+(psi_star_dot_PR_map./Rpos_PR_map(:,1:size_r)-grad_Phi_tor_PR_map).*Bphi_map);
verif_dot_product=(BpolX_PR_map.*EX_map+BpolZ_PR_map.*EZ_map+(Ephi_PR_map(:,1:size_r)).*Bphi_map);

% Ephi_PR_map=(BX_map.*EX_map+BZ_map.*EZ_map)./Bphi_map;
% verif_dot_product=(-BX_map.*EX_map-BZ_map.*EZ_map+Ephi_PR_map.*Bphi_map);

%%
figure(1)
set(gca,'FontSize',22);
imagesc(verif_dot_product,[-250 250]);
% imagesc(scale_X,scale_Z,verif_dot_product',[-1 1]*40);axis xy;
% xlabel('X (m)');
% ylabel('Z (m)');
% xlim([-0.4 0.5]);
% ylim([-0.6 0.6]);
colorbar;
% axis xy square
hold on;

disp('Total cumulated errors');
disp(sum(sum(verif_dot_product)));
%%
calculate_Bstar_phi_rank;
BX_XZ_map=BstarX_XZ_map+BHpolX_initial_XZsmall_map;
BZ_XZ_map=BstarZ_XZ_map+BHpolZ_initial_XZsmall_map;

%%
figure(2);
verif_dot_product_XZ=EX_XZ_map.*BX_XZ_map+EZ_XZ_map.*BZ_XZ_map+Ephi_XZ_map.*Bphi_XZsmall_map;
% verif_dot_product_XZ=EX_XZ_map.*BHpolX_initial_XZsmall_map+EZ_XZ_map.*BHpolZ_initial_XZsmall_map-gPhi_XZ_map.*Bphi_XZsmall_map;
% verif_dot_product_XZ=EX_XZ_map.*BstarX_XZ_map+EZ_XZ_map.*BstarZ_XZ_map+Efield_3_XZ_map.*Bfield_3_XZmap;
imagesc(verif_dot_product_XZ',[-50 50]);axis xy
colorbar;

verif_dot_product_XZ(isnan(verif_dot_product_XZ))=0;
disp('Total cumulated errors');
disp(sum(sum(abs(verif_dot_product_XZ))));

psi_star_dot_PR_map_recalc=-(Rpos_PR_map(:,1:size_r)./Bphi_map).*(BstarX_PR_map.*EX_PR_map+BstarZ_PR_map.*EZ_PR_map);


% BstarX_PR_map=BpolX_PR_map-BHpol_X_PR_map(:,1:size_r);
% BstarZ_PR_map=BpolZ_PR_map-BHpol_Z_PR_map(:,1:size_r);

%%
figure(3);
verif_dot_product_PR=(BstarX_PR_map.*EX_PR_map+BstarZ_PR_map.*EZ_PR_map+psi_star_dot_PR_map.*Bphi_map./Rpos_PR_map(:,1:size_r));
imagesc((verif_dot_product_PR)',[-25 25]);axis xy
colorbar;

verif_dot_product_PR(isnan(verif_dot_product_PR))=0;
disp('Total cumulated errors');
disp(sum(sum(abs(verif_dot_product_PR))));
