
close all

REINIT_ALL_TOKAMAK_DATA=1;

if REINIT_ALL_TOKAMAK_DATA==1
    clear all
    load('../data_tokamak/motions_map_dimensions.mat')
    load('../data_tokamak/XZsmall_fields_tokamak_pre_collapse.mat')
    scale_X_small=scale_X;
    scale_Z_small=scale_Z;
    theta_map=theta_XZsmall_map;
    psi_map=psi_norm_XZsmall_map;

    initialize_folder_names;
    initialize_collapse_map_calculation_context
    rescaling_to_XZsmall_maps
    DT_INTERPOLATION_METHOD='quadratic'     % by default
    CALCULATE_VD_DATA_FILE=0;
end

load('../B_maps/B0311.mat');
load('../E_maps/E0311.mat');

frame_rank=32

psi_star_omega_map_half(:,:)=psi_star_2D_evol_lin(frame_rank,:,:);
% Using symmetry to reconstruct a poloidal turn
psi_star_omega_map_rank=zeros(size_r,NB_THETA);
psi_star_omega_map_rank(:,1:round(0.5*NB_THETA))=psi_star_omega_map_half(:,:);
psi_star_omega_map_rank(:,round(0.5*NB_THETA):NB_THETA)=psi_star_omega_map_half(:,round(0.5*NB_THETA):-1:1);

psi_star_omega_map=zeros(size_r,NB_THETA);
psi_star_PR_map=psi_star_omega_map';


% delta_phi_rank=phi_rank+1;
% psi_star_omega_map_rank_next(:,:)=rotate_map_phi(psi_star_omega_map_rank,delta_phi_rank);
% psi_star_PR_map_rank_next=psi_star_omega_map_rank_next';
% 
% delta_phi_rank=phi_rank-1;
% psi_star_omega_map_rank_prev(:,:)=rotate_map_phi(psi_star_omega_map_rank,delta_phi_rank);
% psi_star_PR_map_rank_prev=psi_star_omega_map_rank_prev';

        dr_data=reshape(dr_X_PR_map(:,1:size_r),NP*size_r,1);
        drX_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,dr_data(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
        drX_XZ_map=drX_XZ_map';
        drX_XZ_map(isnan(drX_XZ_map))=0;
        
        dr_data=reshape(dr_Z_PR_map(:,1:size_r),NP*size_r,1);
        drZ_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,dr_data(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
        drZ_XZ_map=drZ_XZ_map';
        drZ_XZ_map(isnan(drZ_XZ_map))=0;






phi_index=1

if PHI_OMEGA_RATIO==1
    phi_rank=phi_index
else
    phi_rank=PHI_OMEGA_RATIO*phi_index-1
end

delta_phi_rank=phi_rank;
psi_star_omega_map_phi(:,:)=rotate_map_phi(psi_star_omega_map_rank,delta_phi_rank);
psi_star_PR_map=psi_star_omega_map_phi';

% calculate_Btot_phi_rank;
% 
% 
% BHpolX_initial_XZsmall_map=zeros(sizeX,sizeZ);
% BHpolX_initial_XZsmall_map(:,:)=BHpol_X_map(Xinf:Xsup,Zinf:Zsup);
% 
% BHpolZ_initial_XZsmall_map=zeros(sizeX,sizeZ);
% BHpolZ_initial_XZsmall_map(:,:)=BHpol_Z_map(Xinf:Xsup,Zinf:Zsup);
% 
% BstarX_XZ_map=BpolX_XZ_map-BHpolX_initial_XZsmall_map;
% BstarZ_XZ_map=BpolZ_XZ_map-BHpolZ_initial_XZsmall_map;
% 
% BX_map(:,:)=bX_map_phi(phi_index,:,:);
% BZ_map(:,:)=bZ_map_phi(phi_index,:,:);
% Bphi_map(:,:)=sqrt(1-(BX_map.^2+BZ_map.^2));
% Btot_map(:,:)=Btot_map_phi(phi_index,:,:);
% BX_map=BX_map.*Btot_map;
% BZ_map=BZ_map.*Btot_map;
% % Bphi_map=Btor_PR_map(:,1:size_r);
% Bphi_map=Bphi_map.*Btot_map;

calculate_Bstar_phi_rank;

% B_data=reshape(BstarX_PR_map(:,1:size_r),NP*size_r,1);
% BstarX_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,B_data(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
% BstarX_XZ_map=BstarX_XZ_map';
% BstarX_XZ_map(isnan(BstarX_XZ_map))=0;
% 
% B_data=reshape(BstarZ_PR_map(:,1:size_r),NP*size_r,1);
% BstarZ_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,B_data(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
% BstarZ_XZ_map=BstarZ_XZ_map';
% BstarZ_XZ_map(isnan(BstarZ_XZ_map))=0;

% Bmin=min(min(BstarX_XZ_map));
% Bmax=max(max(BstarX_XZ_map));
% BstarX_XZ_map(BstarX_XZ_map<=0.8*Bmin)=0.8*Bmin;
% BstarX_XZ_map(BstarX_XZ_map>=0.8*Bmax)=0.8*Bmax;
% Bmin=min(min(BstarZ_XZ_map));
% Bmax=max(max(BstarZ_XZ_map));
% BstarX_XZ_map(BstarX_XZ_map<=0.8*Bmin)=0.8*Bmin;
% BstarX_XZ_map(BstarX_XZ_map>=0.8*Bmax)=0.8*Bmax;


Bstar_radial_phi1=BstarX_XZ_map.*drX_XZ_map+BstarZ_XZ_map.*drZ_XZ_map;
psi_XZ_map1=psi_XZ_map;
BstarX_XZ_map1=BstarX_XZ_map;
BstarZ_XZ_map1=BstarZ_XZ_map;




phi_index=round(0.5*Nomega)

if PHI_OMEGA_RATIO==1
    phi_rank=phi_index
else
    phi_rank=PHI_OMEGA_RATIO*phi_index-1
end

delta_phi_rank=phi_rank;
psi_star_omega_map_phi(:,:)=rotate_map_phi(psi_star_omega_map_rank,delta_phi_rank);
psi_star_PR_map=psi_star_omega_map_phi';

% calculate_Btot_phi_rank;
% 
% BHpolX_initial_XZsmall_map=zeros(sizeX,sizeZ);
% BHpolX_initial_XZsmall_map(:,:)=BHpol_X_map(Xinf:Xsup,Zinf:Zsup);
% 
% BHpolZ_initial_XZsmall_map=zeros(sizeX,sizeZ);
% BHpolZ_initial_XZsmall_map(:,:)=BHpol_Z_map(Xinf:Xsup,Zinf:Zsup);
% 
% BstarX_XZ_map=BpolX_XZ_map-BHpolX_initial_XZsmall_map;
% BstarZ_XZ_map=BpolZ_XZ_map-BHpolZ_initial_XZsmall_map;
% 
% BX_map(:,:)=bX_map_phi(phi_index,:,:);
% BZ_map(:,:)=bZ_map_phi(phi_index,:,:);
% Bphi_map(:,:)=sqrt(1-(BX_map.^2+BZ_map.^2));
% Btot_map(:,:)=Btot_map_phi(phi_index,:,:);
% BX_map=BX_map.*Btot_map;
% BZ_map=BZ_map.*Btot_map;
% % Bphi_map=Btor_PR_map(:,1:size_r);
% Bphi_map=Bphi_map.*Btot_map;

calculate_Bstar_phi_rank;

% Bmin=min(min(BstarX_XZ_map));
% Bmax=max(max(BstarX_XZ_map));
% BstarX_XZ_map(BstarX_XZ_map<=0.8*Bmin)=0.8*Bmin;
% BstarX_XZ_map(BstarX_XZ_map>=0.8*Bmax)=0.8*Bmax;
% Bmin=min(min(BstarZ_XZ_map));
% Bmax=max(max(BstarZ_XZ_map));
% BstarX_XZ_map(BstarX_XZ_map<=0.8*Bmin)=0.8*Bmin;
% BstarX_XZ_map(BstarX_XZ_map>=0.8*Bmax)=0.8*Bmax;

% B_data=reshape(BstarX_PR_map(:,1:size_r),NP*size_r,1);
% BstarX_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,B_data(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
% BstarX_XZ_map=BstarX_XZ_map';
% BstarX_XZ_map(isnan(BstarX_XZ_map))=0;
% 
% B_data=reshape(BstarZ_PR_map(:,1:size_r),NP*size_r,1);
% BstarZ_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,B_data(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
% BstarZ_XZ_map=BstarZ_XZ_map';
% BstarZ_XZ_map(isnan(BstarZ_XZ_map))=0;

Bstar_radial_phi2=BstarX_XZ_map.*drX_XZ_map+BstarZ_XZ_map.*drZ_XZ_map;
psi_XZ_map2=psi_XZ_map;
BstarX_XZ_map2=BstarX_XZ_map;
BstarZ_XZ_map2=BstarZ_XZ_map;


figure(1)
set(gca,'FontSize',24);
imagesc(scale_X,scale_Z, Bstar_radial_phi2');

colorbar;
axis xy
hold on;


%%
figure(2)
subplot_tight(2,2,1)
set(gca,'FontSize',26);
hold on
imagesc(scale_X,scale_Z ,(psi_XZ_map1)',[0 7]*1e-3);
contour(scale_X_small,scale_Z_small,theta_map',20,'k')
contour(scale_X_small,scale_Z_small,psi_map',20,'k')
colormap(flipud(gray))
xlim([-0.35 0.5])
ylim([-0.45 0.45])
axis xy
% axis square

xlabel('X')
ylabel('Z')
t_h=title('$$\psi_* (\varphi=0)$$','interpreter','latex');
set(t_h,'Interpreter','latex');
colorbar('North')


subplot_tight(2,2,2)
set(gca,'FontSize',26);
hold on
imagesc(scale_X,scale_Z , (psi_XZ_map2)',[0 7]*1e-3);
contour(scale_X_small,scale_Z_small,theta_map',20,'k')
contour(scale_X_small,scale_Z_small,psi_map',20,'k')
colormap(flipud(gray))
xlim([-0.35 0.5])
ylim([-0.45 0.45])
axis xy
% axis square

xlabel('X')
ylabel('Z')
t_h=title('$$\psi_* (\varphi=\pi)$$','interpreter','latex');
set(t_h,'Interpreter','latex');
colorbar('North')

subplot_tight(2,2,3)
set(gca,'FontSize',26);
hold on
imagesc(scale_X,scale_Z ,(Bstar_radial_phi1)',[-0.03 0.03]);
% contour(scale_X_small,scale_Z_small,theta_map',20,'k')
% contour(scale_X_small,scale_Z_small,psi_map',20,'k')
QUIV_FACTOR=29;
quiver(scale_X(1:QUIV_FACTOR:end),scale_Z(1:QUIV_FACTOR:end),...
    BstarX_XZ_map1(1:QUIV_FACTOR:end,1:QUIV_FACTOR:end)',BstarZ_XZ_map1(1:QUIV_FACTOR:end,1:QUIV_FACTOR:end)',3,'k','LineWidth',2)
colormap(jet)
xlim([-0.35 0.5])
ylim([-0.45 0.45])
axis xy
% axis square

xlabel('X')
ylabel('Z')
t_h=title('$$B_* (\varphi=0)$$','interpreter','latex');
set(t_h,'Interpreter','latex');
colorbar('North')


subplot_tight(2,2,4)
set(gca,'FontSize',26);
hold on
imagesc(scale_X,scale_Z , (Bstar_radial_phi2)',[-0.03 0.03]);
QUIV_FACTOR=29;
quiver(scale_X(1:QUIV_FACTOR:end),scale_Z(1:QUIV_FACTOR:end),...
    BstarX_XZ_map2(1:QUIV_FACTOR:end,1:QUIV_FACTOR:end)',BstarZ_XZ_map2(1:QUIV_FACTOR:end,1:QUIV_FACTOR:end)',3,'k','LineWidth',2)
% contour(scale_X_small,scale_Z_small,theta_map',20,'k')
% contour(scale_X_small,scale_Z_small,psi_map',20,'k')
colormap(jet)
xlim([-0.35 0.5])
ylim([-0.45 0.45])
axis xy
% axis square

xlabel('X')
ylabel('Z')
t_h=title('$$B_* (\varphi=\pi)$$','interpreter','latex');
set(t_h,'Interpreter','latex');
colorbar('North')