clear all
load('./data_tokamak/physics_constants.mat')
load('./data_tokamak/tokamak_PR_map.mat')
load('./data_tokamak/tokamak_map_dimensions.mat')
load('./data_tokamak/psi_star_evol.mat')
load('./data_tokamak/pressure_profile.mat')
load('./data_tokamak/q_profile.mat')
load('./data_tokamak/motions_map_dimensions.mat')
load('./data_tokamak/XZsmall_fields_tokamak_pre_collapse.mat')
rho_scale=radial_r_value_flux/max(radial_r_value_flux);
load('./data_tokamak/flux_geometry.mat', 'psi_XZ_map')

%%
[XXsmall ZZsmall]=meshgrid(scale_X,scale_Z);
finesse_data_X=reshape((Rpos_PR_map(:,1:size_r)-R0),NP*size_r,1);
finesse_data_Z=reshape(Z_PR_map(:,1:size_r),NP*size_r,1);

frame_rank=280
minus_phi_rank=1

psi_star_omega_map_half(:,:)=psi_star_2D_evol_lin(round(frame_rank/10)+1,:,:);

% Using symmetry to reconstruct a poloidal turn
psi_star_omega_map=zeros(size_r,NB_THETA);
psi_star_omega_map(:,1:round(0.5*NB_THETA))=psi_star_omega_map_half(:,:);
psi_star_omega_map(:,round(0.5*NB_THETA):NB_THETA)=psi_star_omega_map_half(:,round(0.5*NB_THETA):-1:1);
psi_star_omega_map_rank=psi_star_omega_map;

psi_star_map_phi_rank=zeros(size_r,NP);
psi_star_map_phi_rank(:,:)=[psi_star_omega_map(:,minus_phi_rank:NB_THETA)  psi_star_omega_map(:,1:minus_phi_rank-1)];
psi_star_map_phi_rank=psi_star_map_phi_rank';
psidata=reshape(psi_star_map_phi_rank(:,1:size_r),NB_THETA*size_r,1);

psistarmap=griddata(finesse_data_X,finesse_data_Z,psidata,XXsmall,ZZsmall,'cubic');
psistarmap=psistarmap';
psistarmap(isnan(psistarmap))=0;

%%
figure(1)
hold on; grid on
% contour(scale_X+R0,scale_Z,(psistarmap.*(psistarmap<-0.00257))',40);
% imagesc(scale_X+R0,scale_Z,(psistarmap.*(psistarmap<-0.00257))');

plot([X_axis+R0 X_axis+R0],[-R0 R0],'m--')
plot([0 X_axis+2*R0],[Z_axis Z_axis],'m--')

psi_XZ_map=max(psi_XZ_map,0);

contour(X_scale+R0,Z_scale,psi_XZ_map',6,'k--','linewidth',1)


contour(X_scale+R0,Z_scale,psi_XZ_map',[0 0],'k','linewidth',4)
contour(scale_X+R0,scale_Z,psistarmap',(-0.00257:0.001:0.0001),'g','linewidth',3);
contour(scale_X+R0,scale_Z,psistarmap',[0 0],'b','linewidth',3);
contour(scale_X+R0,scale_Z,psistarmap',-[0.00256 0.00256],'r','linewidth',3);

contour(scale_X+R0,scale_Z,(psistarmap.*(psistarmap<-0.002575))',(-0.003:0.00005:-0.002575));

% contour(scale_X,scale_Z,psi_XZsmall_map',[psi_q1 psi_q1],'g','linewidth',1.5)

set(gca,'FontSize',26)
xlabel('R (m)')
ylabel('Z (m)')
% colorbar
axis image
xlim([-0.6 0.55]+R0)
ylim([-1.0 0.85])

% title('Poloidal angle \theta')

%%
xlim([-0.2 0.28]+R0)
ylim([-0.3 0.3])

ah=annotation('textbox', [0.4,0.4,0.1,0.1],...
           'String', '1');
set(ah,'fontsize',30);
set(ah,'FontWeight','bold');
set(ah,'LineStyle','none');

set(ah,'color',[0 0.7 0]);

ah=annotation('textbox', [0.47,0.78,0.1,0.1],...
           'String', '2');
set(ah,'fontsize',30);
set(ah,'FontWeight','bold');
set(ah,'LineStyle','none');

set(ah,'color',[0 0.0 0.7]);


ah=annotation('textbox', [0.65,0.48,0.1,0.1],...
           'String', '3');
set(ah,'fontsize',30);
set(ah,'FontWeight','bold');
set(ah,'LineStyle','none');

set(ah,'color',[0.9 0.0 0.0]);