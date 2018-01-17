reset_data_analysis_environment;
load('NBI_Phot_data.mat')
scale_X_P=scale_X;
scale_Z_P=scale_Z;


filename=strcat(DATA_FOLDER,'XZsmall_fields_tokamak_pre_collapse.mat');
load(filename);

initialize_xi_map_calculation_context;
rescaling_to_XZsmall_maps;
NBI_Phot_XZ_map_small=interp2(scale_X_P,scale_Z_P,NBI_Phot_XZ_map',XX_small,ZZ_small);
NBI_Phot_XZ_map_small=NBI_Phot_XZ_map_small';
calculate_gradPhot;


%%
figure(1)
set(gca,'fontsize',20)
VECSP=32;

hold on
% contourc(scale_X_P,scale_Z_P,NBI_Phot_XZ_map',30);
% contour(scale_X,scale_Z,NBI_Phot_XZ_map_small',20);
contour(scale_X,scale_Z,NBI_Phot_XZ_map_small',30,'linewidth',3);
contour(scale_X,scale_Z,psi_norm_XZsmall_map',10,'k','linewidth',0.5)
contour(scale_X,scale_Z,psi_norm_XZsmall_map',[psi_rank_q1 psi_rank_q1],'b.','linewidth',5)

colormap summer
% colormap('summer')
[scXX scZZ]=meshgrid(scale_X,scale_Z);

quiver(scXX(1:VECSP:end,1:VECSP:end-VECSP),scZZ(1:VECSP:end,1:VECSP:end-VECSP),gradP_X(1:VECSP:end-VECSP,1:VECSP:end)',gradP_Z(1:VECSP:end-VECSP,1:VECSP:end)',0.8,'r')
% contour(scale_X_P,scale_Z_P,NBI_Phot_XZ_map',30,'linewidth',3);
% colormap hot

contour(scale_X,scale_Z,NBI_Phot_XZ_map_small',30,'linewidth',0.5);
colorbar

xlim([-0.15 0.25])
ylim([-0.25 0.25])
xlabel('X (m)')
ylabel('Z (m)')
