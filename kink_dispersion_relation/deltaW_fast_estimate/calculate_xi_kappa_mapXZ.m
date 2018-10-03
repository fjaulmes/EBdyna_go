
%  reset_data_analysis_environment;
initialize_xi_map_calculation_context;
rescaling_to_XZsmall_maps;

Bphi_PR_map=Btor_PR_map(:,1:size_r);

load('../E_maps/E0111.mat');
load('../B_maps/B0111.mat');

frame_rank=12
phi_index=1



psi_star_dot_omega_map(:,:)=psi_star_dot_evol(frame_rank,:,:);

psi_star_dot_PR_map_phi=zeros(NB_THETA,NB_THETA,size_r);
for index=1:NB_THETA
    %         phi_rank=round(2*phi_index-1);
    psi_star_dot_PR_map_phi(index,:,:)=rotate_map_phi(psi_star_dot_omega_map,index)';
end

Btot_PR_map=squeeze(Btot_map_phi(phi_index,:,:));
BpolX_PR_map=squeeze(bX_map_phi(phi_index,:,:)).*Btot_PR_map;
BpolZ_PR_map=squeeze(bZ_map_phi(phi_index,:,:)).*Btot_PR_map;

EX_PR_map=squeeze(Efield_X_map_phi(phi_index,:,:));
EZ_PR_map=squeeze(Efield_Z_map_phi(phi_index,:,:));
grad_Phi_tor_PR_map=squeeze(grad_Phi_tor_map_phi(phi_index,:,:));
psi_star_dot_PR_map=squeeze(psi_star_dot_PR_map_phi(phi_index,:,:));

Ephi_PR_map=psi_star_dot_PR_map./Rpos_PR_map(:,1:size_r)-grad_Phi_tor_PR_map;

Ephi_PR_map_recalc=BpolX_PR_map.*EX_PR_map+BpolZ_PR_map.*EZ_PR_map;
Ephi_PR_map_recalc=-Ephi_PR_map_recalc./Bphi_PR_map;
Bphi_PR_map_recalc=BpolX_PR_map.*EX_PR_map+BpolZ_PR_map.*EZ_PR_map;
Bphi_PR_map_recalc=-Bphi_PR_map_recalc./Ephi_PR_map;

vE_X_PR_map=EZ_PR_map.*Bphi_PR_map-Ephi_PR_map.*BpolZ_PR_map;
vE_Z_PR_map=Ephi_PR_map.*BpolX_PR_map-EX_PR_map.*Bphi_PR_map;
vE_phi_PR_map=EX_PR_map.*BpolZ_PR_map-EZ_PR_map.*BpolX_PR_map;




%% now get XZ map for B
B_data=reshape(BpolX_PR_map,NP*size_r,1);
BX_XZ_map=gridfit(finesse_data_X,finesse_data_Z,B_data,scale_X,scale_Z,'smoothness',0.1);
BX_XZ_map=BX_XZ_map'; 

B_data=reshape(BpolZ_PR_map,NP*size_r,1);
BZ_XZ_map=gridfit(finesse_data_X,finesse_data_Z,B_data,scale_X,scale_Z,'smoothness',0.1);
BZ_XZ_map=BZ_XZ_map'; 

B_data=reshape(Btot_PR_map,NP*size_r,1);
Btot_XZ_map=gridfit(finesse_data_X,finesse_data_Z,B_data,scale_X,scale_Z,'smoothness',0.1);
Btot_XZ_map=Btot_XZ_map'; 


% calculate curvature
calculate_gradB_phi_rank;

BxgradB_X=-Bphi_XZsmall_map.*gradB_Z;
BxgradB_Z=Bphi_XZsmall_map.*gradB_X;
BxgradB_phi=BX_XZ_map.*gradB_Z-BZ_XZ_map.*gradB_X;

kappa_X_XZ_map=BxgradB_Z.*Bphi_XZsmall_map-BxgradB_phi.*BZ_XZ_map;
kappa_Z_XZ_map=BxgradB_phi.*BX_XZ_map-BxgradB_X.*Bphi_XZsmall_map;
kappa_phi_XZ_map=BxgradB_X.*BZ_XZ_map-BxgradB_Z.*BX_XZ_map;


%% now get XZ map for vE
vE_data=reshape(vE_X_PR_map(:,1:size_r),NP*size_r,1);
vE_X_XZ_map=gridfit(finesse_data_X,finesse_data_Z,vE_data,scale_X,scale_Z,'smoothness',0.1);
vE_X_XZ_map=vE_X_XZ_map'; 

vE_data=reshape(vE_Z_PR_map(:,1:size_r),NP*size_r,1);
vE_Z_XZ_map=gridfit(finesse_data_X,finesse_data_Z,vE_data,scale_X,scale_Z,'smoothness',0.1);
vE_Z_XZ_map=vE_Z_XZ_map'; 

vE_data=reshape(vE_phi_PR_map(:,1:size_r),NP*size_r,1);
vE_phi_XZ_map=gridfit(finesse_data_X,finesse_data_Z,vE_data,scale_X,scale_Z,'smoothness',0.1);
vE_phi_XZ_map=vE_phi_XZ_map'; 

vE_norm_XZ_map=sqrt(vE_X_XZ_map.^2+vE_Z_XZ_map.^2+vE_phi_XZ_map.^2);


xi_X_norm_XZ_map= vE_X_XZ_map./vE_norm_XZ_map;
xi_Z_norm_XZ_map= vE_Z_XZ_map./vE_norm_XZ_map;
xi_phi_norm_XZ_map= vE_phi_XZ_map./vE_norm_XZ_map;

%%
xi_X_norm_XZ_map=xi_X_norm_XZ_map.*psi_XZ_map_mask;
xi_Z_norm_XZ_map=xi_Z_norm_XZ_map.*psi_XZ_map_mask;
xi_phi_norm_XZ_map=xi_phi_norm_XZ_map.*psi_XZ_map_mask;

xi_norm_XZ_map=sqrt(xi_X_norm_XZ_map.^2+xi_Z_norm_XZ_map.^2+xi_phi_norm_XZ_map.^2);


imagesc(scale_X,scale_Z,xi_X_norm_XZ_map');


%%

% now finally evaluate deltaW from Phot !

