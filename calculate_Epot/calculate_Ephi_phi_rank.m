phi_rank=11;
NB_PHI_DATA=65;
DPHI=(2*pi)/(NB_PHI-1);

%evaluate grad phi

E_potential_PR_map=zeros(NP,size_r);
E_potential_PR_map_prev=zeros(NP,size_r);
E_potential_PR_map_next=zeros(NP,size_r);
E_potential_PR_map(:,:)=E_potential_PR_map_phi(phi_rank,:,:);
E_potential_PR_map_prev(:,:)=E_potential_PR_map_phi(phi_rank-1,:,:);
E_potential_PR_map_next(:,:)=E_potential_PR_map_phi(phi_rank+1,:,:);
E_potential_PR_map_DPHI=E_potential_PR_map_next-E_potential_PR_map_prev;
E_potential_PR_data=reshape(E_potential_PR_map_DPHI(:,:),NP*size_r,1);
%E_potential_XZ_map_DPHI=griddata(finesse_data_X,finesse_data_Z,E_potential_PR_data,XX_zoom,ZZ_zoom,'cubic');
E_potential_XZ_map_DPHI=dtinterp(finesse_mesh,finesse_mesh_dtri,E_potential_PR_data(IDTRI),XX_zoom,ZZ_zoom,DT_INTERPOLATION_METHOD);
E_potential_XZ_map_DPHI=E_potential_XZ_map_DPHI';   
E_potential_XZ_map_DPHI(isnan(E_potential_XZ_map_DPHI))=0;

grad_Phi_tor_num=(E_potential_XZ_map_DPHI./Rpos_XZ_zoom_map)/(2*DPHI);



%psi star dot at mid point

psi_star_dot_PR_map=zeros(NP,size_r);
psi_star_dot_PR_map(:,:)=psi_star_dot_PR_map_phi(phi_rank,:,:);

psi_star_dot_data=reshape(psi_star_dot_PR_map(:,1:size_r),NP*size_r,1);
%psi_star_dot_XZ_zoom_map=griddata(finesse_data_X,finesse_data_Z,psi_star_dot_data,XX_zoom,ZZ_zoom,'cubic');
psi_star_dot_XZ_zoom_map=dtinterp(finesse_mesh,finesse_mesh_dtri,psi_star_dot_data(IDTRI),XX_zoom,ZZ_zoom,DT_INTERPOLATION_METHOD);

psi_star_dot_XZ_zoom_map=psi_star_dot_XZ_zoom_map';
psi_star_dot_XZ_zoom_map(isnan(psi_star_dot_XZ_zoom_map))=0;
psi_star_dot_XZ_zoom_map=smooth_small_map(psi_star_dot_XZ_zoom_map);

%psi star at mid point

psi_star_PR_map=zeros(NP,size_r);
psi_star_PR_map(:,:)=psi_star_PR_map_phi(phi_rank,:,:);

psi_star_data=reshape(psi_star_PR_map(:,1:size_r),NP*size_r,1);
%psi_star_XZ_zoom_map=griddata(finesse_data_X,finesse_data_Z,psi_star_data,XX_zoom,ZZ_zoom,'cubic');
psi_star_XZ_zoom_map=dtinterp(finesse_mesh,finesse_mesh_dtri,psi_star_data(IDTRI),XX_zoom,ZZ_zoom,DT_INTERPOLATION_METHOD);

psi_star_XZ_zoom_map=psi_star_XZ_zoom_map';
psi_star_XZ_zoom_map(isnan(psi_star_XZ_zoom_map))=0;
psi_star_XZ_zoom_map=smooth_small_map(psi_star_XZ_zoom_map);


% compute back Ephi

Ephi_XZ_zoom_map_num=psi_star_dot_XZ_zoom_map./Rpos_XZ_zoom_map-grad_Phi_tor_num;

close all
run('calculate_Bstar_RZ');
figure(2);
imagesc(X_scale_zoom,Z_scale_zoom,Ephi_XZ_zoom_map_num');
axis xy
brighten(0.4);colorbar;
phi=(1/(NB_PHI-1))*(phi_rank-1)*(2*pi);
run('define_reference_potential_axis');
