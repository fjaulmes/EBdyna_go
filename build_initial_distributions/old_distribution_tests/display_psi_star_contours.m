
psi_star_PR_map(:,:)=psi_star_PR_map_phi(phi_rank,:,:);

psi_data=reshape(E_potential_PR_map(:,:),NP*size_r,1);
psi_star_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,psi_data(IDTRI),XX_zoom,ZZ_zoom,DT_INTERPOLATION_METHOD);
psi_star_XZ_map=Epot_XZ_map';

imagesc(X_scale_zoom,Z_scale_zoom,psi_star_XZ_map',[-1000 1000]);
