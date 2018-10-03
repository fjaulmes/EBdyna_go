

Epot_data=reshape(E_potential_PR_map(:,:),NP*size_r,1);
Epot_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,Epot_data(IDTRI),XX_zoom,ZZ_zoom,DT_INTERPOLATION_METHOD);
Epot_XZ_map=Epot_XZ_map';

imagesc(X_scale_zoom,Z_scale_zoom,Epot_XZ_map',[-1000 1000]);
