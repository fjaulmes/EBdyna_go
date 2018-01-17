

inf_X=min(find(radial_XZ_map(1:mid_X,mid_Z)/Nradial-1));
sup_X=max(find(radial_XZ_map(mid_X:end,mid_Z)/Nradial-1))+mid_X-1;
inf_Z=1;
sup_Z=NZ;


FILENAME=strcat(DATA_FOLDER,'tokamak_map_dimensions.mat')
save (FILENAME,'NX','NZ','sup_X','sup_Z','inf_X','inf_Z','mid_X','mid_Z','DX','Rpos_map','Rpos','X_scale','Z_scale');
