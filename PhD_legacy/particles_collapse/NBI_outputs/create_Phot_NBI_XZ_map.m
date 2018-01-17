NBI_density_map=EKIN_BIN_SIZE*PICH_BIN_SIZE*squeeze(sum(sum(dist_EpitchRZ_ini(1:end,:,:,:),2),1));
NBI_density_map=(1e6)*NBI_density_map;

alphas_pos_R=pos_X_gc+R0;

dist_EkinRZ_ini=zeros(length(R_values),length(Z_values));

for rb=1:length(R_values)
    for zb=1:length(Z_values)
        PART_POP=find((alphas_pos_R>=R_BINS(rb)).*(alphas_pos_R<R_BINS(rb+1)).*(pos_Z_gc>=Z_BINS(zb)).*(pos_Z_gc<Z_BINS(zb+1)));
        dist_EkinRZ_ini(rb,zb) =mean(alphas_Ekin(PART_POP));  % cm3
    end
end

NBI_temp_map=(2/3)*dist_EkinRZ_ini;

NBI_Phot_map=(NBI_temp_map.*NBI_density_map*eV);

imagesc(R_values,Z_values,NBI_Phot_map');
colorbar
imagesc(R_values,Z_values,NBI_density_map');

% now fitting on the small ZX mesh

[XXm ZZm]=meshgrid(R_values-R0,Z_values)
P_data=reshape(NBI_Phot_map',size(NBI_Phot_map,1)*size(NBI_Phot_map,2),1);
NBI_Phot_XZ_map=gridfit(XXm,ZZm,P_data,scale_X,scale_Z,'smoothness',0.1);
NBI_Phot_XZ_map=NBI_Phot_XZ_map'; 


save NBIopp_simple_Phot_data.mat scale_X scale_Z NBI_Phot_XZ_map NBI_density_map

figure(2)
contour(scale_X,scale_Z,NBI_Phot_XZ_map',30);