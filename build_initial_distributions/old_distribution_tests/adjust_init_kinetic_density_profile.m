volume_flux_psi=volume_flux*0;
volume_flux_psi(2:end)=volume_flux(2:end)-volume_flux(1:end-1)

alphas_psi=interp2(scale_X,scale_Z,psi_norm_XZsmall_map',alphas_pos_x,alphas_pos_z);



PSI_BIN_SIZE=radial_bin_size;
PSI_BIN_SIZE_HALF=round(0.5*PSI_BIN_SIZE);



for psi_pos=1:N_radial_bins
    psi_pop=(alphas_psi>=(psi_pos-0.5)*PSI_BIN_SIZE).*(alphas_psi<(psi_pos+0.5)*PSI_BIN_SIZE);
    RES_POP=find(psi_pop);
    n_parts(psi_pos)=length(RES_POP);
    volume_psi_zone=sum(volume_flux_psi((psi_pos-0.5)*PSI_BIN_SIZE+1:(psi_pos+0.5)*PSI_BIN_SIZE))
    n_parts(psi_pos)=n_parts(psi_pos)/volume_psi_zone;

end

n_parts=density_part_ratio*n_parts;

density_part_ratio=density_part_ratio*Ne0*mean(Ne_psi_bins(1:end-2))/mean(n_parts(1:end-2))

% n_part_initial_profile=interp1(psi_bins_pos,n_parts,1:Nradial);
% 
% volume_flux_psi(1)=0;
% volume_flux_psi(2:Nradial)=volume_flux(2:end)-volume_flux(1:end-1);
% 
% n_part_initial_profile=n_part_initial_profile./volume_flux_psi;

