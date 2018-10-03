volume_flux_psi=volume_flux*0;
volume_flux_psi(2:end)=volume_flux(2:end)-volume_flux(1:end-1);

alphas_psi=interp2(scale_X,scale_Z,psi_norm_XZsmall_map',alphas_pos_x,alphas_pos_z);



PSI_BIN_SIZE=radial_bin_size;
PSI_BIN_SIZE_HALF=round(0.5*PSI_BIN_SIZE);

P_alphas_profile=P_initial_profile-2*(Ne_profile.*Te_profile);


for psi_pos=1:N_radial_bins
    psi_pop=(alphas_psi>=(psi_pos-0.5)*PSI_BIN_SIZE).*(alphas_psi<(psi_pos+0.5)*PSI_BIN_SIZE);
    RES_POP=find(psi_pop);
    n_parts(psi_pos)=length(RES_POP);
    number_parts(psi_pos)=n_parts(psi_pos);
    T_parts(psi_pos)=(2/3)*mean(alphas_Ekin(RES_POP));
    volume_psi_zone=sum(volume_flux_psi((psi_pos-0.5)*PSI_BIN_SIZE+1:(psi_pos+0.5)*PSI_BIN_SIZE));
    n_parts(psi_pos)=n_parts(psi_pos)/volume_psi_zone;
    P_alphas_profile_rough(psi_pos)=mean(P_alphas_profile((psi_pos-0.5)*PSI_BIN_SIZE+1:(psi_pos+0.5)*PSI_BIN_SIZE));
end

n_alphas_profile=P_alphas_profile_rough./(T_parts*eV);
number_alphas_profile=n_alphas_profile.*volume_psi_zone;

density_part_ratio=mean(n_alphas_profile(1:20))/mean(n_parts(1:20));
particles_weight=mean(number_alphas_profile(1:20))/mean(number_parts(1:20))

n_parts=density_part_ratio*n_parts;

frac_fusion=mean(P_alphas_profile_rough)/mean(Ne_profile.*Te_profile)

%%
figure(4)
hold on;
plot(P_alphas_profile_rough,'b');
plot(n_parts.*T_parts*eV,'r');