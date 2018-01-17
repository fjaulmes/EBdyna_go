
load('../data_tokamak/flux_geometry.mat')
alphas_psi=interp2(scale_X,scale_Z,psi_norm_XZsmall_map',alphas_pos_x,alphas_pos_z);

Te0=Te_profile(1)/eV

clear Ekin_DT nDT


PSI_BIN_SIZE=8
PSI_BIN_SIZE_HALF=round(0.5*PSI_BIN_SIZE);

psi_bins_inf=(1:PSI_BIN_SIZE:Nradial-PSI_BIN_SIZE);
psi_bins_sup=(1+PSI_BIN_SIZE:PSI_BIN_SIZE:Nradial);
psi_bins_lim=(1:PSI_BIN_SIZE:Nradial)

psi_bins_pos=(1+PSI_BIN_SIZE_HALF:PSI_BIN_SIZE:Nradial-PSI_BIN_SIZE_HALF)

psi_bins=psi_scale(psi_bins_pos);
NB_PSI_BINS=length(psi_bins)



for psi_pos=1:NB_PSI_BINS
    psi_pop=(alphas_psi>=(psi_pos-0.5)*PSI_BIN_SIZE).*(alphas_psi<(psi_pos+0.5)*PSI_BIN_SIZE);
    RES_POP=find(psi_pop);
    nDT(psi_pos)=length(RES_POP);
    Ekin_DT(psi_pos)=mean(alphas_Ekin(RES_POP));
    
end

nDT_initial_profile=interp1(psi_bins_pos,nDT,1:Nradial);
Ekin_DT_initial_profile=interp1(psi_bins_pos,Ekin_DT,1:Nradial);
Ti_DT_initial_profile=(2/3)*Ekin_DT_initial_profile;



volume_flux_psi(1)=0;
volume_flux_psi(2:Nradial)=volume_flux(2:end)-volume_flux(1:end-1);

nDT_initial_profile=nDT_initial_profile./volume_flux_psi;

P_DT_final_profile=nDT_initial_profile.*Ti_DT_initial_profile*eV;