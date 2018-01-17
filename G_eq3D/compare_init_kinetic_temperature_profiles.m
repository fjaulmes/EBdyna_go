
close all

load('../data_tokamak/flux_geometry.mat')
volume_flux_psi(1)=0;
volume_flux_psi(2:Nradial)=volume_flux(2:end)-volume_flux(1:end-1);
load('../data_tokamak/volume_flux_geometry.mat')



alphas_psi=interp2(scale_X,scale_Z,psi_norm_XZsmall_map',alphas_pos_x,alphas_pos_z);

Te0=Te_profile(1)/eV

clear Ekin_D nD volume_psi_zone


PSI_BIN_SIZE=12
PSI_BIN_SIZE_HALF=round(0.5*PSI_BIN_SIZE);

psi_bins_inf=(1:PSI_BIN_SIZE:Nradial-PSI_BIN_SIZE);
psi_bins_sup=(1+PSI_BIN_SIZE:PSI_BIN_SIZE:Nradial);
psi_bins_lim=(1:PSI_BIN_SIZE:Nradial)

psi_bins_pos=(1+PSI_BIN_SIZE_HALF:PSI_BIN_SIZE:Nradial-PSI_BIN_SIZE_HALF)

psi_bins=psi_scale(psi_bins_pos);
NB_PSI_BINS=length(psi_bins)-1
psi_bins_pos=psi_bins_pos(1:NB_PSI_BINS)

% N_radial_bins=length(psi_bins_pos)


for psi_pos=1:NB_PSI_BINS
    psi_pop=(alphas_psi>=(psi_pos-0.5)*PSI_BIN_SIZE).*(alphas_psi<(psi_pos+0.5)*PSI_BIN_SIZE);
    RES_POP=find(psi_pop);
    nD(psi_pos)=length(RES_POP);
    Ekin_D(psi_pos)=mean(alphas_Ekin(RES_POP));
    volume_psi_zone(psi_pos)=sum(volume_flux_psi((psi_pos-0.5)*PSI_BIN_SIZE+1:(psi_pos+0.5)*PSI_BIN_SIZE));

end

nD=nD./volume_psi_zone;
density_part_ratio=Ne0/nD(1)

nD_initial_profile=Ne_profile/density_part_ratio;
nD_initial_profile(0.5*PSI_BIN_SIZE:(NB_PSI_BINS-9)*PSI_BIN_SIZE)=interp1(psi_bins_pos,nD,0.5*PSI_BIN_SIZE:(NB_PSI_BINS-9)*PSI_BIN_SIZE);
Ekin_D_initial_profile=interp1(psi_bins_pos,Ekin_D,1:Nradial);
Ti_D_initial_profile=(2/3)*Ekin_D_initial_profile;



% volume_flux_psi(1)=0;
% volume_flux_psi(2:Nradial)=volume_flux(2:end)-volume_flux(1:end-1);

% nD_initial_profile=nD_initial_profile./volume_flux_psi;
nD_initial_profile=nD_initial_profile*density_part_ratio;

P_D_final_profile=nD_initial_profile.*Ti_D_initial_profile*eV;

close all

figure(1)
hold on
plot(Ti_D_initial_profile,'r')
plot(Te_profile/eV,'b')

figure(2)
hold on
plot(nD_initial_profile,'r')
plot(Ne_profile,'b')

