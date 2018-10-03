
load('initial_DT_MB_distribution5.mat', 'density_part_ratio')
density_part_ratio=0.25*density_part_ratio;
frac_fusion=0.007
density_part_ratio=frac_fusion*density_part_ratio

load('initialG_alphas_vA_all_pre_collapse.mat')

alphas_psi=interp2(scale_X,scale_Z,psi_norm_XZsmall_map',alphas_pos_x,alphas_pos_z);
volume_flux_psi(1)=0;
volume_flux_psi(2:Nradial)=volume_flux(2:end)-volume_flux(1:end-1);


clear Ekin_D nD volume_psi_zone


PSI_BIN_SIZE=10
PSI_BIN_SIZE_HALF=round(0.5*PSI_BIN_SIZE);

psi_bins_inf=(1:PSI_BIN_SIZE:Nradial-PSI_BIN_SIZE);
psi_bins_sup=(1+PSI_BIN_SIZE:PSI_BIN_SIZE:Nradial);
psi_bins_lim=(1:PSI_BIN_SIZE:Nradial);

psi_bins_pos=(1+PSI_BIN_SIZE_HALF:PSI_BIN_SIZE:Nradial-PSI_BIN_SIZE_HALF);

psi_bins=psi_scale(psi_bins_pos);
NB_PSI_BINS=length(psi_bins)-1
psi_bins_pos=psi_bins_pos(1:NB_PSI_BINS);


for psi_pos=1:NB_PSI_BINS
    psi_pop=(alphas_psi>=(psi_pos-0.5)*PSI_BIN_SIZE).*(alphas_psi<(psi_pos+0.5)*PSI_BIN_SIZE);
    RES_POP=find(psi_pop);
    n_alphas(psi_pos)=length(RES_POP);
    Ekin_alphas(psi_pos)=mean(alphas_Ekin(RES_POP));
    volume_psi_zone(psi_pos)=sum(volume_flux_psi((psi_pos-0.5)*PSI_BIN_SIZE+1:(psi_pos+0.5)*PSI_BIN_SIZE));

end
Ekin_alphas(NB_PSI_BINS)=0;
Ekin_alphas(NB_PSI_BINS-1)=0;

n_alphas_initial_profile=0*Ne_profile/density_part_ratio;
Ekin_alphas_initial_profile=0*Te_profile;
n_alphas=n_alphas./volume_psi_zone;
n_alphas_initial_profile(0.5*PSI_BIN_SIZE:(NB_PSI_BINS)*PSI_BIN_SIZE)=interp1(psi_bins_pos,n_alphas,0.5*PSI_BIN_SIZE:(NB_PSI_BINS)*PSI_BIN_SIZE);
Ekin_alphas_initial_profile=interp1(psi_bins_pos,Ekin_alphas,1:Nradial);
Talphas_initial_profile=(2/3)*Ekin_alphas_initial_profile;



% volume_flux_psi(1)=0;
% volume_flux_psi(2:Nradial)=volume_flux(2:end)-volume_flux(1:end-1);

% nD_initial_profile=nD_initial_profile./volume_flux_psi;
n_alphas_initial_profile=n_alphas_initial_profile*density_part_ratio;


close all

figure(1)
hold on
plot(Talphas_initial_profile,'b')
% plot(Te_profile/eV,'b')

figure(2)
hold on
plot(n_alphas_initial_profile,'b')
% plot(Ne_profile,'b')

n_alphas_initial_profile=n_alphas_initial_profile;
T_alphas_initial_profile=Talphas_initial_profile*eV;








load('alphas_vA_all_collapse_Glisa_fc0p8h2_G160414.mat')

alphas_psi=interp2(scale_X,scale_Z,psi_norm_XZsmall_map',alphas_pos_x,alphas_pos_z);




for psi_pos=1:NB_PSI_BINS
    psi_pop=(alphas_psi>=(psi_pos-0.5)*PSI_BIN_SIZE).*(alphas_psi<(psi_pos+0.5)*PSI_BIN_SIZE);
    RES_POP=find(psi_pop);
    n_alphas(psi_pos)=length(RES_POP);
    Ekin_alphas(psi_pos)=mean(alphas_Ekin(RES_POP));
end
Ekin_alphas(NB_PSI_BINS)=0;
Ekin_alphas(NB_PSI_BINS-1)=0;

n_alphas_final_profile=0*Ne_profile/density_part_ratio;
Ekin_alphas_final_profile=0*Te_profile;
n_alphas=n_alphas./volume_psi_zone;
n_alphas_final_profile(0.5*PSI_BIN_SIZE:(NB_PSI_BINS)*PSI_BIN_SIZE)=interp1(psi_bins_pos,n_alphas,0.5*PSI_BIN_SIZE:(NB_PSI_BINS)*PSI_BIN_SIZE);
Ekin_alphas_final_profile=interp1(psi_bins_pos,Ekin_alphas,1:Nradial);
T_alphas_final_profile=(2/3)*Ekin_alphas_final_profile;



% volume_flux_psi(1)=0;
% volume_flux_psi(2:Nradial)=volume_flux(2:end)-volume_flux(1:end-1);

% nD_initial_profile=nD_initial_profile./volume_flux_psi;
n_alphas_final_profile=n_alphas_final_profile*density_part_ratio;



figure(1)
hold on
plot(T_alphas_final_profile,'r')
% plot(Te_profile/eV,'b')

figure(2)
hold on
plot(n_alphas_final_profile,'r')
% plot(Ne_profile,'b')

T_alphas_final_profile=T_alphas_final_profile*eV;

save alphas_kinetic_density_profiles.mat n_alphas_initial_profile T_alphas_initial_profile n_alphas_final_profile T_alphas_final_profile