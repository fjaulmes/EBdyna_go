%%
clear all
SAVENAME='initial_DT_MB_distribution5.mat';
radial_bin_size=12;
Npart_min=11.6;
build_initial_DT_distribution_MB

clear all
SAVENAME='initial_DT_MB_distribution6.mat';
radial_bin_size=10;
Npart_min=8.6;
build_initial_DT_distribution_MB

clear all
SAVENAME='initial_DT_MB_distribution7.mat';
radial_bin_size=8;
Npart_min=11.6;
build_initial_DT_distribution_MB

clear all
SAVENAME='initial_DT_MB_distribution8.mat';
radial_bin_size=6;
Npart_min=11;
build_initial_DT_distribution_MB

%%

clear all
radial_bin_size=6;
energy_bin_size=61;
Npart_min=10.2;
SAVENAME='initial_alphas_vA_distribution5.mat';
build_initial_alpha_distribution_vA
n1=Nalphas_simulated
save n1

clear all
radial_bin_size=8;
energy_bin_size=57;
% Npart_min=round(0.012*N_energy_bins);
Npart_min=11.2;
SAVENAME='initial_alphas_vA_distribution6.mat';
build_initial_alpha_distribution_vA
n2=Nalphas_simulated
save n2


clear all
radial_bin_size=10;
energy_bin_size=55;
% Npart_min=round(0.013*N_energy_bins);
Npart_min=10.2;
SAVENAME='initial_alphas_vA_distribution7.mat';
build_initial_alpha_distribution_vA
n3=Nalphas_simulated
save n3

clear all
radial_bin_size=12;
energy_bin_size=53;
% Npart_min=round(0.014*N_energy_bins);
Npart_min=21;
SAVENAME='initial_alphas_vA_distribution8.mat';
build_initial_alpha_distribution_vA
n4=Nalphas_simulated
save n4
