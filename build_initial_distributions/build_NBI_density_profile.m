initialize_folder_names;
close all;

filename=strcat(DATA_FOLDER,'physics_constants.mat');
load(filename);
filename=strcat(DATA_FOLDER,'pressure_profile.mat');
load(filename);
filename=strcat(DATA_FOLDER,'flux_geometry.mat');
load(filename);
filename=strcat('../METIS_profiles.mat');
load(filename);
filename=strcat(DATA_FOLDER,'q_profile.mat');
load(filename);

if ~exist('rho_tor_scale')
    tor_flux_scale=q_initial_profile*0;
    for n=2:size_r
        tor_flux_scale(n)=trapz(psi_scale(1:n)',q_initial_profile(1:n));
    end
    rho_tor_scale=sqrt(tor_flux_scale/max(tor_flux_scale));
end

mHe=mD;
ZHe=1

% approximate value of <T> for NBI particles
% used to find out approximate radial density profile of the beam
TNBI=50*1e3 % in eV


if exist('METISdata')
    PNBI_profile_radial_ini=interp1(METISdata.rho_tor_metis,METISdata.pnbi_profile,rho_tor_scale);
    NNBI_profile_radial_ini=PNBI_profile_radial_ini/(TNBI*eV);
    NDCORE=NNBI_profile_radial_ini(1)
    ROTATION_profile_radial_ini=interp1(METISdata.rho_tor_metis,METISdata.Omega_profile,rho_tor_scale);
    CORE_ROTATION=ROTATION_profile_radial_ini(1)  %rad/s
    MEAN_Q1_ROTATION=-mean(ROTATION_profile_radial_ini(1:size_r));
    NBI_density_profile=NNBI_profile_radial_ini;
else
    disp('METIS info is required to build the profiles')
end


plot(radial_r_value_flux,NBI_density_profile)

filename=strcat(DATA_FOLDER,'NBI_density_profile.mat');
save (filename,'NBI_density_profile','TNBI','NDCORE','mHe','ZHe','ROTATION_profile_radial_ini','CORE_ROTATION','MEAN_Q1_ROTATION');
