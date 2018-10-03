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
if ~exist('rho_pol_scale')
    psi_norm_scale=1+(psi_scale/psi_global);
    rho_pol_scale=sqrt(psi_norm_scale);
end
% mW=183.84*1.660538921*1e-27
% mC=12.011*1.660538921*1e-27
mHe=mD;
ZHe=1


if exist('METISdata')
    N_profile_radial_ini=interp1(OMFITdata.rho_profile_ne,OMFITdata.ne_values,rho_pol_scale);
    ROTATION_profile_radial_ini=interp1(METISdata.rho_tor_metis,METISdata.Omega_profile,rho_tor_scale);
    CORE_ROTATION=ROTATION_profile_radial_ini(1)  %rad/s
    D_density_profile=N_profile_radial_ini;
    Ti_profile=interp1(OMFITdata.rho_profile_Te,OMFITdata.Ti_values,rho_pol_scale);
 else
    CORE_ROTATION=1.2e5;
    D_density_profile=(Ne_profile/Ne_profile(1));
end



filename=strcat(DATA_FOLDER,'D_profile.mat');
save (filename,'D_density_profile','Ti_profile','mHe','mD','ZHe','ROTATION_profile_radial_ini','CORE_ROTATION');
