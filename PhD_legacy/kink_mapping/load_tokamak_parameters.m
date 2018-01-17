% loading previously calculated equilibrium data

filename=strcat(DATA_FOLDER,'physics_constants.mat')
load(filename);
filename=strcat(DATA_FOLDER,'tokamak_PR_map.mat')
load(filename);
filename=strcat(DATA_FOLDER,'B_fields.mat')
load(filename);
filename=strcat(DATA_FOLDER,'pressure_profile.mat')
load(filename);
filename=strcat(DATA_FOLDER,'psi_profiles.mat')
load(filename);
filename=strcat(DATA_FOLDER,'q_profile.mat')
load(filename);
filename=strcat(DATA_FOLDER,'flux_geometry.mat')
load(filename);
filename=strcat(DATA_FOLDER,'tokamak_parameters.mat')
load(filename);
filename=strcat(DATA_FOLDER,'tokamak_map_dimensions.mat')
load(filename);



Btor=mean(Btor_PR_map(1:NP-1,:),1);
Bpol_initial=mean(Bpol_PR_map(1:NP-1,:),1);
BHpol_profile=mean(BHpol_PR_map(1:NP-1,:),1);
q_initial=q_initial_profile;
radial_r_value=radial_r_value_flux;

