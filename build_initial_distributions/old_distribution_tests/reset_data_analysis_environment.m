
close all;
pause(0.2);


clear all;
clc
format compact
initialize_folder_names;

filename=strcat(DATA_FOLDER,'physics_constants.mat');
load(filename);
filename=strcat(DATA_FOLDER,'XZsmall_fields_tokamak_pre_collapse.mat');
load(filename);
filename=strcat(DATA_FOLDER,'motions_map_dimensions.mat');
load(filename);
filename=strcat(DATA_FOLDER,'psi_star_evol.mat');
if exist(filename,'file')
    load(filename);
end
% psi_star_map_phi_evol=interp1(1:1001,psi_star_2D_evol_interp,1:101,'cubic');
filename=strcat(DATA_FOLDER,'B_fields.mat');
load(filename,'Btor_PR_map');
Bavg=mean(Btor_PR_map(:,1))

filename=strcat(DATA_FOLDER,'flux_geometry.mat');
load(filename, 'dr_X_PR_map');
load(filename, 'dr_Z_PR_map');

filename=strcat(DATA_FOLDER,'pressure_profile.mat');
load(filename);

REINIT_SIMULATION_DATA=1;

% psi_star_map_phi_evol=interp1(1:1001,psi_star_2D_evol_interp,1:101,'cubic');

% initialize completely the maps for pre-collapse calculations
FAST_SAWTOOTH=1;
tau_cr=4e-4;

