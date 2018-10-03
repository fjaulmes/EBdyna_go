
close all;
pause(0.2);
format compact

REINIT_ALL_MAPS=1;
REINIT_SIMULATION_DATA=1;


clear all
initialize_folder_names;
filename=strcat(DATA_FOLDER,'physics_constants.mat');
load(filename);
filename=strcat(DATA_FOLDER,'tokamak_map_dimensions.mat');
load(filename);

mid_X_large=mid_X;
mid_Z_large=mid_Z;    
NB_PHI=129;
DPHI=2*pi/(NB_PHI-1);
NB_PHI_DATA_HALF=round(0.5*(NB_PHI-1));
calculate_pre_collapse_drift_speed_maps(0);
