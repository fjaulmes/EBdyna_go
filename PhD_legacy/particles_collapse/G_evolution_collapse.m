function [exitcode] = G_evolution_collapse(PROCESS_NUMBER_STRING)



close all;

REINIT_ALL_MAPS=0;
REINIT_SIMULATION_DATA=1;

if REINIT_ALL_MAPS==1
    clearvars -EXCEPT PROCESS_NUMBER_STRING
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
    calculate_pre_collapse_drift_speed_maps;
    clearvars -EXCEPT PROCESS_NUMBER_STRING
    REINIT_SIMULATION_DATA=1;
end

if REINIT_SIMULATION_DATA==1
    clearvars -EXCEPT PROCESS_NUMBER_STRING
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
    load(filename);
    clear psi_star_2D_evol_interp;  % no need for this high precision
    
    filename=strcat(DATA_FOLDER,'flux_geometry.mat');
    load(filename, 'dr_X_PR_map');
    load(filename, 'dr_Z_PR_map');

    filename=strcat(DATA_FOLDER,'Epot_psi_star_dot_evol.mat');
    load(filename,'Epot_evol');
 
    filename=strcat(DATA_FOLDER,'psi_profiles_kadomtsev.mat');
    load(filename, 'psi_star_initial');
 
    REINIT_SIMULATION_DATA=1;
end

NB_PROCESS=64

PROCESS_NUMBER=str2num(PROCESS_NUMBER_STRING)

GT_many_evolution_collapse_pll;


exitcode = 0;


