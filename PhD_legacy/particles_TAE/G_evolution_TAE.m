function [exitcode] = G_evolution_TAE(PROCESS_NUMBER_STRING)

REINIT_ALL_MAPS=0;
REINIT_SIMULATION_DATA=1;

if REINIT_ALL_MAPS==1
    clearvars -EXCEPT PROCESS_NUMBER_STRING
    initialize_folder_names_test;
    filename=strcat(DATA_FOLDER,'physics_constants.mat');
    load(filename);
    filename=strcat(DATA_FOLDER,'tokamak_map_dimensions.mat');
    load(filename);
    filename=strcat(DATA_FOLDER,'TAE_data.mat');
    load(filename);

    DPHI=TAE_angle/(NB_PHI-1);
    size_r=pTAE_sup-pTAE_inf+1;

    mid_X_large=mid_X;
    mid_Z_large=mid_Z;    
    calculate_TAE_drift_speed_maps;
    clear all;
    REINIT_SIMULATION_DATA=1;
end

if REINIT_SIMULATION_DATA==1
    clearvars -EXCEPT PROCESS_NUMBER_STRING
    clc
    format compact
    initialize_folder_names;
    
    filename=strcat(DATA_FOLDER,'physics_constants.mat');
    load(filename);
    filename=strcat(DATA_FOLDER,'XZsmall_fields_tokamak_TAE.mat');
    load(filename);
    filename=strcat(DATA_FOLDER,'motions_map_dimensions.mat');
    load(filename);
    %filename=strcat(DATA_FOLDER,'psi_profiles_kadomtsev.mat');
    %load(filename,'psi_star_initial');

    filename=strcat(DATA_FOLDER,'tokamak_PR_map.mat');
%     load(filename);
    load(filename, 'Rpos_PR_map');
    load(filename, 'radial_r_value_flux');
    load(filename,'BX_PR_map')
    load(filename, 'BZ_PR_map')
    filename=strcat(DATA_FOLDER,'B_fields.mat');
    load(filename,'Bpol_PR_map');
    load(filename,'Btor_PR_map');

    filename=strcat(DATA_FOLDER,'TAE_data.mat');
    load(filename);
 
    REINIT_SIMULATION_DATA=1;
   
end

PROCESS_NUMBER=str2num(PROCESS_NUMBER_STRING)

G_few_evolution_TAE_pll;


exitcode = 0;


