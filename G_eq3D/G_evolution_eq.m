function [exitcode] = G_evolution_eq(PROCESS_NUMBER_STRING)
%% Close all and clear
close all;
drawnow
clc
format compact
clearvars -EXCEPT PROCESS_NUMBER_STRING
%% Switches
% Recalculate maps/ load data
REINIT_ALL_MAPS=0; %ONLY SET TO 1 IN NOT PARALLEL
REINIT_SIMULATION_DATA=1;

% RMP-field
APPLY_RMP_FIELD=0
CALCULATE_TRUE_PPHI=0
%% Recalculate maps
if REINIT_ALL_MAPS==1
    calculate_pre_collapse_drift_speed_maps;
    REINIT_SIMULATION_DATA=1;
end

if REINIT_SIMULATION_DATA==1
    initialize_folder_names;
    
    filename=strcat(DATA_FOLDER,'physics_constants.mat');
    load(filename);
    filename=strcat(DATA_FOLDER,'XZsmall_fields_tokamak_pre_collapse.mat');
    load(filename);
    filename=strcat(DATA_FOLDER,'motions_map_dimensions.mat');
    load(filename);
%     filename=strcat(DATA_FOLDER,'psi_star_evol.mat');
%     load(filename, 'psi_star_2D_evol_interp');
%     psi_star_map_phi_evol=interp1(1:1001,psi_star_2D_evol_interp,1:101,'cubic');
    filename=strcat(DATA_FOLDER,'flux_geometry.mat');
    load(filename, 'dr_X_PR_map');
    load(filename, 'dr_Z_PR_map');
end

TOROIDAL_FIELD_DIRECTION=sign(mean(mean(Bphi_XZsmall_map)))

%% PROCESSORS
NB_PROCESS=32 % Number of processors
if nargin==0
    PROCESS_NUMBER_STRING='1'; % For test purposes
end
PROCESS_NUMBER=str2num(PROCESS_NUMBER_STRING)

%% MAIN
GT_many_evolution_eq_pll;


exitcode = 0;
