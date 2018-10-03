
close all;
pause(0.2);


clear all;
clc
format compact
initialize_folder_names;
DATA_FOLDER='../data_tokamak/'

filename=strcat(DATA_FOLDER,'physics_constants.mat');
load(filename);
filename=strcat(DATA_FOLDER,'XZsmall_fields_tokamak_pre_collapse.mat');
load(filename);
filename=strcat(DATA_FOLDER,'motions_map_dimensions.mat');
load(filename);
% filename=strcat(DATA_FOLDER,'psi_star_evol.mat');
% load(filename, 'psi_star_2D_evol_interp');
% psi_star_map_phi_evol=interp1(1:1001,psi_star_2D_evol_interp,1:101,'cubic');

filename=strcat(DATA_FOLDER,'flux_geometry.mat');
load(filename, 'dr_X_PR_map');
load(filename, 'dr_Z_PR_map');
filename=strcat(DATA_FOLDER,'tokamak_PR_map.mat');
load(filename);
filename=strcat(DATA_FOLDER,'B_fields.mat');
load(filename,'Btor_PR_map');
Bavg=mean(Btor_PR_map(:,1))

filename=strcat(DATA_FOLDER,'pressure_profile.mat');
load(filename);

REINIT_SIMULATION_DATA=1;

% psi_star_map_phi_evol=interp1(1:1001,psi_star_2D_evol_interp,1:101,'cubic');

% initialize completely the maps for pre-collapse calculations
FAST_SAWTOOTH=1;
tau_cr=4e-4;

load(strcat(DATA_FOLDER,'q_profile.mat'));
mid_X_axis=mid_X+round(X_axis/DX)
mid_Z_axis=mid_Z+round(Z_axis/DX)
%Bavg=Btot_XZ_map(mid_X_axis, mid_Z_axis)
