%% this script intiailizes the folder architecture for the EBdyna solver
restoredefaultpath

addpath('./d3d_184342_2500_15R_00')

% main folder for running the solver
addpath(genpath('./EBdyna/G_eq'))
% addpath('./G_eq/general_functions')
% addpath('./G_eq/init_data')
% addpath('./G_eq/trajectories')

% convert the data to an optimal structure for EBdyna
addpath('./EBdyna/init_functions')
init_physics_constants; % some convenient variables to have in context

% vizualize the input and output data
addpath('./EBdyna/visu_functions')
