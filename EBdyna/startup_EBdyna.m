%% this script intiailizes the folder architecture for the EBdyna solver
restoredefaultpath

% adds database acces
addpath('/sw/CDB/0.5-dev/src/matlab')
addpath('/compass/Shared/Common/IT/data_access/MDS_JET')

% main folder for running the solver
addpath(genpath('./G_eq'))
% addpath('./G_eq/general_functions')
% addpath('./G_eq/init_data')
% addpath('./G_eq/trajectories')

% convert the data to an optimal structure for EBdyna
addpath('./init_functions')
init_physics_constants; % some convenient variables to have in context

% vizualize the input and output data
addpath('./visu_functions')
