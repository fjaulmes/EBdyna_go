format long;
format COMPACT
%warning off MATLAB:griddata:DuplicateDataPoints;

filename=strcat(DATA_FOLDER,'physics_constants.mat');
load(filename);
% filename=strcat(DATA_FOLDER,'psi_star_evol.mat');
% load(filename);
filename=strcat(DATA_FOLDER,'tokamak_map_dimensions.mat');
load(filename);
filename=strcat(DATA_FOLDER,'tokamak_PR_map.mat');
load(filename);
filename=strcat(DATA_FOLDER,'B_fields.mat');
load(filename);
filename=strcat(DATA_FOLDER,'flux_geometry.mat');
load(filename);
filename=strcat(DATA_FOLDER,'psi_profiles.mat');
load(filename);
% filename=strcat(DATA_FOLDER,'psi_profiles_kadomtsev.mat');
% load(filename);
filename=strcat(DATA_FOLDER,'grad_flux_geometry.mat');
load(filename);
filename=strcat(DATA_FOLDER,'pressure_profile.mat');
load(filename);
filename=strcat(DATA_FOLDER,'q_profile.mat');
load(filename);

% define arbitrary background density
Ne0=2e19
Ne_profile=Ne0+Ne_profile*0;

%a priori no need to calculate the currents or the pressure
EVALUATE_CURRENTS=0
EVALUATE_PRESSURE=0

% finesse_mesh=[finesse_data_X  finesse_data_Z];
% [finesse_mesh, IDTRI, J] = unique(finesse_mesh,'last','rows');
% finesse_mesh_dtri=DelaunayTri(finesse_mesh);

% finesse_mesh_extended=[finesse_data_X_extended  finesse_data_Z_extended];
% [finesse_mesh_extended, IDTRI_EXT, J] = unique(finesse_mesh_extended,'last','rows');
% finesse_mesh_extended_dtri=DelaunayTri(finesse_mesh_extended);




