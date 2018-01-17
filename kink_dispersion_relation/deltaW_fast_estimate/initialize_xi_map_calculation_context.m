format long;
format COMPACT
%warning off MATLAB:griddata:DuplicateDataPoints;

filename=strcat(DATA_FOLDER,'Epot_psi_star_dot_evol.mat');
load(filename);
filename=strcat(DATA_FOLDER,'tokamak_map_dimensions.mat');
load(filename);
filename=strcat(DATA_FOLDER,'psi_star_evol.mat');
load(filename);


radial_q1_size=psi_rank_q1+8


RR=Rpos_PR_map(:,1:size_r)-R0;

finesse_data_X=reshape((Rpos_PR_map(:,1:radial_q1_size)-R0),NP*radial_q1_size,1);
finesse_data_Z=reshape(Z_PR_map(:,1:radial_q1_size),NP*radial_q1_size,1);

NR=size_r+10
finesse_data_X_extended=reshape((Rpos_PR_map(:,1:NR)-R0),NP*NR,1);
finesse_data_Z_extended=reshape(Z_PR_map(:,1:NR),NP*NR,1);

PR_DT_INTERPOLATION_METHOD='linear';    % for the cast to PR maps
% getting the reference field in PR coords
BHpol_X_PR_map_size_r(:,:)=BHpol_X_PR_map(:,1:size_r);
BHpol_Z_PR_map_size_r(:,:)=BHpol_Z_PR_map(:,1:size_r);
BHpol_X_PR_map_size_r(:,1)=BHpol_X_PR_map_size_r(:,2);
BHpol_Z_PR_map_size_r(:,1)=BHpol_Z_PR_map_size_r(:,2);

finesse_mesh=[finesse_data_X  finesse_data_Z];
[finesse_mesh, IDTRI, J] = unique(finesse_mesh,'last','rows');
finesse_mesh_dtri=DelaunayTri(finesse_mesh);

finesse_mesh_extended=[finesse_data_X_extended  finesse_data_Z_extended];
[finesse_mesh_extended, IDTRI_EXT, J] = unique(finesse_mesh_extended,'last','rows');
finesse_mesh_extended_dtri=DelaunayTri(finesse_mesh_extended);

initialize_XZ_maps_dimensions;
[XX_small ZZ_small]=meshgrid(scale_X,scale_Z);
X_scale_data=reshape(XX_small,sizeX*sizeZ,1);
Z_scale_data=reshape(ZZ_small,sizeX*sizeZ,1);
XZ_mesh=[X_scale_data Z_scale_data];
XZ_mesh_dtri=DelaunayTri(XZ_mesh);

NB_THETA=NP;

PHI_OMEGA_RATIO=2;
NB_PHI=round((NP-1)/PHI_OMEGA_RATIO)+1;
DPHI=(2*pi/(NB_PHI-1));
DOMEGA=DPHI/PHI_OMEGA_RATIO;

CALCULATE_PR_DATA_FILE=1;
CALCULATE_BSTAR_DATA=0;
psi_max=size_r-1

