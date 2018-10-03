

filename=strcat(DATA_FOLDER,'physics_constants.mat');
load(filename);
filename=strcat(DATA_FOLDER,'psi_star_evol.mat');
load(filename);
filename=strcat(DATA_FOLDER,'tokamak_map_dimensions.mat');
load(filename);
filename=strcat(DATA_FOLDER,'tokamak_PR_map.mat');
load(filename);
filename=strcat(DATA_FOLDER,'B_fields.mat');
load(filename);
filename=strcat(DATA_FOLDER,'flux_geometry.mat');
load(filename);
filename=strcat(DATA_FOLDER,'psi_profiles_kadomtsev.mat');
load(filename);

increase_XZ_resolution;
interpolate_XZ_zoom_maps;


X_axis_offset=round(abs(X_axis)/DX_zoom)+1;
Z_axis_offset=round(abs(Z_axis)/DX_zoom)+1;

B3_XZ_zoom_map=Bphi_XZ_zoom_map./Rpos_XZ_zoom_map;
finesse_data_X=reshape((Rpos_PR_map(:,1:size_r)-R0),NP*size_r,1);
finesse_data_Z=reshape(Z_PR_map(:,1:size_r),NP*size_r,1);

initial_data=reshape(Bstar_PR_map(:,1:size_r),NP*size_r,1);

Bstar_initial_map=griddata(finesse_data_X,finesse_data_Z,initial_data,XX_zoom,ZZ_zoom,'cubic');
Bstar_initial_map=Bstar_initial_map';
Bstar_initial_map(isnan(Bstar_initial_map))=0;

Domega=(pi)/(Nomega-1);

% save Epot_calculation_context.mat