function [const,maps,dim]=load_distr_maps(par)
%% Load in main files
filename=strcat(par.paths.DATA_FOLDER,'physics_constants.mat');
const=load(filename);
filename=strcat(par.paths.DATA_FOLDER,'XZsmall_fields_tokamak_pre_collapse.mat');
maps=load(filename,'size_X','size_Z','Bphi_XZsmall_map','BpolX_initial_XZsmall_map','BpolZ_initial_XZsmall_map','Btot_XZ_map','psi_global','psi_XZsmall_map','psi_norm_XZsmall_map','theta_XZsmall_map');
filename=strcat(par.paths.DATA_FOLDER,'motions_map_dimensions.mat');
dim=load(filename,'a','mid_X','mid_Z','DX','scale_X','scale_Z','NB_PSI','R0','Z_axis','mid_Xzero','X_axis','psi_scale','Raxis');
filename=strcat(par.paths.DATA_FOLDER,'q_profile.mat');
dim2=load(filename,'q_initial_profile'); 
dim=combine_structs(dim,dim2);
filename=strcat(par.paths.DATA_FOLDER,'flux_geometry.mat');
maps2=load(filename,'X_PR_map','Z_PR_map');
maps=combine_structs(maps,maps2);
%% Remove map name in fields, which I find unnessecary
fnames=fieldnames(maps);
for i=1:length(fnames)
    if strendswith(fnames{i},'small_map')
        maps.(fnames{i}(1:end-9))=maps.(fnames{i});
        maps=rmfield(maps,fnames{i});
    elseif strendswith(fnames{i},'_map')
        maps.(fnames{i}(1:end-4))=maps.(fnames{i});
        maps=rmfield(maps,fnames{i});
    end
end
%% Split two theta maps, to prevent interpolation error at 2*pi border
[maps]=split_theta_map(maps);
end