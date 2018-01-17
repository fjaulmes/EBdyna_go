function [const,maps,dim] = initialize_sim_data(par)
%initialize_sim_data Loads in data files for simulation
%   Returns structs with constants, B-maps and dimension parameters

%% Load in main files
filename=strcat(par.paths.DATA_FOLDER,'physics_constants.mat');
const=load(filename);
filename=strcat(par.paths.DATA_FOLDER,'XZsmall_fields_tokamak_pre_collapse.mat');
maps=load(filename,'size_X','size_Z','BpolX_initial_XZsmall_map','BpolZ_initial_XZsmall_map','Btot_XZ_map','psi_global','psi_XZsmall_map','psi_norm_XZsmall_map','theta_XZsmall_map');
filename=strcat(par.paths.DATA_FOLDER,'motions_map_dimensions.mat');
dim=load(filename,'mid_X','mid_Z','mid_Zzero','DX','scale_X','scale_Z','NB_PSI','R0');

%% Load PR-maps
filename=strcat(par.paths.DATA_FOLDER,'flux_geometry.mat');
map2=load(filename, 'dr_X_PR_map','dr_Z_PR_map');
maps=combine_structs(maps,map2);

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

%% Add RMP-field to data
if par.APPLY_RMP_FIELD
    %% Load RMP-field (theta,psi,phi!)
    n3D_path=strcat(par.paths.DATA_FOLDER,'BS_AUG_theta_psi_phi_n=2,2016-04-26.mat');
    if par.CALCULATE_PPHI_3D
        maps.n3D=load(n3D_path,'BR','BZ','Bphi','Aphi','dAR_dphi','dAZ_dphi','dAphi_dphi');
    else
        maps.n3D=load(n3D_path,'BR','BZ','Bphi');
    end
    
    %% Get psi-values in which RMP-map is made
    fnames=fieldnames(maps.n3D);
    dim.n3D.size_3D=size(maps.n3D.(fnames{1})); % original size 3D-map
    
    if dim.NB_PSI~=dim.n3D.size_3D(2)
        ratio_rmp_fin=(dim.NB_PSI-1)/(dim.n3D.size_3D(2)-1);
        if mod(ratio_rmp_fin,2)~=0
            error('Ratio of 3D field map and finesse not found. Used for psi interpolation');
        end
        dim.n3D.ind_fin_2_3D=@(ind) 1+(ind-1)/ratio_rmp_fin; % Function to find index in 3D-field map when one has the index in the finesse data
    else
        dim.n3D.psi_scale=dim.psi_scale;
    end
    theta_scale =linspace(0,2*pi,dim.n3D.size_3D(1));
    phi_scale   =linspace(0,2*pi,dim.n3D.size_3D(3));
    
    % Define incrediments
    dim.Dpsi  =1;
    dim.Dtheta=theta_scale(2)  -theta_scale(1);
    dim.Dphi  =phi_scale(2)    -phi_scale(1);
end

end