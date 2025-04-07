function [ const,maps,dim] = reinitialize_data(par)
%reinitialize_data Loads in data files for simulation
%   Returns structs with constants, B-maps and dimension parameters
if ~isfield(par,'APPLY_RMP_FIELD')
    par.APPLY_RMP_FIELD=false;
end
%% Load in main files
filename=strcat(par.paths.DATA_FOLDER,'physics_constants.mat');
const=load(filename);
filename=strcat(par.paths.DATA_FOLDER,'XZsmall_fields_tokamak_pre_collapse.mat');
maps=load(filename);
filename=strcat(par.paths.DATA_FOLDER,'motions_map_dimensions.mat');
dim=load(filename);
%     filename=strcat(DATA_FOLDER,'psi_star_evol.mat');
%     load(filename, 'psi_star_2D_evol_interp');
%     psi_star_map_phi_evol=interp1(1:1001,psi_star_2D_evol_interp,1:101,'cubic');

%% Load PR-maps
filename=strcat(par.paths.DATA_FOLDER,'flux_geometry.mat');
map2=load(filename, 'dr_X_PR_map');
maps=combine_structs(maps,map2);
map2=load(filename, 'dr_Z_PR_map');
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

%% Add RMP-field to data
if par.APPLY_RMP_FIELD
    %% Load RMP-field
    RMP_path=strcat(par.paths.DATA_FOLDER,'BS_AUG_theta_psi_phi_n=2,2016-04-26.mat');
    RMP_field=load(RMP_path);
    if par.CALCULATE_PPHI
        RMP_field=remove_fields(RMP_field,{'AR','AZ'} );
    else
        RMP_field=remove_fields(RMP_field,{'AR','AZ','Aphi','dAR_dphi','dAZ_dphi','dAphi_dphi'} );
    end
    fnames=fieldnames(RMP_field);
    size_RMP=size(RMP_field.(fnames{1})); % original size RMP-map
	
    %% Load in correct struct
    for i=1:length(fnames)
        %permute to phi,theta,psi
        maps.RMP_field.(fnames{i})=permute(RMP_field.(fnames{i}),[3 1 2]); 
    end
    
    %% Get psi-values in which RMP-map is made
    if length(dim.scale_psi)~=size_RMP(2)
%         error('Check if psi data correspond with coordinates in RMP-map')
        ratio_rmpmap=round((length(dim.scale_psi)-1)/(size_RMP(2)-1));
        dim.RMPmap.psi_scale=dim.psi_scale(1:ratio_rmpmap:end);
    else
        dim.RMPmap.psi_scale=dim.psi_scale;
    end
    dim.RMPmap.psi_index=1:length(dim.RMPmap.psi_scale); % Norm is the one to interpolate in since psi in not equidistant
    dim.RMPmap.theta_scale=linspace(0,2*pi,size_RMP(1));
    dim.RMPmap.phi_scale=linspace(0,2*pi,size_RMP(3));
    
    % Define incrediments
    dim.Dpsi  =diff(dim.RMPmap.psi_index);      dim.Dpsi  =dim.Dpsi(1);
    dim.Dtheta=diff(dim.RMPmap.theta_scale);    dim.Dtheta=dim.Dtheta(1);
    dim.Dphi  =diff(dim.RMPmap.phi_scale);      dim.Dphi  =dim.Dphi(1);
end
%% Split two theta maps, to prevent interpolation error at 2*pi border
[maps]=split_theta_map(maps);
end