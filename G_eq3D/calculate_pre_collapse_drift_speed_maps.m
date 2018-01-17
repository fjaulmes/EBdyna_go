function []=calculate_pre_collapse_drift_speed_maps(ISSQUAREMAP)
%% Load data
initialize_folder_names;    %Find the folders and load:
filename=strcat(DATA_FOLDER,'physics_constants.mat');
load(filename);
filename=strcat(DATA_FOLDER,'tokamak_map_dimensions.mat');
load(filename);
filename=strcat(DATA_FOLDER,'flux_geometry.mat');
load(filename);
filename=strcat(DATA_FOLDER,'tokamak_PR_map.mat');
load(filename);
filename=strcat(DATA_FOLDER,'B_fields.mat');
load(filename);
filename=strcat(DATA_FOLDER,'psi_star_evol.mat');
if exist(filename, 'file')
    load(filename);
end
filename=strcat(DATA_FOLDER,'q_profile.mat');
load(filename);

%% Other values
mid_X_large=mid_X;
mid_Z_large=mid_Z;
NB_PHI=513; % Number of toroidal positions incl. 0 and 2pi
DPHI=2*pi/(NB_PHI-1);   %Angle between toroidal positions
NB_PHI_DATA_HALF=round(0.5*(NB_PHI-1)); %Number of half the toroidal positions

% original mapping parameters

mid_Xaxis_large=mid_X_large+round(X_axis/DX);

%% Size of small map

SIZEX_LARGE=size(psi_XZ_map,1)
SIZEZ_LARGE=size(psi_XZ_map,2)

if ISSQUAREMAP
    % case of circular flux surfaces
    minX=40;
    maxX=SIZEX_LARGE-40;
    minZ=40;
    maxZ=SIZEZ_LARGE-40;
else
    % needs to be tuned for each equilibrium
    disp('dimensions of the small maps for 2D interpolation during simulations:')
    minX=max(mid_X_large-round(1.2*Nradial),3);
    maxX=min(mid_X_large+round(1.2*Nradial),SIZEX_LARGE-2);
    minZ=max(mid_Z-round(1.2*elongation*Nradial),3);
    maxZ=min(mid_Z+round(1.2*elongation*Nradial),SIZEZ_LARGE-2);
end


% Width and middle in index numbers
size_X=2*ceil(0.5*(maxX-minX));
size_Z=2*ceil(0.5*(maxZ-minZ));

mid_X=ceil(0.5*(maxX-minX))+1;
mid_Z=ceil(0.5*(maxZ-minZ))+1;

%% Centering the small map on the magnetic axis
X_pos_axis=round(X_axis/DX);
% taking half value for better centering
X_pos_center=0

Z_pos_center=0 % mandatory with current implementation
% later upgrade code to use mid_Zzero instead of mid_Z for index
% calculation

if ((maxX+X_pos_center)<SIZEX_LARGE) && ((maxZ+Z_pos_center)<SIZEZ_LARGE) && (ISSQUAREMAP==0)
    
    minX=minX+X_pos_center
    maxX=maxX+X_pos_center
    minZ=minZ+Z_pos_center
    maxZ=maxZ+Z_pos_center
    size_X=2*ceil(0.5*(maxX-minX));
    size_Z=2*ceil(0.5*(maxZ-minZ));

    mid_X=ceil(0.5*(maxX-minX))+1;
    mid_Xzero=mid_X-X_pos_center;
    mid_Z=ceil(0.5*(maxZ-minZ))+1;
    mid_Zzero=mid_Z-Z_pos_center;
    
    scale_X=DX*((1:size_X)-mid_X+X_pos_center);
    scale_Z=DX*((1:size_Z)-mid_Z+Z_pos_center);
    
elseif (ISSQUAREMAP==1)
    disp('using square map configuration')
    mid_X=ceil(0.5*(maxX-minX))+1;
    mid_Xzero=mid_X;
    scale_X=DX*((1:size_X)-mid_X);
    scale_Z=DX*((1:size_Z)-mid_Z);
else
    disp('problem with small maps dimensions : please edit script so that maps fit!')
    return
end
rescaling_to_XZsmall_maps;

Fmirror_coef=(eV/mHe);     % expressing mu in eV.T^-1

%% Making of maps
Btot_XZ_map=sqrt(BpolX_initial_XZsmall_map.^2+BpolZ_initial_XZsmall_map.^2+Bphi_XZsmall_map.^2);
% unit B-field
bX_XZ_map=BpolX_initial_XZsmall_map./Btot_XZ_map;
bZ_XZ_map=BpolZ_initial_XZsmall_map./Btot_XZ_map;
bphi_XZ_map=Bphi_XZsmall_map./Btot_XZ_map;

bX_XZ_map(isnan(bX_XZ_map))=0;
bZ_XZ_map(isnan(bZ_XZ_map))=0;
bphi_XZ_map(isnan(bphi_XZ_map))=0;

%% For vD investigation
gradB_X=zeros(size_X,size_Z);
gradB_Z=zeros(size_X,size_Z);

for (x=3:size_X-2)
    for (z=3:size_Z-2)
        if radial_XZsmall_map(x,z)<Nradial-1
            gradB_X(x,z)=(1/12)*(-Btot_XZ_map(x+2,z)+Btot_XZ_map(x-2,z))+(2/3)*(Btot_XZ_map(x+1,z)-Btot_XZ_map(x-1,z));
            gradB_Z(x,z)=(1/12)*(-Btot_XZ_map(x,z+2)+Btot_XZ_map(x,z-2))+(2/3)*(Btot_XZ_map(x,z+1)-Btot_XZ_map(x,z-1));
        end
    end
end
gradB_X(isnan(gradB_X))=0;
gradB_Z(isnan(gradB_Z))=0;

gradB_X=gradB_X/DX;
gradB_Z=gradB_Z/DX;

vD_X_XZ_map=-Bphi_XZsmall_map.*gradB_Z;
vD_Z_XZ_map=Bphi_XZsmall_map.*gradB_X;
vD_phi_XZ_map=BpolX_initial_XZsmall_map.*gradB_Z-BpolZ_initial_XZsmall_map.*gradB_X;

vD_X_XZ_map=vD_X_XZ_map./(Btot_XZ_map.^3);
vD_Z_XZ_map=vD_Z_XZ_map./(Btot_XZ_map.^3);
vD_phi_XZ_map=vD_phi_XZ_map./(Btot_XZ_map.^3);

gBX_B_XZ_map=gradB_X./(Btot_XZ_map);
gBZ_B_XZ_map=gradB_Z./(Btot_XZ_map);

vD_X_XZ_map(isnan(vD_X_XZ_map))=0;
vD_Z_XZ_map(isnan(vD_Z_XZ_map))=0;
vD_phi_XZ_map(isnan(vD_phi_XZ_map))=0;
gBX_B_XZ_map(isnan(gBX_B_XZ_map))=0;
gBZ_B_XZ_map(isnan(gBZ_B_XZ_map))=0;
%
%
% Fmirror_XZ_map=(gradB_X.*BpolX_initial_XZsmall_map+gradB_Z.*BpolZ_initial_XZsmall_map)./(Btot_XZ_map);
% Fmirror_XZ_map(isnan(Fmirror_XZ_map))=0;
%
% Fmirror_XZ_map=Fmirror_coef*Fmirror_XZ_map;

%Fmirror_coef=eV/mHe
% Fmirror_XZ_map=smooth_small_map(Fmirror_XZ_map);


%% Rename NP and NR from finesse
NB_PSI=Nradial;
NB_THETA=NP;
% NB_THETA=257;
DTHETA=2*pi/(NB_THETA-1);

% Don't know if this is still correct, since Nradial changed in more
% detailed FINESSE maps.
if ~exist('size_r')
    size_r=Nradial
end
scale_psi=1:size_r; % just the index number

scale_phi=2*pi*(0:NB_PHI-1)/(NB_PHI-1); % from 0 to 2pi
scale_theta=2*pi*(0:NB_THETA-1)/(NB_THETA-1); % from 0 to 2pi

% [scale_phi_3D scale_theta_3D scale_psi_3D] =meshgrid(scale_theta,scale_phi,1:size_r);
% psi_scale=(psi_scale-1).*psi_global;
% run('calculate_rotB_vDcurv')

% extending boundary to smooth things out at the edge
if psi_scale(1)<psi_scale(end) % increasing to 0
    MIN_PSI_VALUE=min(psi_scale)
    MAX_PSI_VALUE=max(psi_scale)
    psi_XZsmall_map(psi_XZsmall_map>-psi_scale(end-2))=-psi_scale(end-2); 
else                           % decreasing to 0
    MIN_PSI_VALUE=min(psi_scale)
    MAX_PSI_VALUE=max(psi_scale) %
    psi_XZsmall_map(psi_XZsmall_map<-psi_scale(end-2))=-psi_scale(end-2); 
end
% psi_XZsmall_map=max(psi_XZsmall_map,MIN_PSI_VALUE);

if max(psi_XZsmall_map(:))>MAX_PSI_VALUE  % to prevent extrapolation
    psi_norm_XZsmall_map=interp1(psi_scale,1:NB_PSI,min(psi_XZsmall_map,MAX_PSI_VALUE));
else
    psi_norm_XZsmall_map=interp1(psi_scale,1:NB_PSI,psi_XZsmall_map);
end
% Check for if a NaN occured by interpolation whilst psi was higher.
psi_norm_XZsmall_map(isnan(psi_norm_XZsmall_map))=1;

filename=strcat(DATA_FOLDER,'flux_geometry.mat');
load(filename, 'radial_r_value_flux')
r_XZsmall_map=interp1(psi_scale,radial_r_value_flux,psi_XZsmall_map);

%%
filename=strcat(DATA_FOLDER,'flux_geometry.mat');
load(filename, 'dr_X_PR_map');
load(filename, 'dr_Z_PR_map');
[XX ZZ]=meshgrid(scale_X,scale_Z);

FILENAME=strcat(FINESSE_FOLDER,'finesse_data.mat');
load(FILENAME);

FILENAME=strcat(FINESSE_FOLDER,'finesse_tokamak_parameters.mat');
load(FILENAME);

% get XZ maps for radial vector
dr_data=reshape((dr_X_PR_map'),NP*Nradial,1);

finesse_data_X=finesse_data(:,1)*a;
finesse_data_Z=finesse_data(:,2)*a;

dr_X=griddata(finesse_data_X,finesse_data_Z,dr_data,XX,ZZ,'cubic');
dr_X(isnan(dr_X))=0;
dr_X_XZ_map=dr_X';

dr_data=reshape((dr_Z_PR_map'),NP*Nradial,1);
dr_Z=griddata(finesse_data(:,1),finesse_data(:,2),dr_data,XX,ZZ,'cubic');
dr_Z(isnan(dr_Z))=0;
dr_Z_XZ_map=dr_Z';



%% Save maps and new names for dimensions
FILENAME=strcat(DATA_FOLDER,'XZsmall_fields_tokamak_pre_collapse.mat')
save (FILENAME,'q_initial_XZsmall_map','bX_XZ_map','bZ_XZ_map','bphi_XZ_map','size_X','size_Z','Btot_XZ_map','Bphi_XZsmall_map','BpolX_initial_XZsmall_map','BpolZ_initial_XZsmall_map','Rpos_XZsmall_map','radial_XZsmall_map','theta_XZsmall_map','psiH_XZsmall_map','psi_norm_XZsmall_map','psi_XZsmall_map','psi_global',...
 'vD_X_XZ_map','vD_Z_XZ_map','vD_phi_XZ_map','r_XZsmall_map','dr_X_XZ_map' ,'dr_Z_XZ_map'  );

Raxis=R0+X_axis;
FILENAME=strcat(DATA_FOLDER,'motions_map_dimensions.mat')
save (FILENAME,'psi_scale','scale_phi','scale_theta','scale_psi','mid_Xaxis_large','mid_Xzero','mid_Zzero','NB_PSI','NB_THETA','NB_PHI','DTHETA','DPHI','NB_PHI_DATA_HALF','DX','R0','a','XX_small','ZZ_small','Z_PR_map','scale_X','scale_Z','size_r','X_axis','Z_axis','mid_X','mid_Z','Raxis','Nradial');
