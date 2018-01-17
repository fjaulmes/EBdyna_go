function [ struct,data_tokamak ] = BS_AUG_flux_coordinates( struct)
%BS_AUG_flux_coordinates Loads the flux coordinates in a struct. Determines
%also R,Z,phi and X,Y,Z-coordinates. Note Z is relative to center of plasma!
global par
nr_psi=par.size_total(2);
nr_theta=par.size_total(1);
nr_phi=par.size_total(3);
indexes=par.indexes;
%% Input
narginchk(0,1);

if floor(log2(nr_psi-1))~=log2(nr_psi-1) || floor(log2(nr_theta-1))~=log2(nr_theta-1)
    error('Please provide nr_psi and nr_theta as 2^L+1 so it can exactly match FINESSE')
end

%% Load FINESSE files
load('../data_tokamak/motions_map_dimensions.mat', 'scale_X', 'scale_Z', 'R0')
load('../data_tokamak/flux_geometry.mat','Z_PR_map', 'X_PR_map','psi_scale');

z_axis=-0.031281806367984

%% Load Z-position shift from experimental data 
% (but the data file should be available for this process)
if ~exist('z_axis')
	%% Load Z-position from experimental data
	warning('Using data file to read Z-offset. Revise this section of the code!!!')
	filename='../data_tokamak/30382P01.CDF';
	if exist(filename,'file')
		z_axis=ncread(filename,'YAXIS');        % position magnetic axis (cm)
		T_data=ncread(filename,'TIME');         % time
		time_eq=2.5; % time equilibrium 
		[~,index_time]=min(abs(T_data-time_eq));
	%     index_time=108;
		z_axis=z_axis(index_time)*0.01;         % position magnetic axis (m)
	else
		error('File with Z-axis data not found in data_tokamak folder')
	%     Z_axis=0.074397392;
	end
end

%% Determine coordinates
NP=size(X_PR_map,1); %Number of poloidal/theta values from FINESSE
NR=size(X_PR_map,2); %Number of radial/psi values from FINESSE

if nr_psi>NR || nr_theta>NP
    error('More psi / flux points requested than found from data_tokamak (FINESSE)')
end

theta_indexes=round(linspace(1,NP,nr_theta)); %indexes of theta for BS3D
psi_indexes=round(linspace(1,NR,nr_psi)); %indexes of psi for BS3D

%% Flux coordinates (theta,psi,phi)
%phi-scale: toroidal angle
phi=linspace(0,2*pi,nr_phi)';

if nargin==0
    % theta
    theta_finesse=linspace(0,2*pi,NP);
    theta=theta_finesse(theta_indexes);  % theta's used in BS3D
    %Psi-scale, starts at 0.3 and ends at 0 [2*pi*Wb]
    psi=psi_scale(psi_indexes);
    
    % Grid:
    [~,~,struct.phi]=meshgrid(psi,theta,phi); %3D-mesh
    
    struct.R=R0+repmat(X_PR_map(theta_indexes,psi_indexes),[1 1 nr_phi]); %R-values in (theta,psi,phi)
    struct.Z=   repmat(Z_PR_map(theta_indexes,psi_indexes),[1 1 nr_phi]); %Z-values in (theta,psi,phi) (respect to magnetic axis)
else
    % Find theta and psi in indices that are used (e.g. 3:4 of 33)
    [theta_ind,psi_ind,phi_ind]=ind2sub([nr_theta nr_psi nr_phi],indexes); 
    % Find indeces in 2D matrix
    theta_ind_XZ=theta_indexes(theta_ind); 
    psi_ind_XZ=psi_indexes(psi_ind);
    % Get as single index
    in_2D_ind=sub2ind(size(X_PR_map),theta_ind_XZ,psi_ind_XZ);
    % Fill in to get R and Z
    struct.R=R0+X_PR_map(in_2D_ind)';
    struct.Z=   Z_PR_map(in_2D_ind)';
    struct.phi=phi(phi_ind);
end

%% Cartesian coordinates
struct.X= struct.R.*cos(struct.phi);     %x in Cartesian (m)
struct.Y=-struct.R.*sin(struct.phi);    %y in Cartesian (m)
%% Data of tokamak in struct
data_tokamak.z_axis=z_axis;
data_tokamak.size=[NP NR];
end