function [ struct,data_tokamak ] = BS_AUG_toroidal_coordinates( struct)
%BS_AUG_flux_coordinates Loads the flux coordinates in a struct. Determines
%also R,Z,phi and X,Y,Z-coordinates. Note Z is relative to center of plasma!
global par
% Either use nr_X and nr_Z as number of X and Z on map or as fractions on
% map. Please note fractions might result in incorrect error it doesn't fit
% precisely

indexes=par.indexes;
nr_X=par.size_total(1);
nr_Z=par.size_total(2);
nr_phi=par.size_total(3);
%% Input
narginchk(0,1);

%% Load FINESSE files
dim=load('../data_tokamak/motions_map_dimensions.mat', 'scale_X', 'scale_Z', 'R0');

%% Load Z-position from experimental data
warning('Using data file to read Z-offset. Revise this section of the code!!!')
filename='../data_tokamak/30382P01.CDF';
if exist(filename,'file')
    z_axis=ncread(filename,'YAXIS');        % position magnetic axis (cm)
    T_data=ncread(filename,'TIME');         % time
    time_eq=2.5;                            % time equilibrium 
    [~,index_time]=min(abs(T_data-time_eq));
%     index_time=108;
    z_axis=z_axis(index_time)*0.01;         % position magnetic axis (m)
else
    error('File with Z-axis data not found in data_tokamak folder')
%     z_axis=0.074397392;
end

%% Determine coordinates
NX=numel(dim.scale_X) ;%Number of x values from FINESSE
NZ=numel(dim.scale_Z) ; %Number of poloidal/theta values from FINESSE

if nr_X<1 && mod(NX,nr_X)==0
    nr_X=NX*nr_X;
end
if nr_Z<1 && mod(NZ,nr_Z)==0
    nr_Z=NZ*nr_Z;
end

if mod(NX,nr_X)~=0 || mod(NZ,nr_Z)~=0
    error('More X / Z points requested than found from data_tokamak (FINESSE)')
end

%% Toroidal coordinates (R,Z,phi)
phi=linspace(0,2*pi,nr_phi)';
% Grid:
r=dim.R0+linspace(dim.scale_X(1),dim.scale_X(end),nr_X)';
z=linspace(dim.scale_Z(1),dim.scale_Z(end),nr_Z)';
if nargin==0
    [struct.Z,struct.R,struct.phi]=meshgrid(z,r,phi); %3D-mesh
else
    [R_ind,Z_ind,phi_ind]=ind2sub([nr_X nr_Z nr_phi],indexes);
    struct.R=r(R_ind);
    struct.Z=z(Z_ind);
    struct.phi=phi(phi_ind);
end
%% Cartesian coordinates
struct.X= struct.R.*cos(struct.phi);     %x in Cartesian (m)
struct.Y=-struct.R.*sin(struct.phi);    %y in Cartesian (m)
%% Data of tokamak in struct
data_tokamak.z_axis=z_axis;
end