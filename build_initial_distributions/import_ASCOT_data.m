%% Import data from text file.
% Script for importing data from the following text file:
%
%    C:\Users\fabien\EBdyna_go\JET_sim\build_initial_distributions\ascot_ion_distrib.dat
%
% To extend the code to different selected data or a different text file,
% generate a function instead of a script.

% Auto-generated by MATLAB on 2017/10/12 11:59:38

%% Initialize variables.
filename = 'C:\Users\fabien\EBdyna_go\JET_sim\build_initial_distributions\ascot_ion_distrib.dat';
delimiter = ' ';

%% Format string for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN, 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
clear   alphas_Ekin alphas_pos_x  alphas_pos_z alphas_pitch alphas_weight
format compact

load('../data_tokamak/physics_constants.mat')
load('../data_tokamak/motions_map_dimensions.mat');
load('../data_tokamak/XZsmall_fields_tokamak_pre_collapse.mat','psi_XZsmall_map');

ZHe=1
mHe=mD

MIN_EKIN=10000
alphas_Ekin = dataArray{:, 3} /eV ; % eV
alphas_pos_x = dataArray{:, 1};
alphas_pos_z = dataArray{:, 2};

PARTPOP=find((alphas_Ekin>MIN_EKIN));
alphas_Ekin=alphas_Ekin(PARTPOP);

alphas_pos_x = dataArray{:, 1};
alphas_pos_z = dataArray{:, 2};
alphas_pitch = dataArray{:, 4};
alphas_weight = dataArray{:, 5};

alphas_pos_x=alphas_pos_x(PARTPOP);
alphas_pos_z=alphas_pos_z(PARTPOP);
alphas_pitch=alphas_pitch(PARTPOP);
alphas_weight=alphas_weight(PARTPOP);


ZAXIS_POS_ASCOT=0.5*(max(alphas_pos_z)+min(alphas_pos_z))
disp('shifting vertically the distribution to have it match the position of axis in FINESSE equilibrium')
alphas_pos_z=-alphas_pos_z+ZAXIS_POS_ASCOT-Z_axis;
ZMIN_POS=scale_Z(4)
ZMAX_POS=scale_Z(end-3)
alphas_pos_z(alphas_pos_z<ZMIN_POS)=alphas_pos_z(alphas_pos_z<ZMIN_POS)*0.95;
alphas_pos_z(alphas_pos_z>ZMAX_POS)=alphas_pos_z(alphas_pos_z>ZMAX_POS)*0.95;

XAXIS_POS_ASCOT=0.5*(max(alphas_pos_x)+min(alphas_pos_x))
disp('shifting horizontally the distribution to have it match the position of axis in FINESSE equilibrium')
alphas_pos_x=alphas_pos_x-XAXIS_POS_ASCOT+R0;

psi_part=interp2(scale_X,scale_Z,psi_XZsmall_map',alphas_pos_x-R0,alphas_pos_z);
PARTPOP=find(psi_part<-0.1);
alphas_Ekin=alphas_Ekin(PARTPOP);
alphas_pos_x=alphas_pos_x(PARTPOP);
alphas_pos_z=alphas_pos_z(PARTPOP);
alphas_pitch=alphas_pitch(PARTPOP);
alphas_weight=alphas_weight(PARTPOP);


Nalphas_simulated=length(alphas_pos_x);

% Deuterium particles
alphas_vtot=sqrt(2*alphas_Ekin*(eV/mD));
alphas_vpll=alphas_pitch.*alphas_vtot;

alphas_pos_phi=rand(Nalphas_simulated,1)*2*pi;

load('../data_tokamak/tokamak_PR_map.mat');
alphas_pos_x=alphas_pos_x-R0;


%% Allocate an integer number of times the particles so that we increase toroidal representation
NUMBER_TOR_MUL=10

if NUMBER_TOR_MUL>1
    
    alphas_pos_x=repmat(alphas_pos_x,NUMBER_TOR_MUL,1);
    alphas_pos_z=repmat(alphas_pos_z,NUMBER_TOR_MUL,1);
    alphas_weight=repmat(alphas_weight,NUMBER_TOR_MUL,1);
    alphas_Ekin=repmat(alphas_Ekin,NUMBER_TOR_MUL,1);
    alphas_pitch=repmat(alphas_pitch,NUMBER_TOR_MUL,1);
    
    alphas_weight=alphas_weight/NUMBER_TOR_MUL;
    
    % add a little random part to all data to avoid artificial pixelized
    % features
    
    DX=0.01;
    DZ=0.015;
    DEkin=100;
    Dpitch=0.02;
    
    Nalphas_simulated=length(alphas_pos_x)
    
    alphas_pos_x=alphas_pos_x+(rand(Nalphas_simulated,1)-0.5)*DX;
    alphas_pos_z=alphas_pos_z+(rand(Nalphas_simulated,1)-0.5)*DZ;
    alphas_Ekin=alphas_Ekin+(rand(Nalphas_simulated,1)-0.5)*DEkin;
    alphas_Ekin(alphas_Ekin<0)
    alphas_pitch=alphas_pitch+(rand(Nalphas_simulated,1)-0.5)*Dpitch;
    alphas_pitch=max(alphas_pitch,-1);
    alphas_pitch=min(alphas_pitch,1);
   
    % recalculate everything
    alphas_vtot=sqrt(2*alphas_Ekin*(eV/mD));
    alphas_vpll=alphas_pitch.*alphas_vtot;
    
    
    alphas_pos_phi=rand(Nalphas_simulated,1)*2*pi;
    
end

%% use description as for structures in EBdyna_go

SAVE_STRUCTS=0

if SAVE_STRUCTS==1
    x=zeros(length(alphas_pos_x),3);
    x(:,1)=alphas_pos_x;
    x(:,2)=alphas_pos_z;
    x(:,3)=alphas_pos_phi;
    
    
    vperp=sqrt(alphas_vtot.^2-alphas_vpll.^2);
    % 1.1.4.2.2. Determine v by direction of b
    [N,b]=make_normal_vector(x,Nalphas_simulated); % Vector of perpendicular velocity
    v=bsxfun(@times,vpll,b)+bsxfun(@times,vperp,N);
    
end
%% Clear temporary variables
clearvars filename delimiter formatSpec fileID dataArray ans;
SAVENAME='ASCOT_dist_initial_data.mat'

save(SAVENAME,'Nalphas_simulated','mHe','ZHe','alphas_pos_x','alphas_pos_z','alphas_pos_phi','alphas_vpll','alphas_Ekin','alphas_weight')


%%

 contour(scale_X,scale_Z,psi_XZsmall_map',[0 0])
hold on
plot(alphas_pos_x,alphas_pos_z,'b.')
