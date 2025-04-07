function [ ] = BS_COMPASS_TFR_main(PROCESS_NUMBER,N_PROCESS)
%%BS_COMPASS_TFR_main Determines the magnetic field of TFR field
% MAIN FILE
%   Script to determine field of coils in AUG. Can store field of each coil seperately in RMP_map-struct.
%   calls BS_COMPASS_RMP_load_coils to load in all coils from data in folder
%   define the number of psi, theta and phi points
%   Uses a standard current of 4.8kAt
global par
%#ok<*COLND>
%% PROCESS NUMBERS
switch nargin
    case 0
        par.PROCESS_NUMBER=1;
        par.N_PROCESS=1;
    case 1
        par.N_PROCESS=20;
        par.PROCESS_NUMBER=str2double(PROCESS_NUMBER);
    case 2
        par.N_PROCESS=str2double(N_PROCESS);
        par.PROCESS_NUMBER=str2double(PROCESS_NUMBER);
end

%% Parameters
TF_map.I=(5/5)*0.1995e6;   par.I=TF_map.I;       % Coil current [At] (default),
par.determine_vector_potential=true;    % Determines derivatives in phi of vector potential. Needs auxillary points for nice finite difference calculation.
par.example=0;                          % No plotting / example
par.coord_sys='toroidal';                   % Switch between coordinal systems
par.nr_coils=1;                         % #coils calculated : Switch between 1 or 16 coils
par.apply_Z_corection_from_METIS=0;     % take into account a shift from METIS data if FIESTA / FINESSE eq was calculated with a shit of center for Fourier transform of  separatrix
par.INPUTFILE='./input/CU_coil_3DXYZ_coords.mat';
par.REDUCE_DATA=0;

%% Load Coils
par.NB_TF_COILS=16
[TF_coils] = BS_COMPASS_TF_load_coils(par.INPUTFILE,par.NB_TF_COILS,par.REDUCE_DATA);
TF_coils=TF_coils(1:par.nr_coils);

%% Define grid where to determine B

[TF_map,Z_correction]=BS_COMPASS_get_grid(TF_map);
if par.apply_Z_corection_from_METIS==0
    Z_correction=0;
else
    disp('shifiting the Z coordinates in order to have proper Z=0 position in mid-plane (based on ZMA_METIS - ZMA_FINESSE)')
    Z_correction
end

%% Calculate the fields of each coil
[TF_map] = BS_TFR_calc_coils_individual(par,TF_map,TF_coils,Z_correction);

% Remove the link between A and B-positions / fields
remove_link={'index_link'};
par=remove_fields(par,remove_link);

%% Determine phi and R-components fields from original Cartesian
remove_obsolete={'BX','BY','AX','AY'}; % To remove obsolete fields
for i=1:length(TF_map.TF_coil)  % For each coil
    TF_map.TF_coil(i).Aphi  =-TF_map.X./TF_map.R.*TF_map.TF_coil(i).AY      +   TF_map.Y./TF_map.R.*TF_map.TF_coil(i).AX;
    TF_map.TF_coil(i).Bphi  =-TF_map.X./TF_map.R.*TF_map.TF_coil(i).BY      +   TF_map.Y./TF_map.R.*TF_map.TF_coil(i).BX;
    
    TF_map.TF_coil(i).AR    = TF_map.X./TF_map.R.*TF_map.TF_coil(i).AX      +   TF_map.Y./TF_map.R.*TF_map.TF_coil(i).AY;
    TF_map.TF_coil(i).BR    = TF_map.X./TF_map.R.*TF_map.TF_coil(i).BX      +   TF_map.Y./TF_map.R.*TF_map.TF_coil(i).BY;
    
    % Store neccesary fields in seperate variable and pre-allocate
    if i==1
        field(length(TF_map.TF_coil))=remove_fields(TF_map.TF_coil(i),remove_obsolete); %#ok<*AGROW>
    end
    field(i)=remove_fields(TF_map.TF_coil(i),remove_obsolete); 
end

%% Save output
if ~par.example
    save_name_individual=strcat('./output/BS_COMPASS_TF3D_',par.coord_sys,'_individual_',datestr(now,'yyyy-mm-dd'),'_process',num2str(par.PROCESS_NUMBER),'.mat');
    save(save_name_individual,'-v7.3','field','par');
    disp(['Saved individual coil file: ',num2str(par.PROCESS_NUMBER),' of ',num2str(par.N_PROCESS)])
    return
end
end

%% Function to create 2,2-subplot with certain name, X,Z-plane
function [ha1,ha2,ha3,ha4]=make_figure(name)
hf=figure('name',name,'tag','BS_COMPASS_RMP_map_RMP_coils_figure');
ha1=subplot(2,2,1,'parent',hf); axis(ha1,'equal')
ha2=subplot(2,2,2,'parent',hf); axis(ha2,'equal')
ha3=subplot(2,2,3,'parent',hf); axis(ha3,'equal')
ha4=subplot(2,2,4,'parent',hf); axis(ha4,'equal')
hold(ha1,'on'), hold(ha2,'on'),hold(ha3,'on'),hold(ha4,'on')
for ha=[ha1 ha2 ha3 ha4]
    xlabel(ha,'$x$ [m]','interpreter','latex')
    ylabel(ha,'$z$ [m]','interpreter','latex')
end
end

%% Function to add coil shape in figure
function plot_coils(ha1,TF)

for i=1:16
    h=plot3(ha1,TF(i).X,TF(i).Y,TF(i).Z,'k','displayname',['TF coil ',num2str(i)]);
    if i~=1
        hasbehavior(h,'legend',false)
    end
end
end