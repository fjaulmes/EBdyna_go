function [exitcode] = G_eq(PROCESS_NUMBER,NB_PROCESS)
%G_eq Main file for simulation
%   This file sets parameters and loads particles 
%   Fully functional for parallel computation or standalone.
%   Calls GT_eq, and thereby starts the simulation 
warning('on','all');



%% Global parameters
global maps dim par time

% IT4I: run on the IT4I cluster
% IPP.CR: run on the IPP CAS servers
par.mach='IPP.CR';


%% Comment this section if you want to use same maps / parameters / dim

maps=[];
par=[];     % Clear any possible parameter
time=0;
par.shot_no=184342;
par.tokamak='d3d';
par.time_id='2500'
par.nbi_id='15R'
par.sim_id='00'
par.shot_name=[num2str(par.shot_no)];
par.ID=[par.shot_name,'_',par.time_id,'_',par.nbi_id,'_',par.sim_id];
par.sim_folder_string=['d3d_',par.ID];

par.mach='IT4I';
folders=initialize_folder_names();

% input and output file names
par.DO_NOT_USE_PREC      = false;
par.USE_SPECIFIC_SAVENAME= true;
par.SAVENAME_STATS       = [folders.GEQ_OUTPUT,'stats_NUR_markers_D3D.mat'];
par.SAVENAME             = [folders.GEQ_OUTPUT,'full_NUR_markers_D3D.mat'];
par.DISTNAME             = [folders.GEQ_INPUT,'184342_2500_test_FIDASIM_markers.mat'];
par.VESSEL_FILENAME      = [folders.DATA_COMMON_TOKAMAK,'RZ_wall_mask.mat'];
par.USE_VESSEL_LIMIT     = 1;
par.remove_thermalized_markers = 0;   % not a slowing down simulation
par.TEST_DIST=false;                   % Use the test distribution and parameters?     

% Hard set the mode (after maybe precession in 2D has been made)
% [] - DEFAULT (first 2 then 3)
% 1  - TEST
% 2  - PRECESSION       Will be performed before 3,5 or 6
% 3  - FULL (after precession)
% 4  - FULL (without precession)
% 5  - POINCARE (need hard setting of the `mode')
% 6  - SAWTOOTH simulation
par.mode=3;


define_basic_parameters();


par.TESTING_CX=0;
par.TESTING_SD=0;

try
    dist=load(par.DISTNAME);
catch
    error('THE DISTRIBUTION FILE COULD NOT BE LOADED !!!')
end

if isfield(dist,'input') % refers to an already formatted input dist data from a previous simulation
    input=dist.input;
    if isfield(input,'SD_MARKERS_END')
        par.SD_MARKERS_END       = input.SD_MARKERS_END;
        disp(['Updating SD_MARKERS_END to value = ' num2str(par.SD_MARKERS_END)] );
        disp('the simulation will split sim between sd markers and ionization markers');
    end
else
    %         disp('all markers treated equally....')
end

    
%% PROCESSORS
switch nargin
    case 0 % Local runs / test
        par.PROCESS_NUMBER=1;         % This 'box'
        par.NB_PROCESS=1;            % Number of 'boxes'
        load_worker=false;              % Load maps in shared RAM
        par.shared_memory=false;        % Make use of maps in shared RAM
	case 2 
		if strcmp(par.mach,'IT4I')
			par.PROCESS_NUMBER=str2double(PROCESS_NUMBER);  % This 'box'
			par.NB_PROCESS=str2double(NB_PROCESS);          % Number of 'boxes'
			load_worker=false;                              % No shared RAM
			par.shared_memory=false;                         % No shared RAM
        end
    otherwise
        error('Input arguments don''t meet requirements')
end

if par.NB_PROCESS>1
    % standard names are needed for parallel runs
	par.USE_SPECIFIC_SAVENAME= false; 
end

   %% Time paramters and savenames etc. sim_parameters;            
    sim_parameters;            % Parameters (e.g. number of time steps and simulation length )

    %% Load data

    if isempty(maps)
        [maps,dim]            =initialize_maps_ram();     % Calculate maps in own RAM / standalone session
    else
        warning('Re-using maps / data from previous run!')
    end

    par.offset_midplane=[0.0, dim.Z_axis]; % Polynomials of midplane (R) where first is slope, second is offset (ax+b)


% Add function for offset midplane (function Z=Z(R))
dim.func_Z_cross=@(R) polyval(par.offset_midplane,R);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  DO NOT MESS WITH WHATS BELOW %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start the script as interpolations, if an error occurs, detach any shared RAM

  
% Initialize particles
[x,v,input,output,ejected]=initialize_particles;

% Force turn off of centrifugal effect
% because of little knowledge of rotation profile
par.APPLY_FC=0;
input.Fc_field=x(:,1)*0;

% MAIN
% Start simulation
[exitcode]=GT_eq(x,v,input,output,ejected);



