function sim_parameters
%%initialize_eq_sim_parameters_struct Initializes the parameters (namely
%%time)
% Initializes time and savenames
global par
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Folder/file names

% Save names:
% Raw   - if failed, raw output
% Stats - Precession data
% Full  - final, evaluated result (after evaluate_output)

par.SAVENAME_RAW        =['./output/',par.ID_NAME, '_raw_process',num2str(par.PROCESS_NUMBER),'.mat'];      % Position data and intermediate storage
if ~par.USE_SPECIFIC_SAVENAME
   % use old way of naming prec file automatically
   if isfield(par,'SAVENAME_STATS') || isfield(par,'SAVENAME') 
	warning('Overriding SAVENAMEs with default naming convention!');
   end
   par.SAVENAME_STATS      =[par.folders.GEQ_OUTPUT,par.ID_NAME,'_prec_process',num2str(par.PROCESS_NUMBER),'.mat'];      % Precession data (e.g. trapped or passing)
   par.SAVENAME            =[par.folders.GEQ_OUTPUT,par.ID_NAME,'_full_process',num2str(par.PROCESS_NUMBER),'.mat'];      % Name of full data, with post-determined parameters (e.g. magnetic moment, toroidal angular momentum)
end

par.SAVENAME_RAW        =['./output/',par.ID_NAME, '_raw_process',num2str(par.PROCESS_NUMBER),'.mat'];      % Position data and intermediate storage

% the above makes little sense for 1 processor
% par.SAVENAME_STATS      =['./input/poincarre_stats.mat'];      % Precession data (e.g. trapped or passing)

loc_ID=strfind(par.ID_NAME,par.ID);
if isempty(loc_ID)
	error('ID is not contained within ID_NAME');
end
SAVENAME_STATS_xx   =[par.folders.GEQ_OUTPUT,par.ID_NAME(1:(loc_ID+1)),'xx_prec_process',num2str(par.PROCESS_NUMBER),'.mat'];            % Precession data from a previous run USES FIRST 2 OF THE ID-string
if ~par.TEST_DIST
    dist_fields = whos('-file',par.DISTNAME); dist_fields={dist_fields.name};
    ex_prec=any(strcmp(dist_fields,'prec'));
    
    ex_a=exist(par.SAVENAME_STATS,'file');
	if par.DO_NOT_USE_PREC
		ex_a = false;
		warning('previous prec file was ignored according to par.DO_NOT_USE_PREC value!')
	end
    ex_b=exist(SAVENAME_STATS_xx,'file');
    prec_data_avail= ex_prec | ex_a | ex_b;
    
    if ex_prec
        disp('PREC-info from distribution file used.')
    elseif ex_a && ex_b
        warning('Ignoring ID_xx-precession file since ID-precession file has been found')
    elseif ex_a && ~ex_b
        disp('ID-precession file found')
    elseif ex_b
        disp('ID_xx-precession file found')
    end
end
%% Determine mode (type of simulation)
if isfield(par,'mode')
    warning(['Mode has been enforced to: ',num2str(par.mode)])
elseif par.TEST_DIST
    par.mode=1;
elseif par.GET_PREC_INFO && ~prec_data_avail   % A precession file is needed, therefore mode 2
    par.mode=2;    
else
    par.mode=3;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulation time parameters
% General
par.NR_FUND_IN_LOOP = 5;  % Nr. fundamental steps in one loop

% Simulation dependent time parameters:
par.DISPLAY_TIMESTAMP_LOSSES = 1;
switch par.mode
    case 1 % TEST                      
        par.dt                      =0.5*1e-9;                          % Fundamental time step
        t_sim                       =20*1e-6;                          % Length of simulation
        par.NB_TIME_STAMPS          =10; 
        PC_SAVE                     =Inf;   
        par.NB_STAMPS_saved         =par.NB_TIME_STAMPS;
        
    case 2 % PRECESSION
        par.dt                      =0.5*1e-10;                        % Fundamental time step
        t_sim                       =8*1e-4;                          % Length of simulation
        par.NB_TIME_STAMPS          =800;                           % # stored values (stamps)
        PC_SAVE                     =Inf;                           % percentage after which a intermediate save is done (for debugging)
        par.NB_STAMPS_saved         =2;
        
        % Load the distribution
        par.LOADNAME                =par.DISTNAME;                  % Load in distribution file
        par.APPLY_3D=false;                                         % Turn off 3D simulation
        par.CALCULATE_PPHI_3D       =false;
        par.APPLY_SAWTOOTH          =false;                         % Tur off ST simulation
    case 3 % FULL                         
        par.dt                          =1*1e-9;                        % Fundamental time step
        t_sim                           =1*1e-6;                        % Length of simulation : 60 mus
        par.NB_TIME_STAMPS              =100;                           % # stored values (stamps)
        PC_SAVE                         =Inf;                           % percentage after which a intermediate save is done (for debugging)
        par.NB_STAMPS_saved             =100;        
        par.DISPLAY_TIMESTAMP_LOSSES    = 0;
   case 4
        error('Mode is obsolete')
    case 5 % POINCARE
        par.dt                      =2e-9;                          % Fundamental time step
        par.NR_FUND_IN_LOOP         =10;                            % # Fundamental time steps per loop
        t_sim                       =0.4;                           % Length of simulation
        par.NB_TIME_STAMPS          =1e3;                           % # stored values (stamps)
        PC_SAVE                     =10;        % percentage after which a intermediate save is done (for debugging)        
        
        par.poincare_type='tor';            % Plot in psi-theta (sfl) space or in R,Z (tor) space
    otherwise
        error('Failed to specify simulation time for this parameter set')
end




%% Determine distribution file
if ~par.TEST_DIST
    if ~par.GET_PREC_INFO || ~prec_data_avail || ex_prec
        par.LOADNAME=par.DISTNAME;
    elseif ex_a
        par.LOADNAME=par.SAVENAME_STATS;    % Load in the precession file (made previously in mode 2)
    elseif ex_b
        par.LOADNAME=SAVENAME_STATS_xx;    % Load in the precession file (made previously in mode 2)
    end
end
%% Correct time parameters so they match
time_per_loop               =par.NR_FUND_IN_LOOP*par.dt;        % fundamental time step in one loop / time step
par.NB_TIME_STEPS           =ceil(t_sim/time_per_loop);         % # loops
if par.NB_TIME_STEPS<par.NB_TIME_STAMPS
    warning('Truncating the number of time stamps, since the number of loops is less. Either change the par.NR_FUND_IN_LOOP or the simulation time')
    par.NB_TIME_STAMPS=par.NB_TIME_STEPS;
end

par.TIME_STAMP_PRECISION    =floor(par.NB_TIME_STEPS/par.NB_TIME_STAMPS);       % # loops before storing values
par.NB_TIME_STAMPS          =floor(par.NB_TIME_STEPS/par.TIME_STAMP_PRECISION); % Definite number of stamps (might be made higher due to longer simulation length

time_per_stamp              =time_per_loop*par.TIME_STAMP_PRECISION;    % time per stamp (time per loop times # loops before stamp)
par.time_scale              =time_per_stamp*(1:par.NB_TIME_STAMPS);     % time scale of stamps 

NB_SAVES                    =floor(100/PC_SAVE);                                % # of saves
par.RECORD_PRECISION        =floor(par.NB_TIME_STEPS/NB_SAVES);                 % # loops between 2 save points


if par.TIME_STAMP_PRECISION==0
    error('Cannot stamp more frequently than once per loop. Consider descreasing the number of stamps or the number of fundamental time steps in a loop.')
end
if par.NB_TIME_STAMPS<NB_SAVES && par.mode~=5
    error('First save has no output data and is therefore obsolete / error sensitive')
end
if ~isfield(par,'NB_STAMPS_saved')
    par.NB_STAMPS_saved         =par.NB_TIME_STAMPS;
end

if mod(par.NB_TIME_STAMPS,par.NB_STAMPS_saved)~=0
    par.NB_STAMPS_saved = par.NB_TIME_STAMPS/ceil(par.NB_TIME_STAMPS/par.NB_STAMPS_saved) ;
    warning('Altered NB_STAMPS_saved in order to match the original time stamping')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% locating maps file
filename=strcat(par.folders.DATA_SHOT,'XZsmall_fields_tokamak_pre_collapse.mat');
if ~exist(filename, 'file')
    warning('XZsmall_fields_tokamak_pre_collapse not found in expected location. Trying to find proper folder...')
    warning('folder not found')
    error('Manual error before recalculating B-field maps. Are you in the execution folder where this script is located?')
end

end