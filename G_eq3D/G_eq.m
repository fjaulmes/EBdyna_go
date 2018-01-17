function [exitcode] = G_eq(PROCESS_NUMBER,NB_PROCESS,load_worker)
%G_eq Main file for simulation
%   This file sets parameters and loads particles 
%   Fully functional for parallel computation or standalone.
%   Calls GT_eq, and thereby starts the simulation 
warning('on','all');

%% Global parameters
global maps dim par time

if isstruct(par) 
    par=remove_fields(par,{'mode'}); % Remove the mode-field if parameters already exist
end

%% Comment this section if you want to use same maps / parameters / dim
if isstruct(par) && isfield(par,'shared_memory') && par.shared_memory && ~isempty(maps)
    detach_shared_mem;
end
maps=[];
par=[];     % Clear any possible parameter
time=0;

% Please enter the distribution name
par.DISTNAME            ='./input/poincare_distribution_0.04eV_rand_edge.mat'; 

%% PROCESSORS
switch nargin
    case 0 % Local runs / test
        par.PROCESS_NUMBER=1;         % This 'box'
        par.NB_PROCESS=1;            % Number of 'boxes'
        par.TEST_DIST=false;            % Use the test distribution and parameters?       
        load_worker=false;              % Load maps in shared RAM
        par.shared_memory=false;        % Make use of maps in shared RAM
	case 2 % Deployment on LISA
        par.PROCESS_NUMBER=str2double(PROCESS_NUMBER);  % This 'box'
        par.NB_PROCESS=str2double(NB_PROCESS);          % Number of 'boxes'
        par.TEST_DIST=false;                            % Use the test distribution and parameters?
        load_worker=false;                              % Load maps in shared RAM
        par.shared_memory=true;                         % Make use of maps in shared RAM
    case 3 % Loading of RAM
        par.PROCESS_NUMBER=0;         % This 'box'
        par.NB_PROCESS=0;            % Number of 'boxes'
		par.TEST_DIST=false;                            % Use the test distribution and parameters?
        par.shared_memory=true;                        % Make use of maps in shared RAM
        if isempty(PROCESS_NUMBER) && isempty(NB_PROCESS) && load_worker
            disp('LOADING RAM WITH SHARED DATA!')
        elseif load_worker
            error('One cannot load shared maps and do a certain job (not the use case)')
        else
            warning('Continuing this box w/o loading a shared RAM. Consider using this function with only 2 input parameters')
            exitcode=G_eq(PROCESS_NUMBER,NB_PROCESS);   % Please consider not using the third parameter
            return
        end
    otherwise
        error('Input arguments don''t meet requirements')
end

%% Simulation Parameters
% Please enter the ID of the simulation:
% 3xxx - LISA / large with cubic        interpolation
% 2xxx - LISA / large with quadratic    interpolation
% 1xxx - LISA / large with linear       interpolation
% 0xxx      local run
%
% x1xx - trial      with particles distributed evenly in real space  (mid-plane LFS)
% x2xx - trial      with particles distributed evenly in q or psi and theta
% x3xx - poincare   with particles distributed random in psi and theta
% x4xx - NBI steady state distribution
% x5xx - homogeneous distribution
%
% xx00 - Test of the code
% xx01 - Nr. 01 of this simulation (last 2 digets are counting)

% Please provide a unique ID (see above table) and give the use case a comment (e.g. goal of the simulation)
par.ID='19xx';
par.comment='simulation in 19xx D edge test distribution';

% Please provide the interpolation scheme   (advice: 3)
% 0     -   self-written linear, with vector potential
% 1     -   self-written linear
% 2     -   interp linear               (ML standard)
% 3     -   ba_interp linear            (File Exchange)
% 4     -   griddedInterpolant linear   (makes use of lookup table)
% 5     -   interp cubic                (ML standard)
% 6     -   ba_interp cubic             (File Exchange)
% 7     -   griddedInterpolant cubic    (makes use of lookup table)
% 8     -   interp2 spline              (ML standard)
% 9     -   griddedInterpolant spline   (makes use of lookup table - RAM intensive!)
par.interp_scheme=3;

switch par.interp_scheme
	case {2,5,8}
		error('The interpolation scheme isn''t allowed, since it is too slow');
	case {9}
		error('The interpolation scheme isn''t allowed, since the RAM required is very high and cannot be in shared RAM');
end

% Please enter scheme for the particle pusher (FabienB for optimized BORIS, BORIS with vector-operations, "old" Fabien approximation)
% BORIS     - one-lined BORIS scheme without electric field  -- Better with interp schemes: 3 or 6
% FabienB   - BORIS scheme based on older Fabien code which might work faster sometimes -- Better with interp schemes: 1 2 4 5 7 8 9
switch par.interp_scheme
    case {3,6}
        par.scheme='BORIS';
    otherwise
        par.scheme='FabienB';
end

% Name to save files
par.ID_NAME=strcat('G_eq_',par.ID);     % Name to save files
par.coord_syst='flux';              % Switch for 3D coordinate system. Use 'flux': (theta,psi,phi) or 'toroidal': (R,Z,phi)

% 3D field
par.superimpose_2D_3D	=true; 					% Superimpose (if possible) the 2D and 3D fields, making a 2D interpolation obsolete and speeding up EBdyna_go. Disable to avoid numerical differences.

par.APPLY_RMP           =true;                  % Resonant Magnetic Perturbations (RMP)
if par.APPLY_RMP
    par.RMP_file        =['../data_tokamak/RMP_n=2_odd_flux_2018-01-07.mat']; 
end
par.APPLY_TFR           =false;                 % Toroidal Field Ripple
if par.APPLY_TFR
    par.TFR_file        =['../data_tokamak/','TFR_',par.coord_syst,'_2016-10-20.mat'];
end

if par.APPLY_RMP || par.APPLY_TFR
    par.APPLY_3D=true;
else
    par.APPLY_3D=false;
end
par.CALCULATE_PPHI_3D   =false...               % Determine pphi based on addition of (local) evolution
    & (par.APPLY_3D);

% Sawtooth
par.APPLY_SAWTOOTH      = false;
if par.APPLY_SAWTOOTH
    par.st.t_reconnection = 1e-4;
    par.st.t_relaxation   = 1.8e-4;
    par.st.n_reconnection = 1e3;
    par.st.n_relaxation   = 1e2;
    par.st.DR=1e-9;
    par.st.DZ=1e-9;
    par.st.stable_time    = 10; % The first time point in the sawtooth that should be taken into account! Condition: (k-kr)*c > -1 and no numeric oscillations
    par.st.rec_end_time   = 3; % Parameter from which the reconnection should get kappa_c from interpolation rather than formula.
    if ~any(par.interp_scheme==[3,6]); error('Sawtooth only compatible with ba_interp-interpolation'); end;
    par.st.show_parameters = false & nargin==0;
end

% Simulation type / distribution
par.GET_PREC_INFO       =true;                  % Find precession data in original 2D-field (trapped or not etc.)
par.SAVE_DATA_FILE      =true;                  % Save output files (raw + full + possible prec)




% Hard set the mode (after maybe precession in 2D has been made)
% [] - DEFAULT (first 2 then 3)
% 1  - TEST
% 2  - PRECESSION       Will be performed before 3,5 or 6
% 3  - FULL (after precession)
% 4  - FULL (without precession)
% 5  - POINCARE (need hard setting of the `mode')
% 6  - SAWTOOTH simulation
par.mode=5; % par.Poincare_plot=true; 
par.Poincare_plot=false;

%% Time paramters and savenames etc. sim_parameters;            
sim_parameters;            % Parameters (e.g. number of time steps and simulation length )

%% Load data
% Make maps coarser (use integers)
if par.interp_scheme==0
    par.step.R   =1;
    par.step.Z   =1;
    par.step.phi =1; % Take this one a 2^L number! (e.g. 1,2,4,8,16)
else
    par.step.R   =1;
    par.step.Z   =1;
    par.step.phi =1;
end
if mod(log(par.step.phi)/log(2),1)~=0; error('Phi-step not 2^L!'); end



if isempty(maps)
    if par.shared_memory
        initialize_maps(false,load_worker);               % Try to read shared memory or make it
        if load_worker; exitcode=55000; return; end
    else
        [maps,dim]            =initialize_maps(true,false);     % Calculate maps in own RAM / standalone session
    end
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
try     
% Initialize particles
[x,v,input,output,ejected]=initialize_particles;

% MAIN
% Start simulation
[exitcode]=GT_eq(x,v,input,output,ejected);

% Wrapup 
if par.shared_memory; 	detach_shared_mem;  end % 

catch err_initial
    %% Error catch
    if par.shared_memory; 	detach_shared_mem;  end
    try
        err_f_names=fieldnames(err_initial);
        for i=1:length(err_f_names)
            switch err_f_names{i}
                case 'message'
                    err_thrown.(err_f_names{i})=[err_initial.(err_f_names{i}),' IN PROCESS ',PROCESS_NUMBER];
                otherwise
                    err_thrown.(err_f_names{i})=err_initial.(err_f_names{i});
            end
        end
    catch
        rethrow(err_initial)
    end
    error(err_thrown)
end

end

%% Function to detach the shared RAM
function detach_shared_mem
global mach map_data maps
if ~isempty(maps)
    switch mach
        case 'WINDOWS'
            SharedMemory('detach',map_data.keys.key_maps ,maps);
        case 'LINUX'
            sharedmatrix('detach',map_data.keys.key_maps ,maps);
    end
end
clearvars -global const maps dim
end
