%% Simulation time parameters
tau_sim=2.0e-3;

TIME_BASE=1e-9;             % 'Normal' time scale 1 ns
TPRECISE=2.0;              % Factor to decrease fundamental step
NR_FUND_IN_LOOP = 10;       % Nr. fundamental steps in one loop 

TIME_STAMP_PRECISION=100000;    % Number of loops for 1 timestamp
RECORD_PRECISION=TIME_STAMP_PRECISION*1; % Number of time steps before record  (or display)

DELTA_TIME=TIME_BASE/TPRECISE; %Fundamental time step
NB_TIME_STEPS=round(tau_sim/(NR_FUND_IN_LOOP*DELTA_TIME)) % NB time steps in simulation
NB_TIME_STAMPS=round(NB_TIME_STEPS/TIME_STAMP_PRECISION) % Number of time stamps 
time_scale=(1:NB_TIME_STAMPS)*NR_FUND_IN_LOOP*DELTA_TIME*TIME_STAMP_PRECISION; %time scale of stamps

% Rename DELTA_TIME for convenience
h=DELTA_TIME;

clear NR_FUND_IN_LOOP TIME_BASE

%% Save names
SAVE_DATA_FILE=1;
PROCESS_VAL=(round(PROCESS_NUMBER)-0.5)*10;

DISTNAME=strcat('initial_flatD20keV_pre_collapse_all.mat')

SAVENAME=strcat('post_RMP_flatD20keV_',num2str(PROCESS_NUMBER),'.mat')

%% Correction options
%to mimic collisional background
MIN_ALPHAS_EKIN=0

%algorithm correction
ECOEF=0
PPHI_CORR_FACTOR=0

if SAVE_DATA_FILE==0
    disp('no data file for this simulation')
end

disp('no graphics will be displayed during this simulation')

