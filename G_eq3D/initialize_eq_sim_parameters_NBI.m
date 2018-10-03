
tau_sim=2.5e-4;

TPRECISE=1.6;
TIME_STAMP_PRECISION=10;
RECORD_PRECISION=TIME_STAMP_PRECISION*40;
DELTA_TIME=(1e-9)/TPRECISE;
h=DELTA_TIME;
NB_TIME_STEPS=round(0.1*tau_sim/DELTA_TIME)
% no more than one time stamp in one loop
NB_TIME_STAMPS=round(NB_TIME_STEPS/RECORD_PRECISION)
% ten time steps in one loop
time_scale=(1:NB_TIME_STAMPS)*DELTA_TIME*RECORD_PRECISION*10;

%simulation options
SAVE_DATA_FILE=1;
USE_LAP_PSI=0;
% initial_NBI60keV_transp_D_distribution1
% initial_flatD_10keV_distribution
% PROCESS_VAL=(PROCESS_NUMBER+1)*5
%PROCESS_VAL=PROCESS_NUMBER
PROCESS_VAL=(round(PROCESS_NUMBER)-0.5)*10;

DISTNAME=strcat('initial_NBI60keV_transpM_D_distribution',num2str(PROCESS_NUMBER),'.mat')
SAVENAME=strcat('initial_NBI60keV_precession',num2str(PROCESS_NUMBER),'.mat')

SAVENAME_PRE=strcat('initial_NBI60keV_pre_collapse',num2str(PROCESS_NUMBER),'.mat')
SAVENAME_STATS=strcat('initial_NBI60keV_precession_stats',num2str(PROCESS_NUMBER),'.mat')


save(DISTNAME,'-append','PROCESS_NUMBER');

%to mimic collisional background
MIN_ALPHAS_EKIN=500


if SAVE_DATA_FILE==0
    disp('no data file for this simulation')
end

disp('no graphics will be displayed during this simulation')

