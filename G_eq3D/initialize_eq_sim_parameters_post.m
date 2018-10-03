
tau_sim=5.0e-4;

TPRECISE=2.0;
TIME_STAMP_PRECISION=10;
RECORD_PRECISION=TIME_STAMP_PRECISION*25;
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


%DISTNAME=strcat('initial_alphas_MB_D_distribution',num2str(PROCESS_NUMBER),'.mat')
DISTNAME=strcat('initial_NBI60keV_transp_D_distribution',num2str(PROCESS_NUMBER),'.mat')

%INPUTNAME=strcat('MB_D_fc2h1p6_',num2str(PROCESS_NUMBER),'.mat')
%SAVENAME=strcat('finalG_MB_D_precession',num2str(PROCESS_NUMBER),'.mat');
INPUTNAME=strcat('NBI60kev_fc1p6h2_',num2str(PROCESS_NUMBER),'.mat')
SAVENAME=strcat('finalG_NBI60kev_precession',num2str(PROCESS_NUMBER),'.mat');

save(DISTNAME,'-append','PROCESS_NUMBER');

if SAVE_DATA_FILE==0
    disp('no data file for this simulation')
end

disp('no graphics will be displayed during this simulation')

