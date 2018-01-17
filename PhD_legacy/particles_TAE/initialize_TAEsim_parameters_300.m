
%Adjusting to correct the mode energy
WTAE_AVG=0.98*WTAE_AVG

%include a switch to remove the perturbation from the orbit trajectories
NO_PERTURB_ORBITS=0

gamma_TAE_min=-10*omega_TAE;
gamma_TAE_max=10*omega_TAE;

TAE_angle=2*pi/nTAE
NB_PHI
DPHI=TAE_angle/(NB_PHI-1);
size_r=size_r_TAE;
scale_phi=TAE_angle*(0:NB_PHI-1)/(NB_PHI-1);
scale_psi=pTAE_inf:pTAE_sup;

TAE_amplitude=1.0;
    

MINIMUM_TAE_AMPLITUDE=0.005

NB_OSCILLATIONS=8;

TPRECISE=4;
GYRO_ORBIT_PRECISION=20;
TIME_STAMP_PRECISION=1000*TPRECISE;

TIME_GO_SIZE=round(20*GYRO_ORBIT_PRECISION);

%we would like to have about 10 time steps between 2 frames
TAE_PERIOD=(2*pi)/omega_TAE
SIMULATION_TIME=TAE_PERIOD*NB_OSCILLATIONS
TIME_STEP_REF_SIZE=TAE_PERIOD/(NB_FRAME)/20/10

DELTA_TIME=TIME_STEP_REF_SIZE/TPRECISE;
h=DELTA_TIME;
NB_TIME_STEPS=round(0.05*SIMULATION_TIME/DELTA_TIME)
INTER_FRAME_TIME=TAE_PERIOD/NB_FRAME
% one time stamp in one loop
NB_TIME_STAMPS=round(NB_TIME_STEPS/TIME_STAMP_PRECISION)
NB_GYRO_STAMPS=round(NB_TIME_STEPS/GYRO_ORBIT_PRECISION)
TIME_STAMP_SIZE=round(20*NB_TIME_STEPS/NB_TIME_STAMPS)
% ten time steps in one loop
time_scale=(1:NB_TIME_STAMPS)*DELTA_TIME*TIME_STAMP_PRECISION*20;

TAE_ANGLE=(2*pi)/nTAE


% WTAE_AVG=mean(W_TAE_oscill_evol)

NB_PROCESS=64


LOAD_DATA_FILE=1;

if LOAD_DATA_FILE==0
    EKIN0_FAC=200
    SAVENAME=strcat('fewG_Ekin',num2str(EKIN0_FAC*8));
    SAVENAME=strcat(SAVENAME,'_TAEn8_250714_o');
    SAVENAME=strcat(SAVENAME,num2str(NB_OSCILLATIONS),'h');
    SAVENAME=strcat(SAVENAME,num2str(TPRECISE))	
else
    INPUTFILE='initial_MB_H200_pre_collapse_all.mat'
    SAVENAME=strcat('fewG_MB_H200_',num2str(PROCESS_NUMBER),'_TAEn',num2str(nTAE),'s_o');
    SAVENAME=strcat(SAVENAME,num2str(NB_OSCILLATIONS),'h');
    SAVENAME=strcat(SAVENAME,num2str(TPRECISE))
    SAVENAME=strcat(SAVENAME,'.mat')
    load('initial_MB_H200_precession_stats_all.mat','r_avg');
end

%NB_PART_RESCALE=particles_weight
NB_PART_RESCALE=4.22e12
%NB_PART_RESCALE=6e14


%simulation options
USE_VD=0
EVOLVE_TAE_AMPLITUDE=1
REINIT_PERP_SPEED=0;
SAVE_DATA_FILE=1;
DISPLAY_OUTPUTS=0;
CALCULATE_VD_POWER=1;
CALCULATE_DELTAE_POWER=1;

CALCULATE_TRUE_PPHI=1;
USE_DELTAE_TH=0
MIN_DIV_DPPHI=0.00001

PPHI_CORR_FACTOR=0.002