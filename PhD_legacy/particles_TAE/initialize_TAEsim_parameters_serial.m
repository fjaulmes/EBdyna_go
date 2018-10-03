
gamma_TAE_min=-10*omega_TAE;
gamma_TAE_max=10*omega_TAE;

TAE_angle=2*pi/nTAE
NB_PHI
DPHI=TAE_angle/(NB_PHI-1);
size_r=size_r_TAE;
scale_phi=TAE_angle*(0:NB_PHI-1)/(NB_PHI-1);
scale_psi=pTAE_inf:pTAE_sup;
TAE_amplitude=50;
% TAE_amplitude=0;
  
NB_PART_RESCALE=1e15
MINIMUM_TAE_AMPLITUDE=0.005





% WTAE_AVG=mean(W_TAE_oscill_evol)

NB_PROCESS=1


LOAD_DATA_FILE=1;

% if LOAD_DATA_FILE==0
%     EKIN0_FAC=200
%     SAVENAME=strcat('fewG_Ekin',num2str(EKIN0_FAC*8));
%     SAVENAME=strcat(SAVENAME,'_TAEn8_250714_o');
%     SAVENAME=strcat(SAVENAME,num2str(NB_OSCILLATIONS),'h');
%     SAVENAME=strcat(SAVENAME,num2str(TPRECISE))	
% else
%     INPUTFILE='initialG_alphas_TAE_co_n12.mat'
%     SAVENAME=strcat('fewG_alphas_co',num2str(PROCESS_NUMBER),'_TAEn',num2str(nTAE),'s_o');
%     SAVENAME=strcat(SAVENAME,num2str(NB_OSCILLATIONS),'h');
%     SAVENAME=strcat(SAVENAME,num2str(TPRECISE))
% end


%simulation options
USE_DELTAE_TH=1
MIN_DIV_DPPHI=0.002


nTAE

TOROIDAL_FIELD_DIRECTION=sign(mean(mean(Bphi_XZsmall_map)))
PSI_CORE_SIGN=sign(psi_scale(1))
% PSI_STAR_SIGN=sign(psi_star_initial(round(0.5*size_r)))

% psi_star_map_phi_evol=interp1(1:1001,psi_star_2D_evol_interp,1:101,'cubic');

% initialize completely the maps for pre-collapse calculations
NB_OSCILLATIONS=8;

TPRECISE=2;
GYRO_ORBIT_PRECISION=10*TPRECISE;
TIME_STAMP_PRECISION=100*TPRECISE;
TIME_GO_SIZE=round(20*GYRO_ORBIT_PRECISION)

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


%simulation options
USE_VD=0
EVOLVE_TAE_AMPLITUDE=1;
REINIT_PERP_SPEED=0;
LOAD_DATA_FILE=1;
SAVE_DATA_FILE=1;
DISPLAY_OUTPUTS=0;
CALCULATE_TRUE_PPHI=1
CALCULATE_VD_POWER=1;
CALCULATE_DELTAE_POWER=1;

if (USE_VD==1)&&(CALCULATE_VD_POWER==0)
    disp('inconsitent simulation parameters for vD !')
    CALCULATE_VD_POWER=1
end

if LOAD_DATA_FILE==0
    SAVENAME=strcat('fewG_Ekin',num2str(EKIN0_FAC*8));
    SAVENAME=strcat(SAVENAME,'_TAEn8_250714_o');
    SAVENAME=strcat(SAVENAME,num2str(NB_OSCILLATIONS),'h');
    SAVENAME=strcat(SAVENAME,num2str(TPRECISE))
else
    INPUTFILE='initialG_alphas_TAE_co_n8m9.mat'
    SAVENAME=strcat('fewG_alphas_coh_TAEn',num2str(nTAE),'s_o');
    SAVENAME=strcat(SAVENAME,num2str(NB_OSCILLATIONS),'h');
    SAVENAME=strcat(SAVENAME,num2str(TPRECISE))
end
