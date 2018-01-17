TOROIDAL_FIELD_DIRECTION=sign(mean(mean(Bphi_XZsmall_map)))
PSI_CORE_SIGN=sign(psi_scale(1))
PSI_STAR_SIGN=sign(psi_star_initial(round(0.5*size_r)))

% initialize completely the maps for pre-collapse calculations
FAST_SAWTOOTH=1.6;
tau_cr=4e-4;

TPRECISE=1.6;
TIME_STAMP_PRECISION=round(5*TPRECISE);
REINIT_PERP_SPEED=0;

DELTA_TIME=(1e-9)/TPRECISE;
h=DELTA_TIME;
NB_TIME_STEPS=round(0.05*(tau_cr/FAST_SAWTOOTH)/DELTA_TIME)
% one time stamp in one loop
NB_TIME_STAMPS=round(NB_TIME_STEPS/TIME_STAMP_PRECISION)
% ten time steps in one loop
time_scale=(1:NB_TIME_STAMPS)*DELTA_TIME*TIME_STAMP_PRECISION*20;

RECORD_TS=round(NB_TIME_STEPS/10)+1;

%simulation options
SAVE_DATA_FILE=1;
CALCULATE_TRUE_PPHI=0;

MIN_ALPHAS_EKIN=0;


%PROCESS_VAL=(round(PROCESS_NUMBER)-0.5)*10

INPUTFILE='initial_NBI60keV_R_pre_collapse_all.mat'
%INPUTFILE=strcat('initial_NBI60keV_pre_collapse',num2str(PROCESS_NUMBER),'.mat')
%INPUTFILE=strcat('initial_NBI60kev_',num2str(PROCESS_NUMBER),'pre_collapse.mat')
%INPUTFILE=strcat('initial_MB_W40c0_pre_collapse',num2str(PROCESS_NUMBER),'.mat')
%INPUTFILE=strcat('initial_flatD_',num2str(PROCESS_VAL),'keV_pre_collapse.mat')
SAVENAME=strcat('NBI60keV_R_fc1p6h1p6')
%SAVENAME=strcat('MB_W40c0_fc1h1')
%SAVENAME=strcat('MB_W40_f1h1p6_',num2str(PROCESS_NUMBER),'.mat')
%SAVENAME=strcat('flatD_',num2str(PROCESS_VAL),'keV_fc1p6h',num2str(TPRECISE),'.mat')


MEAN_Q1_ROTATION=-5e4



if SAVE_DATA_FILE==0
    disp('no data file for this simulation')
end

% if DISPLAY_OUTPUTS==0
    disp('no graphics will be displayed during this simulation')
% else
%    filename=strcat(DATA_FOLDER,'tokamak_PR_map.mat');
%    load(filename);
%    [XXsmall ZZsmall]=meshgrid(scale_X,scale_Z);
%    finesse_data_X=reshape((Rpos_PR_map(:,1:size_r)-R0),NP*size_r,1);
%    finesse_data_Z=reshape(Z_PR_map(:,1:size_r),NP*size_r,1);
% end
