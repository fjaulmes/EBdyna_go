
SAVE_LOG_FILE=0;
SAVE_DATA_FILE=1;
DT_INTERPOLATION_METHOD='quadratic'     % by default
PR_DT_INTERPOLATION_METHOD='linear';    % for the final cast to PR maps
INTEGRATION_PRECISION=1;                % 1 is high, 0 is low (two out of one dl function)
DISPLAY_OUTPUTS=1;                      % remainig problem with contourc function: still has to display contours....
ADDITIONAL_MAP_SMOOTHING=0;
CATCH_BACK_POTENTIAL_INTEGRATION_TO_ZERO=1;



[rmix_value EXPULSION_FRAME_NUMBER]=max(rx_evol_lin);
% EXPULSION_FRAME_NUMBER=round(EXPULSION_FRAME_NUMBER-1)
% for kink simulation
if EXPULSION_FRAME_NUMBER<10
    EXPULSION_FRAME_NUMBER=100
end



if SAVE_DATA_FILE==0
    disp('no data file for this simulation')
end

if DISPLAY_OUTPUTS==0
    disp('no graphics will be displayed during this simulation')
end

disp('**************************************************************');

SMALL_FRAME_NUMBER=25;


NB_PHI=129;
NB_PHI_RANK=17;
PHI_STEP_SIZE=round((NB_PHI-1)/(NB_PHI_RANK-1))
% NB_PHI_DATA_HALF=round(0.5*NB_PHI-1);

if SAVE_LOG_FILE==1
    disp('log information in script_output.txt');
    % open output file with write permission
    fid = fopen('script_output.txt', 'w');
    fprintf(fid, '****************** NB_PHI = %d *********************\n',NB_PHI);
    fprintf(fid, '********************************************************\n');
else
    disp('no log file for this simulation')
end

SIGN_PSI0=1;
if psi_star_final_profile(1)<0
    %psi_star_2D_evol_interp=-psi_star_2D_evol_interp;
    psi_star_2D_evol_lin=-psi_star_2D_evol_lin;
	psi_star_initial=-psi_star_initial;
	psi_star_final=-psi_star_final;
	psi_star_final_profile=-psi_star_final_profile;
    SIGN_PSI0=-1
end

if ISKINK==1
	NB_FRAMES=100
else
	NB_FRAMES=101
end

min_Psih=min(psi_star_initial);