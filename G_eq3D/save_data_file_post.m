% SAVENAME=strcat('final',SAVENAME)
if SAVE_DATA_FILE==1
    save(SAVENAME,'alphas_ejected','time_scale','vpll_output','Xpos_gc_output','Zpos_gc_output','phipos_output','alphas_pos_x','alphas_pos_z','alphas_pos_phi',...
	'v_X','v_Z','v_phi','alphas_vpll','alphas_mm','alphas_psi','alphas_Ekin','alphas_pphi0','Nalphas_simulated');
end

clearvars -EXCEPT PROCESS_NUMBER DISTNAME SAVENAME
clc
format compact
reset_data_analysis_environment_post;

% DISTNAME=strcat('initial_NBI_1MEV_D_distribution',num2str(PROCESS_NUMBER),'.mat')
%DISTNAME=strcat('initial_alphas_MB_D_distribution',num2str(PROCESS_NUMBER),'.mat')
load(DISTNAME);
INPUTNAME=SAVENAME


wrap_precession_file;

%SAVENAME=strcat('final_MB_D_post_collapse',num2str(PROCESS_NUMBER),'.mat')
SAVENAME=strcat('final_NBI60keV_post_collapse',num2str(PROCESS_NUMBER),'.mat')

save_distribution_without_precession_data;
clear alphas_pos_x alphas_pos_z alphas_pos_phi alphas_vpll alphas_psi v_X v_Z v_phi

INPUTNAME
SAVENAME=strcat('final_NBI60keV_precession_stats',num2str(PROCESS_NUMBER),'.mat')


