%SAVENAME=strcat('initial',SAVENAME)
if SAVE_DATA_FILE==1
    save(SAVENAME,'alphas_ejected','alphas_momentum','time_scale','vpll_output','Xpos_gc_output','Zpos_gc_output','phipos_output','alphas_pos_x','alphas_pos_z','alphas_pos_phi',...
	'v_X','v_Z','v_phi','alphas_vpll','alphas_mm','alphas_psi','alphas_Ekin','alphas_pphi0','Nalphas_simulated');
end
save(SAVENAME,'-append','PROCESS_NUMBER');

clearvars -EXCEPT PROCESS_NUMBER SAVENAME DISTNAME SAVENAME_PRE SAVENAME_STATS
clc
format compact
reset_data_analysis_environment;
clear q_XZ_map

%DISTNAME=strcat('initial_NBI_1MEV_D_distribution',num2str(PROCESS_NUMBER),'.mat')
%DISTNAME=strcat('initial_alphas_vA_distribution',num2str(PROCESS_NUMBER),'.mat')
load(DISTNAME,'mHe','ZHe');
INPUTNAME=SAVENAME


wrap_precession_file;

%SAVENAME=strcat('initial_NBI_1MEV_D_pre_collapse',num2str(PROCESS_NUMBER),'.mat')
SAVENAME=SAVENAME_PRE

save_distribution_without_precession_data;
save(SAVENAME,'-append','PROCESS_NUMBER');
clear alphas_pos_x alphas_pos_z alphas_pos_phi alphas_vpll alphas_psi v_X v_Z v_phi

INPUTNAME
%SAVENAME=strcat('initial_NBI_1MEV_D_precession_stats',num2str(PROCESS_NUMBER),'.mat')
SAVENAME=SAVENAME_STATS


