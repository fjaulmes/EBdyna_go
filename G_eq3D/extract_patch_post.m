for n=1:16


	PROCESS_NUMBER=n

    reset_data_analysis_environment_post;
	clear q_XZ_map
	%PROCESS_VAL=(round(PROCESS_NUMBER)-0.5)*10


	%DISTNAME=strcat('initial_flatD_',num2str(PROCESS_VAL),'keV_distribution.mat')
	DISTNAME=strcat('initial_alphas_MB_D_distribution',num2str(PROCESS_NUMBER),'.mat')
	%DISTNAME=strcat('initial_NBI60keV_transp_D_distribution',num2str(PROCESS_NUMBER),'.mat')
	SAVENAME=strcat('finalG_MB_D_precession',num2str(PROCESS_NUMBER),'.mat')


	%DISTNAME=strcat('initial_NBI60keV_transp_D_distribution',num2str(PROCESS_NUMBER),'.mat')
	%SAVENAME=strcat('initial_NBI60keV_precession',num2str(PROCESS_NUMBER),'.mat')
	load(DISTNAME);
	INPUTNAME=SAVENAME


	wrap_precession_file;

	%SAVENAME=strcat('initial_alphas_TAE_pre_collapse',num2str(PROCESS_NUMBER),'.mat')
    %SAVENAME=strcat('initial_flatD_',num2str(PROCESS_VAL),'keV_pre_collapse.mat')
    %SAVENAME=strcat('final_NBI60keV_post_collapse',num2str(PROCESS_NUMBER),'.mat')
	SAVENAME=strcat('final_MB_D_post_collapse',num2str(PROCESS_NUMBER),'.mat')

	save_distribution_without_precession_data;
	clear alphas_pos_x alphas_pos_z alphas_pos_phi alphas_vpll alphas_psi v_X v_Z v_phi

	INPUTNAME
	%SAVENAME=strcat('initial_alphas_TAE_precession_stats',num2str(PROCESS_NUMBER),'.mat')
    SAVENAME=strcat('final_MB_D_precession_stats',num2str(PROCESS_NUMBER),'.mat')
    %SAVENAME=strcat('final_NBI60keV_precession_stats',num2str(PROCESS_NUMBER),'.mat')
	%SAVENAME=strcat('initial_NBI60kev_precession_stats',num2str(PROCESS_NUMBER),'.mat')

	extract_precession_information_post;

end