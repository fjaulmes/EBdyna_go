	if (time_step<NB_TIME_STEPS)
       SAVENAME_TS=strcat(SAVENAME,'_',num2str(time_step),'_',num2str(PROCESS_NUMBER),'.mat')
	   save (SAVENAME_TS,'mHe','ZHe','time','frame_rank_precise','pos_X_gc','pos_Z_gc','alphas_pos_x','alphas_pos_z',...
        'alphas_pos_phi','alphas_omega','alphas_vpll','alphas_mm',...
        'alphas_psi','alphas_pphi0','alphas_psi_value_corr','alphas_Etot','alphas_Ekin','alphas_ejected');

	else
        SAVENAME_TS=strcat(SAVENAME,'_',num2str(PROCESS_NUMBER),'.mat')
	    save (SAVENAME_TS,'mHe','ZHe','time','frame_rank_precise','pos_X_gc','pos_Z_gc','alphas_pos_x','alphas_pos_z',...
        'alphas_eject_posX','alphas_eject_posZ','alphas_eject_vpll','alphas_pos_phi','alphas_omega','alphas_vpll','alphas_mm',...
        'alphas_psi','alphas_pphi0','alphas_psi_value_corr','alphas_Etot','alphas_Ekin','alphas_ejected');
	end
%    else
        %SAVENAME_TS=strcat(SAVENAME,'_',num2str(PROCESS_NUMBER),'.mat')
%        SAVENAME_TS=strcat(SAVENAME,'_',num2str(time_step),'_',num2str(PROCESS_NUMBER),'.mat')

%	    save (SAVENAME_TS,'mHe','ZHe','time','frame_rank_precise','pos_X_gc','pos_Z_gc','alphas_pos_x','alphas_pos_z',...
%        'alphas_eject_posX','alphas_eject_posZ','alphas_eject_vpll','alphas_pos_phi','alphas_omega','v_X','v_Z','v_phi','alphas_vpll','alphas_mm',...
%        'alphas_psi','alphas_pphi0','alphas_psi_value_corr','alphas_Etot','alphas_Ekin','alphas_ejected');
%    end