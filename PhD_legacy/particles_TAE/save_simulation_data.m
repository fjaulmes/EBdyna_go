        toc
        
        tic
        clear outcast recast
        % correct the particles that are out of simulation domain
        % by giving them a position randomly in the initial distribution
        outcast=find(alphas_psi>NB_PSI-2);
        recast=randi(Nalphas_simulated,size(outcast,1),1);
        if (~isempty(outcast))
            reposition_lost_particles_3DG;
            alphas_ejected(outcast)=1;
        end

        disp('------------------------------------------------')
        disp(strcat('time_step # ',num2str(time_step),' ; time = ',num2str(time)));
        disp(strcat('time_stamp # ',num2str(time_stamp),' ; time = ',num2str(time)));
        disp(strcat('number of repositionned particles =  ',num2str(size(outcast,1))));
        disp(strcat('total number of ejected particles =  ',num2str(size(find(alphas_ejected),1))));

        save (SAVENAME ,'time_stamp','Epot_output','Etot_output','vparallel_output',...
             'psi_value_output','phipos_output','theta_output','Ekin_output','Ekin_gc_output',...
             'Xpos_part_output','vphi_output','pphi_output','time_scale','alphas_ejected',...
             'Epart_tot_rel_evol','Epart_tot_vD_rel_evol','power_exchange_evol','WTAE_evol',...
             'Ekin_part_evol','Etot_part_evol','Etot_th_part_evol',...
             'gamma_TAE_evol','gamma_TAE_vD_evol');
        disp('------------------------------------------------')
		
		
		%if (PROCESS_NUMBER>=20)&&(PROCESS_NUMBER<=23)
		%     save(strcat('Ekin_bug_',num2str(PROCESS_NUMBER),'.mat'),Ekin_output,alphas_Ekin);
	    %end
		