

initialize_eq_sim_parameters_RMP;
initialize_eq_sim_maps;

initialize_eq_simulation_arrays;


% benchmarking time between two data recordings
disp('**************************************************************');
tic

time=0;


for time_step=1:NB_TIME_STEPS
   
    if ACTIVATE_CENTRIFUG_TERM==0
		time_step_integration_GT_eq;
		time_step_integration_GT_eq;

		time_step_integration_GT_eq;
		time_step_integration_GT_eq;
		
		time_step_integration_GT_eq;
		time_step_integration_GT_eq;
		
		time_step_integration_GT_eq;
		time_step_integration_GT_eq;

		time_step_integration_GT_eq;
		time_step_integration_GT_eq;
	else
		time_step_integration_GT_eq_centrifug;
		time_step_integration_GT_eq_centrifug;

		time_step_integration_GT_eq_centrifug;
		time_step_integration_GT_eq_centrifug;
		
		time_step_integration_GT_eq_centrifug;
		time_step_integration_GT_eq_centrifug;
		
		time_step_integration_GT_eq_centrifug;
		time_step_integration_GT_eq_centrifug;

		time_step_integration_GT_eq_centrifug;
		time_step_integration_GT_eq_centrifug;
	end
    adapt_speed_Ekin_G;
    clear outcast recast
    % correct the particles that are out of simulation domain
    % by giving them a position randomly in the initial distribution
    outcast=find(alphas_psi>SIMULATION_RADIAL_LIMIT);
    recast=randi(Nalphas_simulated,size(outcast,1),1);
    if (~isempty(outcast))
        reposition_lost_particles_3DG;
        alphas_ejected(outcast)=1;
    end

    if (mod(time_step,TIME_STAMP_PRECISION)==0)
        record_time_stamp;
	end


    if (SAVE_DATA_FILE==1) && (mod(time_step,RECORD_PRECISION*10)==0)
        toc
        tic

        disp('------------------------------------------------')
        disp(strcat('time_step # ',num2str(time_step),' ; time = ',num2str(time)));
        disp(strcat('time_stamp # ',num2str(time_stamp),' ; time = ',num2str(time)));
        disp(strcat('number of repositionned particles =  ',num2str(size(outcast,1))));
        disp(strcat('total number of ejected particles =  ',num2str(size(find(alphas_ejected),1))));

    end

end

toc
disp('**************************************************************');
disp('done');
Nalphas_simulated

save_data_file_RMP;

disp('**************************************************************');
disp('done');

exit
