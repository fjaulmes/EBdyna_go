
close all;
pause(0.2);




initialize_TAEsim_parameters;
initialize_TAE_sim_maps;

if (USE_VD==1)&&(CALCULATE_VD_POWER==0)
    disp('inconsitent simulation parameters for vD !')
    CALCULATE_VD_POWER=1
end

initialize_all_simulation_arrays;

for time_step=1:NB_TIME_STEPS
    
    time_step_integration_G_TAE;
    time_step_integration_G_TAE;
    time_step_integration_G_TAE;
    time_step_integration_G_TAE;
    time_step_integration_G_TAE;
    time_step_integration_G_TAE;
    time_step_integration_G_TAE;
    time_step_integration_G_TAE;
    time_step_integration_G_TAE;
    time_step_integration_G_TAE;

	time_step_integration_G_TAE;
    time_step_integration_G_TAE;
    time_step_integration_G_TAE;
    time_step_integration_G_TAE;
    time_step_integration_G_TAE;
    time_step_integration_G_TAE;
    time_step_integration_G_TAE;
    time_step_integration_G_TAE;
    time_step_integration_G_TAE;
    time_step_integration_G_TAE;

    % TAE growth rate from power transfer and mode energy
    evolve_TAE_amplitude_pll;

    if (mod(time_step-1,TIME_STAMP_PRECISION)==0)
        record_time_stamp;
    end

    if (SAVE_DATA_FILE==1) && (mod(time_step,1000)==0)
        save_simulation_data;
    end



end

save_simulation_data;

disp('**************************************************************');
disp('done');
Nalphas_simulated

