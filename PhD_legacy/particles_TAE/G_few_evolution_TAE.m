
close all;
pause(0.2);

REINIT_ALL_MAPS=0;
REINIT_SIMULATION_DATA=1;

if REINIT_ALL_MAPS==1
    clear all
    initialize_folder_names_test;
    filename=strcat(DATA_FOLDER,'physics_constants.mat');
    load(filename);
    filename=strcat(DATA_FOLDER,'tokamak_map_dimensions.mat');
    load(filename);
    filename=strcat(DATA_FOLDER,'TAE_data.mat');
    load(filename);

    DPHI=TAE_angle/(NB_PHI-1);
    size_r=pTAE_sup-pTAE_inf+1;

    mid_X_large=mid_X;
    mid_Z_large=mid_Z;    
    calculate_TAE_drift_speed_maps;
    clear all;
    REINIT_SIMULATION_DATA=1;
end

if REINIT_SIMULATION_DATA==1
    clear all;
    clc
    format compact
    initialize_folder_names_test;
    
    filename=strcat(DATA_FOLDER,'physics_constants.mat');
    load(filename);
    filename=strcat(DATA_FOLDER,'XZsmall_fields_tokamak_TAE.mat');
    load(filename);
    filename=strcat(DATA_FOLDER,'motions_map_dimensions.mat');
    load(filename);
%     filename=strcat(DATA_FOLDER,'psi_profiles_kadomtsev.mat');
%     load(filename,'psi_star_initial');

    filename=strcat(DATA_FOLDER,'tokamak_PR_map.mat');
%     load(filename);
    load(filename, 'Rpos_PR_map');
    load(filename, 'radial_r_value_flux');
    load(filename,'BX_PR_map')
    load(filename, 'BZ_PR_map')
    filename=strcat(DATA_FOLDER,'B_fields.mat');
    load(filename,'Bpol_PR_map');
    load(filename,'Btor_PR_map');

    filename=strcat(DATA_FOLDER,'TAE_data.mat');
    load(filename);
 
    REINIT_SIMULATION_DATA=1;
   
end

PROCESS_NUMBER=1

initialize_TAEsim_parameters_serial;
initialize_TAE_sim_maps;
% coupling_m_mp1=1.0

if (USE_VD==1)&&(CALCULATE_VD_POWER==0)
    disp('inconsitent simulation parameters for vD !')
    CALCULATE_VD_POWER=1
end

initialize_all_simulation_arrays_serial;

for time_step=1:NB_TIME_STEPS

    %%
    if (mod(time_step,GYRO_ORBIT_PRECISION)==0)
        delta_E=zeros(Nalphas_simulated,1);
        delta_Epot=zeros(Nalphas_simulated,1);
        delta_Ekin_vD=zeros(Nalphas_simulated,1);
        P_VD_PART_TAE=0;
        P_PART_TAE=0;
    end
    
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
    evolve_TAE_amplitude;

    
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

