
close all;
pause(0.2);

REINIT_ALL_MAPS=0;
REINIT_SIMULATION_DATA=1;

if REINIT_ALL_MAPS==1
    clear all
    initialize_folder_names;
    filename=strcat(DATA_FOLDER,'physics_constants.mat');
    load(filename);
    filename=strcat(DATA_FOLDER,'tokamak_map_dimensions.mat');
    load(filename);

    mid_X_large=mid_X;
    mid_Z_large=mid_Z;    
    NB_PHI=129;
    DPHI=2*pi/(NB_PHI-1);
    NB_PHI_DATA_HALF=round(0.5*(NB_PHI-1));
    calculate_post_collapse_drift_speed_maps;
    clear all;
    REINIT_SIMULATION_DATA=1;
end

if REINIT_SIMULATION_DATA==1
    clear all;
    clc
    format compact
    initialize_folder_names;
    
    filename=strcat(DATA_FOLDER,'physics_constants.mat');
    load(filename);
    filename=strcat(DATA_FOLDER,'XZsmall_fields_tokamak_post_collapse.mat');
    load(filename);
    filename=strcat(DATA_FOLDER,'motions_map_dimensions.mat');
    load(filename);
%     filename=strcat(DATA_FOLDER,'psi_star_evol.mat');
%     load(filename, 'psi_star_2D_evol_interp');
%     psi_star_map_phi_evol=interp1(1:1001,psi_star_2D_evol_interp,1:101,'cubic');

    filename=strcat(DATA_FOLDER,'flux_geometry.mat');
    load(filename, 'dr_X_PR_map');
    load(filename, 'dr_Z_PR_map');

    REINIT_SIMULATION_DATA=1;
end

TOROIDAL_FIELD_DIRECTION=sign(mean(mean(Bphi_XZsmall_map)))

PROCESS_NUMBER=6


% initialize completely the maps for pre-collapse calculations

initialize_eq_sim_parameters_post;
initialize_eq_sim_maps;

initialize_eq_simulation_arrays_post;




    
% benchmarking time between two data recordings
disp('**************************************************************');
tic

time=0;


for time_step=1:NB_TIME_STEPS
   
    time_step_integration_GT_eq;
    time_step_integration_GT_eq;

    time_step_integration_GT_eq;
    time_step_integration_GT_eq;
    
    time_step_integration_GT_eq;
    time_step_integration_GT_eq;
    
    time_step_integration_GT_eq;
    time_step_integration_GT_eq;

    time_step_integration_GT_eq;
%     time_step_integration_GT_collapse;
    time_step_integration_GT_eq;

    
%    estimate_half_time_step_positions_eq_G;

%    adapt_speed_pphi_G;
%    v_phi=v_phi_tilde;


    if (mod(time_step,TIME_STAMP_PRECISION)==0)
        record_time_stamp;
    end

    if (SAVE_DATA_FILE==1) && (mod(time_step,RECORD_PRECISION*10)==0)
        toc
        tic
        clear outcast recast
        % correct the particles that are out of simulation domain
        % by giving them a position randomly in the initial distribution
        outcast=find(alphas_psi>SIMULATION_RADIAL_LIMIT);
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
	end



end

toc
disp('**************************************************************');
disp('done');
Nalphas_simulated

save_data_file_post;

extract_precession_information_post;
disp('**************************************************************');
disp('done');
