%% Initialize
initialize_eq_sim_parameters; % Parameters (e.g. time steps)
initialize_eq_sim_maps;       % Maps (e.g.  RMP-field)

initialize_eq_simulation_arrays;

% benchmarking time between two data recordings
disp('**************************************************************');

%% MAIN loop
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
    
    %% Correct values
    %calculating integer time step velocity values
    if CALCULATE_TRUE_PPHI==0
        update_GT_3D_eq;
        v_X_step=(0.5*v_X+0.5*v_X_prev);
        v_Z_step=(0.5*v_Z+0.5*v_Z_prev);
        v_phi_step=(0.5*v_phi+0.5*v_phi_prev);
        v_X=v_X_prev;
        v_Z=v_Z_prev;
        v_phi=v_phi_prev;
    else
        alphas_pphi0=alphas_pphi0_prev+alphas_Delta_pphi;
        alphas_pphi0_prev=alphas_pphi0;
        alphas_Delta_pphi=zeros(Nalphas_simulated,1);
    end
    %applying vphi correction
    if PPHI_CORR_FACTOR~=0
        v_phi_step_recalc=(alphas_pphi0+(ZHe)*alphas_psi_value)./(alphas_Rpos);
        v_phi_step_recalc=(eV/mHe)*v_phi_step_recalc;
        adapt_speed_pphi_G;
    end
    %applying Ekin correction
    if ECOEF~=0
        adapt_speed_Ekin_G;
    end
    %applying vphi correction
    if PPHI_CORR_FACTOR~=0
        update_GT_3D_eq;
        v_X_step=(0.5*v_X+0.5*v_X_prev);
        v_Z_step=(0.5*v_Z+0.5*v_Z_prev);
        v_phi_step=(0.5*v_phi+0.5*v_phi_prev);
        v_X=v_X_prev;
        v_Z=v_Z_prev;
        v_phi=v_phi_prev;
        adapt_speed_pphi_G;
    end
    
    %% Record
    if (mod(time_step,TIME_STAMP_PRECISION)==0)
        record_time_stamp;
    end
    
    %% SAVE DATA FILE intermediately
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
        
        %         filename=strcat('alphas_pre_Grecord_',num2str(time_step),'.mat');
        %        save(strcat(SAVENAME,num2str(time_step)),'time','alphas_pos_x','alphas_pos_z','alphas_pos_phi','alphas_vpll','alphas_mm_part','alphas_psi','alphas_ejected');
        %        disp('------------------------------------------------')
    end
end

toc
disp('**************************************************************');
disp('done');
Nalphas_simulated

save_data_file;

extract_precession_information;
disp('**************************************************************');
disp('done');

exit % Notation: exit prevents exitcode to work