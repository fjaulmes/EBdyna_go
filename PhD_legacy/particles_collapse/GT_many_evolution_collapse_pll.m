


initialize_collapse_parameters;

initialize_sim_maps;
initialize_simulation_arrays;



% benchmarking time between two data recordings
disp('**************************************************************');
tic

time=0;


for time_step=1:NB_TIME_STEPS
    
    time_step_integration_GT_collapse;
    time_step_integration_GT_collapse;
    time_step_integration_GT_collapse;
    time_step_integration_GT_collapse;
    time_step_integration_GT_collapse;
    
    time_step_integration_GT_collapse;
    time_step_integration_GT_collapse;
    time_step_integration_GT_collapse;
    time_step_integration_GT_collapse;
    time_step_integration_GT_collapse;
    
    time_step_integration_GT_collapse;
    time_step_integration_GT_collapse;
    time_step_integration_GT_collapse;
    time_step_integration_GT_collapse;
    time_step_integration_GT_collapse;
    
    time_step_integration_GT_collapse;
    time_step_integration_GT_collapse;
    time_step_integration_GT_collapse;
    time_step_integration_GT_collapse;
    time_step_integration_GT_collapse;

	outcast=find(alphas_psi>Nradial-2);
	if (~isempty(outcast))
		alphas_eject_posX(outcast)=alphas_pos_x(outcast);
		alphas_eject_posZ(outcast)=alphas_pos_z(outcast);
		alphas_eject_vpll(outcast)=alphas_vpll(outcast);
		alphas_pos_x(outcast)=0;
		alphas_pos_z(outcast)=0;
        alphas_ejected(outcast)=1;
        disp(strcat('number of ejected particles (reset to 0,0) =  ',num2str(size(outcast,1))));
	end
    
    
    record_time_stamp;
    
    if (SAVE_DATA_FILE==1) && (mod(time_step,RECORD_TS)==0)
        toc
        tic
        clear outcast recast
        % correct the particles that are out of simulation domain
        % by giving them a position randomly in the initial distribution
        outcast=find(alphas_psi>simulation_size_r+22);
        recast=randi(Nalphas_simulated,size(outcast,1),1);
        if (~isempty(outcast))
            reposition_lost_particles_3DG;
            alphas_ejected(outcast)=1;
        end
        
        disp('------------------------------------------------')
        disp(strcat('time_step # ',num2str(time_step),' ; time = ',num2str(time)));
        disp(strcat('number of repositionned particles =  ',num2str(size(outcast,1))));
        disp(strcat('total number of ejected particles =  ',num2str(size(find(alphas_ejected),1))));
        
        pos_X_gc=(mHe/eV)*(1/ZHe)*(v_Z_step.*bphi-v_phi_step.*bZ)./alphas_Bfield;
        pos_Z_gc=(mHe/eV)*(1/ZHe)*(v_phi_step.*bX-v_X_step.*bphi)./alphas_Bfield;
		pos_X_gc=alphas_pos_x+pos_X_gc;
		pos_Z_gc=alphas_pos_z+pos_Z_gc;

        % filename=strcat('flatD_90keV_collapse_Glisa_fc2h2_G260214record_',num2str(time_step),'.mat')
        save_data_file;
        disp('------------------------------------------------')
    end
    
    
    
end
        pos_X_gc=(mHe/eV)*(1/ZHe)*(v_Z_step.*bphi-v_phi_step.*bZ)./alphas_Bfield;
        pos_Z_gc=(mHe/eV)*(1/ZHe)*(v_phi_step.*bX-v_X_step.*bphi)./alphas_Bfield;
		pos_X_gc=alphas_pos_x+pos_X_gc;
		pos_Z_gc=alphas_pos_z+pos_Z_gc;

toc


if SAVE_DATA_FILE==1
    save_data_file;
end
disp('**************************************************************');
disp('done');
Nalphas_simulated
