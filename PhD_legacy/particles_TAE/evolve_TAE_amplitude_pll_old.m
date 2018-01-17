


%     delta_Epot=alphas_Epot_gc-alphas_Epot_gc_part_prev;
%    alphas_Epot(OUTER_PART)=0;
%    alphas_Epot(INNER_PART)=lininterp3(E_potential_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
%     delta_Epot=alphas_Epot-alphas_Epot_prev;
%     delta_Epot=alphas_Epot_gc-alphas_Epot_gc_part_prev;

    % TAE growth rate from power transfer and mode energy
    if (mod(time_step-1,GYRO_ORBIT_PRECISION)==0)
        if CALCULATE_DELTAE_POWER==1
            % keep energy change from GC calculations [impossible to use GO position]
            delta_Ekin=delta_Ekin_vD;

			% transfer only total energy change to the wave
            % divide by TIME_GO_SIZE for better Gyro orbit stats
            part2wave_power=-delta_E/TIME_GO_SIZE/h;
            %         alphas_Ekin(INNER_PART_GC)=alphas_Ekin_th(INNER_PART_GC);
            if CALCULATE_VD_POWER==0
                part2wave_vD_power=part2wave_power;
            end
        end
        
        if CALCULATE_VD_POWER==1
            part2wave_vD_power=-delta_Ekin_vD/TIME_GO_SIZE/h;
        end
        

        
        if (f_counter<3)
            % no drive initially
            P_VD_PART_TAE=0;
            P_PART_TAE=0;
        else
            if USE_VD==0
                Pexchange_part_record=Pexchange_part_record+part2wave_power;
            else
                Pexchange_part_record=Pexchange_part_record+part2wave_vD_power;
            end
            P_VD_PART_TAE=(NB_PART_RESCALE*eV)*sum(part2wave_vD_power);
            P_PART_TAE=(NB_PART_RESCALE*eV)*sum(part2wave_power);
        end
        
        % TAE growth rate
        gamma_TAE_prev=gamma_TAE_local;
        if (TAE_amplitude>=MINIMUM_TAE_AMPLITUDE) && (f_counter>2)
            % only the vD.E part is transfered to the wave ?
            if USE_VD==0
                gamma_TAE_local=0.5*(P_PART_TAE/W_TAE)/WTAE_AVG;
            else
                gamma_TAE_local=0.5*(P_VD_PART_TAE/W_TAE)/WTAE_AVG;
            end
            gamma_TAE_local=max(gamma_TAE_local,gamma_TAE_min);
            gamma_TAE_local=min(gamma_TAE_local,gamma_TAE_max);
        else
            gamma_TAE_local=0;
        end
        
binary_filename=strcat('pll_gamma',num2str(PROCESS_NUMBER),'.bin');
fileID = fopen(binary_filename,'wb');
% could not separate data types....
% so store time_step as double
fwrite(fileID,time_step,'double');
fwrite(fileID,gamma_TAE_local,'double');
fclose(fileID);
time_step_values(PROCESS_NUMBER)=time_step;
gamma_TAE_values(PROCESS_NUMBER)=gamma_TAE_local;
synchro_flag=0;

% 8 parallel processes
while synchro_flag==0
    pause(0.3)
    for p=1:8
        if (p~=PROCESS_NUMBER) && (time_step_values(p)~=time_step)
            binary_filename=strcat('pll_gamma',num2str(p),'.bin');
            fileID = fopen(binary_filename,'r');
			if fileID~=-1
                FD=fread(fileID,2,'double');
                if ~isempty(FD)
                   time_step_values(p)=FD(1);
                   gamma_TAE_values(p)=FD(2);
                end
			    fclose(fileID);
			end            
        end
    end
    synchro_flag=1;
    for p=1:8       
        if (time_step_values(p)~=time_step)
            synchro_flag=0;
        end
    end
    % don't want to keep accessing file system!
    pause(0.2)
end

        gamma_TAE=sum(gamma_TAE_values);
        %     gamme_TAE_half=max(0.5*(gamme_TAE+gamme_TAE_prev),0);
        %     gamma_TAE_half=0.5*(gamma_TAE+gamma_TAE_prev);
        if (gamma_TAE_prev~=0)
            if EVOLVE_TAE_AMPLITUDE==1
                % multiply by TIME_GO_SIZE for better Gyro orbit stats
                if (TAE_amplitude>=MINIMUM_TAE_AMPLITUDE)
                    TAE_amplitude=max(TAE_amplitude+TAE_amplitude*gamma_TAE*TIME_GO_SIZE*h,0);
                else
                    TAE_amplitude=max(TAE_amplitude+TAE_amplitude*gamma_TAE*TIME_GO_SIZE*h,MINIMUM_TAE_AMPLITUDE);
                end
            end
            % the vD evolution makes little sense here....
            % it only takes the local process particles into account
%             Epart_tot_vD_rel_ini=Epart_tot_vD_rel_ini-P_VD_PART_TAE*TIME_GO_SIZE*h;
            Epart_tot_vD_rel_ini=Epart_tot_vD_rel_ini-2*W_TAE*gamma_TAE*TIME_GO_SIZE*h;
            Epart_tot_rel_ini=Epart_tot_rel_ini-2*W_TAE*gamma_TAE*TIME_GO_SIZE*h;
            % Ekin evolution follows gc equation
            if USE_VD==0
                % alphas_Ekin=max(alphas_Ekin-h*TIME_GO_SIZE*part2wave_power,0);
				alphas_Ekin=alphas_Ekin+delta_Ekin;
            else
                alphas_Ekin=max(alphas_Ekin-h*TIME_GO_SIZE*part2wave_vD_power,0);
            end
        end
        W_TAE=(TAE_amplitude)^2;
%         quasi_lin_flag=1;
    end
%     alphas_Epot_prev=alphas_Epot;
%     alphas_Epot_gc_part_prev=alphas_Epot_gc;
