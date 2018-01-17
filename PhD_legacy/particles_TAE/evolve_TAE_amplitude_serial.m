
adapt_speed_Ekin_G;
% v_phi_next=v_phi_next+0.1*v_phi_corr_integ+0.1*(v_phi_step_recalc-v_phi_step);
v_phi_next=v_phi_next+PPHI_CORR_FACTOR*(v_phi_step_recalc-v_phi_step);
% % 
adapt_speed_Ekin_G;
v_phi_next=v_phi_next+PPHI_CORR_FACTOR*(v_phi_step_recalc-v_phi_step);
% 
adapt_speed_Ekin_G;
v_phi_next=v_phi_next+PPHI_CORR_FACTOR*(v_phi_step_recalc-v_phi_step);

% 			deltaE_Epart=sum(alphas_weight.*delta_E)
%             deltaE_Epart_th=sum(alphas_weight.*delta_E_th)

    % TAE growth rate from power transfer and mode energy
    if (mod(time_step-1,GYRO_ORBIT_PRECISION)==0)
	    disp('-------------------------------------------------');

	    EJECTED_POP=find(alphas_ejected);

	    time_stamp_go=ceil((time_step)/GYRO_ORBIT_PRECISION);

        if CALCULATE_DELTAE_POWER==1
            % keep energy change from GC calculations [impossible to use GO position]
            delta_Ekin=delta_Ekin_vD;

			% transfer only total energy change to the wave
            % divide by TIME_GO_SIZE for better Gyro orbit stats
            part2wave_power=-delta_E/(TIME_GO_SIZE*h);
 			part2wave_power(EJECTED_POP)=0;

            if CALCULATE_VD_POWER==0
                part2wave_vD_power=part2wave_power;
            end
        end
        
        if CALCULATE_VD_POWER==1
            part2wave_vD_power=-delta_Ekin_vD/(TIME_GO_SIZE*h);			
            part2wave_vD_power(EJECTED_POP)=0;
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
            P_VD_PART_TAE=(NB_PART_RESCALE*eV)*sum(alphas_weight.*part2wave_vD_power);
            P_PART_TAE=-((NB_PART_RESCALE*eV)/(TIME_GO_SIZE*h))*sum(alphas_weight.*delta_E);
        end
        
 
        
        % TAE growth rate
        gamma_TAE_prev=gamma_TAE_local;
        if (TAE_amplitude>=MINIMUM_TAE_AMPLITUDE) && (f_counter>2)
            % only the vD.E part is transfered to the wave ?
			gamma_TAE_vD_local=0.5*(P_VD_PART_TAE/W_TAE);
			gamma_TAE_local=0.5*(P_PART_TAE/W_TAE);

            gamma_TAE_local=max(gamma_TAE_local,gamma_TAE_min);
            gamma_TAE_local=min(gamma_TAE_local,gamma_TAE_max);
        else
            gamma_TAE_local=0;
        end
        gamma_TAE_local=1.5*gamma_TAE_local-0.5*gamma_TAE_prev;
        
time_step_values=zeros(NB_PROCESS,1);
time_step_values(PROCESS_NUMBER)=time_step;

% NB_PROCESS parallel processes
% waiting for at least previous synchro
%before wirting anything

gamma_TAE_values(PROCESS_NUMBER)=gamma_TAE_local;
gamma_TAE_vD_values(PROCESS_NUMBER)=gamma_TAE_vD_local;
synchro_flag=0;


        gamma_TAE=sum(gamma_TAE_values);
        gamma_TAE_vD=sum(gamma_TAE_vD_values);
        %     gamme_TAE_half=max(0.5*(gamme_TAE+gamme_TAE_prev),0);
        %     gamma_TAE_half=0.5*(gamma_TAE+gamma_TAE_prev);
        if (gamma_TAE_prev~=0)
            if EVOLVE_TAE_AMPLITUDE==1
                % multiply by TIME_GO_SIZE for better Gyro orbit stats
                if USE_VD==0
                    if (TAE_amplitude>=MINIMUM_TAE_AMPLITUDE)
                        TAE_amplitude=max(TAE_amplitude+TAE_amplitude*gamma_TAE*TIME_GO_SIZE*h,0);
                    else
                        TAE_amplitude=max(TAE_amplitude+TAE_amplitude*gamma_TAE*TIME_GO_SIZE*h,MINIMUM_TAE_AMPLITUDE);
                    end
                else
                    if (TAE_amplitude>=MINIMUM_TAE_AMPLITUDE)
                        TAE_amplitude=max(TAE_amplitude+TAE_amplitude*gamma_TAE_vD*TIME_GO_SIZE*h,0);
                    else
                        TAE_amplitude=max(TAE_amplitude+TAE_amplitude*gamma_TAE_vD*TIME_GO_SIZE*h,MINIMUM_TAE_AMPLITUDE);
                    end
                end
            end
            % the vD evolution makes little sense here....
            % it only takes the local process particles into account
            Epart_tot_vD_rel_ini=Epart_tot_vD_rel_ini-P_VD_PART_TAE*(TIME_GO_SIZE*h/(NB_PART_RESCALE*eV));
            Epart_tot_rel_ini=Epart_tot_rel_ini+(NB_PART_RESCALE*eV)*sum(alphas_weight.*delta_E);
% 			deltaE_Epart1=sum(alphas_weight.*dE)
% 			deltaE_Epart=sum(alphas_weight.*delta_E)
%             deltaE_Epart_th=sum(alphas_weight.*delta_E_th)
        end
        W_TAE_raw=(TAE_amplitude)^2;
        W_TAE=W_TAE_raw*WTAE_AVG;
        quasi_lin_flag=1;

        part2wave_power_evol(time_stamp_go)=P_PART_TAE;
        WTAE_evol(time_stamp_go)=W_TAE/WTAE_AVG;
        Ekin_part_evol(time_stamp_go)=(NB_PART_RESCALE*eV)*sum(alphas_weight.*alphas_Ekin);
        Etot_part_evol(time_stamp_go)=(NB_PART_RESCALE*eV)*sum(alphas_weight.*alphas_Etot);

        
        Etot_th_part_evol(time_stamp_go)=(NB_PART_RESCALE*eV)*sum(alphas_weight.*alphas_Etot_th);
        gamma_TAE_evol(time_stamp_go)=gamma_TAE;
        gamma_TAE_vD_evol(time_stamp_go)=gamma_TAE_vD;
        Epart_tot_vD_rel_evol(time_stamp_go)=Epart_tot_vD_rel_ini;
        Epart_tot_rel_evol(time_stamp_go)=Epart_tot_rel_ini;
        
        alphas_pphi0=alphas_pphi0_prev+delta_pphi_th;

        
        alphas_Etot_prev=alphas_Etot;
        alphas_Etot_th_prev=alphas_Etot_th;
        alphas_pphi0_prev=alphas_pphi0;
        delta_E=zeros(Nalphas_simulated,1);  
        delta_E_th=zeros(Nalphas_simulated,1);  
        delta_pphi_th=zeros(Nalphas_simulated,1);  
    
    end
