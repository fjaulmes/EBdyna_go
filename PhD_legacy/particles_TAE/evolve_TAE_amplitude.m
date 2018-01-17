
adapt_speed_Ekin_G;
v_phi_next=v_phi_next+v_phi_corr_integ+0.01*(v_phi_step_recalc-v_phi_step);

adapt_speed_Ekin_G;
v_phi_next=v_phi_next+0.01*(v_phi_step_recalc-v_phi_step);


%     delta_Epot=alphas_Epot_gc-alphas_Epot_gc_part_prev;
%    alphas_Epot(OUTER_PART)=0;
%    alphas_Epot(INNER_PART)=lininterp3(E_potential_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
%     delta_Epot=alphas_Epot-alphas_Epot_prev;
%     delta_Epot=alphas_Epot_gc-alphas_Epot_gc_part_prev;

    % TAE growth rate from power transfer and mode energy
    if (mod(time_step-1,GYRO_ORBIT_PRECISION)==0)
	    EJECTED_POP=find(alphas_ejected);
        
	    time_stamp_go=ceil((time_step)/GYRO_ORBIT_PRECISION);
		if f_counter>2
           alphas_pphi0_recalc_prev=alphas_pphi0_recalc;
 	       alphas_pphi0_recalc=pphi_record/TIME_GO_SIZE;
		   alphas_pphi0_recalc=0.9*alphas_pphi0+0.1*alphas_pphi0_recalc;
%            alphas_pphi0=alphas_pphi0_recalc;

           delta_ps=alphas_pphi0_recalc-alphas_pphi0_recalc_prev;
		else
		   alphas_pphi0_recalc=(mHe/eV)*(alphas_pos_x+R0).*v_phi_step-(ZHe)*alphas_psi_value_corr;
		end
 		pphi_record=zeros(Nalphas_simulated,1);
      if CALCULATE_DELTAE_POWER==1
            % keep energy change from GC calculations [impossible to use GO position]
%             alphas_Ekin=alphas_Etot-ZHe*alphas_Epot_gc;
%             delta_Ekin=alphas_Ekin-alphas_Ekin_prev;
% 			alphas_Ekin_prev=alphas_Ekin;

			% transfer only total energy change to the wave
            % divide by TIME_GO_SIZE for better Gyro orbit stats
            part2wave_power=-delta_E/TIME_GO_SIZE/h;
			part2wave_power(EJECTED_POP)=0;
            %         alphas_Ekin(INNER_PART_GC)=alphas_Ekin_th(INNER_PART_GC);
            if CALCULATE_VD_POWER==0
                part2wave_vD_power=part2wave_power;
            end
        end
        
        if CALCULATE_VD_POWER==1
            part2wave_vD_power=-(delta_Ekin_vD)/TIME_GO_SIZE/h;	
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
            P_VD_PART_TAE=(NB_PART_RESCALE*eV)*sum(part2wave_vD_power);
            P_PART_TAE=(NB_PART_RESCALE*eV)*sum(part2wave_power);
			
        end
        
        % TAE growth rate
        if USE_VD==0
            gamma_TAE_prev=gamma_TAE_local;
        else
            gamma_TAE_prev=gamma_TAE_vD_local;
        end
        if (TAE_amplitude>=MINIMUM_TAE_AMPLITUDE) && (f_counter>2)
            % only the vD.E part is transfered to the wave ?
			gamma_TAE_vD_local=0.5*(P_VD_PART_TAE/W_TAE);
			gamma_TAE_local=0.5*(P_PART_TAE/W_TAE);

            gamma_TAE_local=max(gamma_TAE_local,gamma_TAE_min);
            gamma_TAE_local=min(gamma_TAE_local,gamma_TAE_max);
        else
            gamma_TAE_local=0;
        end
       

        gamma_TAE=gamma_TAE_local;
        gamma_TAE_vD=gamma_TAE_vD_local;
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
%             Epart_tot_vD_rel_ini=Epart_tot_vD_rel_ini-P_VD_PART_TAE*TIME_GO_SIZE*h;
            Epart_tot_vD_rel_ini=Epart_tot_vD_rel_ini-2*W_TAE*gamma_TAE_vD*TIME_GO_SIZE*h;
            Epart_tot_rel_ini=Epart_tot_rel_ini-2*W_TAE*gamma_TAE*TIME_GO_SIZE*h;
            % Ekin evolution follows gc equation
%             if USE_VD==0
%                 % alphas_Ekin=max(alphas_Ekin-h*TIME_GO_SIZE*part2wave_power,0);
% 				alphas_Ekin=alphas_Ekin+delta_Ekin;
%             if USE_VD==1
%                 alphas_Epot(INNER_PART)=lininterp3(E_potential_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
%                 if ~isempty(OUTER_PART)
%                     alphas_Epot(OUTER_PART)=0;
%                 end
%                 alphas_Ekin=max(alphas_Ekin+delta_Ekin_vD,0);
%                 alphas_Etot=alphas_Ekin+ZHe*alphas_Epot_gc;
%             end
        end
        W_TAE_raw=(TAE_amplitude)^2;
        W_TAE=W_TAE_raw*WTAE_AVG;
%         quasi_lin_flag=1;

        part2wave_power_evol(time_stamp_go)=P_PART_TAE;
        WTAE_evol(time_stamp_go)=W_TAE/WTAE_AVG;
        Ekin_part_evol(time_stamp_go)=Ekin_record/TIME_GO_SIZE;
        Etot_part_evol(time_stamp_go)=Etot_record/TIME_GO_SIZE;
        Etot_th_part_evol(time_stamp_go)=Etot_th_record/TIME_GO_SIZE;
        gamma_TAE_evol(time_stamp_go)=gamma_TAE;
        gamma_TAE_vD_evol(time_stamp_go)=gamma_TAE_vD;
        Epart_tot_vD_rel_evol(time_stamp_go)=Epart_tot_vD_rel_ini;
        Epart_tot_rel_evol(time_stamp_go)=Epart_tot_rel_ini;

        Ekin_record=0;
        Etot_record=0;
        Etot_th_record=0;
        gamma_TAE_record=0;



    end
%     alphas_Epot_prev=alphas_Epot;
%     alphas_Epot_gc_part_prev=alphas_Epot_gc;
