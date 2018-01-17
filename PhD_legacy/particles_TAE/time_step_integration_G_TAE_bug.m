    

    % -------------------------------------------------------------------
    % ----------------------- next time step ----------------------------
    % -------------------------------------------------------------------

    time=time+DELTA_TIME;
    % recording values
    v_X_prev=v_X;
    v_Z_prev=v_Z;
    v_phi_prev_prev=v_phi_prev;
    v_phi_prev=v_phi;
	
    
    %  update particle position tokamak_coordinates;
    alphas_Rpos=alphas_pos_x+R0;
    alphas_pos_x=alphas_pos_x+h*v_X;
    alphas_Rpos_int=0.5*(alphas_Rpos+alphas_pos_x+R0);
    alphas_pos_phi=alphas_pos_phi_wrapped+h*(v_phi./alphas_Rpos_int);
    alphas_pos_z=alphas_pos_z+h*v_Z;
    alphas_Rpos=alphas_pos_x+R0;
    alphas_pos_phi_wrapped=wrap2pi(alphas_pos_phi);
    alphas_pos_phi_TAE=wrapTAEangle(alphas_pos_phi_wrapped,TAE_ANGLE);
    
    % limited map range
    alphas_pos_x=min(alphas_pos_x,scale_X(end-3));
    alphas_pos_z=min(alphas_pos_z,scale_Z(end-3));
    alphas_pos_x=max(alphas_pos_x,scale_X(4));
    alphas_pos_z=max(alphas_pos_z,scale_Z(4));

    % estimating psi and theta values
    interpolate_theta_psi_fromXZ;
    % sorting out which particles are inside mixing region and which are out
	% clear INNER_PART OUTER_PART;
    INNER_PART=find((alphas_psi<=(pTAE_sup-1)).*(alphas_psi>=(pTAE_inf+1)));
    OUTER_PART=find((alphas_psi>(pTAE_sup-1))+(alphas_psi<(pTAE_inf+1))>0);
    [IL3D_1 IL3D_2 IL3D_3 IL3D_4 IL3D_5 IL3D_6 IL3D_7 IL3D_8 slopex slopey slopez] = ...
        build_3Dinterp_indexarrays(scale_phi, scale_theta, scale_psi, DPHI,DTHETA,1,alphas_pos_phi_TAE(INNER_PART), alphas_theta(INNER_PART), alphas_psi(INNER_PART));
    
    % getting field values for this time step and current position	
	Btot_map_phi_prev=Btot_map_phi;
    E_potential_map_phi_prev=E_potential_map_phi;
    update_E_b_fields;
    update_Gfields_TAE;

 
	%then deriving immediately the next speed value
    update_G_3D_TAE;
    v_X_next=v_X;
    v_Z_next=v_Z;
    v_phi_next=v_phi;
	v_X=v_X_prev;
	v_Z=v_Z_prev;
	v_phi=v_phi_prev;
	
	% now we have velocity values at step position
	v_X_step=0.5*(v_X+v_X_next);
	v_Z_step=0.5*(v_Z+v_Z_next);
	v_phi_step=0.5*(v_phi+v_phi_next);


    alphas_Epot_gc_part_prev=alphas_Epot_gc;

    alphas_Epot(INNER_PART)=lininterp3(E_potential_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);

    calculate_gc_pos;
	
    if CALCULATE_VD_POWER==1
        if ~isempty(OUTER_PART_GC)
            alphas_Bfield_gc_prev(OUTER_PART_GC)=interp2_XZ(interp_x_gc,interp_z_gc,Btot_XZ_map,INDEX_LIST_1_GC,INDEX_LIST_2_GC,INDEX_LIST_3_GC,INDEX_LIST_4_GC,OUTER_PART_GC);
        end
        alphas_Bfield_gc_prev(INNER_PART_GC)=lininterp3(Btot_map_phi_prev,IL3D_1_GC,IL3D_2_GC,IL3D_3_GC,IL3D_4_GC,IL3D_5_GC,IL3D_6_GC,IL3D_7_GC,IL3D_8_GC, slopex_gc,slopey_gc,slopez_gc);
    end

    
    if CALCULATE_DELTAE_POWER==1
        % this is Epot prev at the NEW location
        alphas_Epot_gc_prev(INNER_PART_GC)=lininterp3( E_potential_map_phi_prev,IL3D_1_GC,IL3D_2_GC,IL3D_3_GC,IL3D_4_GC,IL3D_5_GC,IL3D_6_GC,IL3D_7_GC,IL3D_8_GC, slopex_gc,slopey_gc,slopez_gc);
        if ~isempty(OUTER_PART_GC)
            % %         alphas_psi_star_prev(OUTER_PART)=zeros(length(OUTER_PART),1);
            alphas_Epot_gc_prev(OUTER_PART_GC)=zeros(length(OUTER_PART_GC),1);
            alphas_Epot_gc(OUTER_PART_GC)=zeros(length(OUTER_PART_GC),1);
        end
        alphas_psi_star_prev=-(kTAE/omega_TAE)*alphas_Epot_gc_prev.*alphas_Rpos_gc;
        
        alphas_Epot_gc(INNER_PART_GC)=lininterp3(E_potential_map_phi,IL3D_1_GC,IL3D_2_GC,IL3D_3_GC,IL3D_4_GC,IL3D_5_GC,IL3D_6_GC,IL3D_7_GC,IL3D_8_GC, slopex_gc,slopey_gc,slopez_gc);
        
        if ~isempty(OUTER_PART_GC)
            alphas_Epot_gc(OUTER_PART_GC)=zeros(length(OUTER_PART_GC),1);
        end
        alphas_psi_star=-(kTAE/omega_TAE)*alphas_Epot_gc.*alphas_Rpos_gc;
    end
    
    % amplitude of the B field and toroidal gradients at local positions
    if ~isempty(OUTER_PART)
        if CALCULATE_TRUE_PPHI==1
            alphas_grad_Phi(OUTER_PART)=zeros(length(OUTER_PART),1);
            alphas_grad_psi_star(OUTER_PART)=zeros(length(OUTER_PART),1);
	    end
    end
    if CALCULATE_TRUE_PPHI==1
        alphas_grad_psi_star(INNER_PART)=lininterp3( grad_psi_star_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
        % no adjustment with R necessary since we will divide by R this value
        % for pphi calculation
        %         alphas_grad_psi_star(INNER_PART)=alphas_grad_psi_star(INNER_PART).*alphas_Rpos(INNER_PART);
    end
    

    
    if CALCULATE_VD_POWER==1
        calculate_vD_E_gc_components;
    end
	

    % Total energy evolution
    if CALCULATE_DELTAE_POWER==1
        v_phi_gc=v_phi.*(bphi_gc./bphi);
        % psi star evolution
        delta_ps(INNER_PART_GC)=(alphas_psi_star(INNER_PART_GC)-alphas_psi_star_prev(INNER_PART_GC));
        delta_Epot(INNER_PART_GC)=(alphas_Epot_gc(INNER_PART_GC)-alphas_Epot_gc_part_prev(INNER_PART_GC));
       % Kinetic energy evolution
        alphas_Ekin_prev=alphas_Ekin;
        delta_E(INNER_PART_GC)=ZHe*((alphas_Epot_gc(INNER_PART_GC)-alphas_Epot_gc_prev(INNER_PART_GC))+(v_phi_gc(INNER_PART_GC)./alphas_Rpos_gc(INNER_PART_GC)).*delta_ps(INNER_PART_GC));
        alphas_Etot(INNER_PART_GC)=alphas_Etot(INNER_PART_GC)+delta_E(INNER_PART_GC);
        alphas_Ekin_th(INNER_PART_GC)=max(alphas_Etot(INNER_PART_GC)-ZHe*(alphas_Epot_gc(INNER_PART_GC)),0);
    end



    % particles -> wave interaction

    part2wave_power(OUTER_PART_GC)=0;
    if (CALCULATE_VD_POWER==1)

        delta_B(INNER_PART_GC)=(alphas_Bfield_gc(INNER_PART_GC)-alphas_Bfield_gc_prev(INNER_PART_GC));
        part2wave_vD_power(OUTER_PART_GC)=0;
        part2wave_vD_power(INNER_PART_GC)=-ZHe*vD0(INNER_PART_GC).*(vDX(INNER_PART_GC).*EX_gc(INNER_PART_GC)+vDZ(INNER_PART_GC).*EZ_gc(INNER_PART_GC)+vDphi(INNER_PART_GC).*Ephi_gc(INNER_PART_GC));
        part2wave_vD_power(INNER_PART_GC)=part2wave_vD_power(INNER_PART_GC)-(alphas_mm_gc(INNER_PART_GC).*delta_B(INNER_PART_GC))/h;
        if USE_VD==1
            Pexchange_part_record(INNER_PART_GC)=Pexchange_part_record(INNER_PART_GC)+(part2wave_vD_power(INNER_PART_GC));
		end

        %         disp('STOP')
        %         pause

        if CALCULATE_DELTAE_POWER==0
            part2wave_power(INNER_PART_GC)=part2wave_vD_power(INNER_PART_GC);
        end
    end
	if CALCULATE_DELTAE_POWER==1
        delta_Ekin(OUTER_PART_GC)=0;
        delta_Ekin(INNER_PART_GC)=delta_E(INNER_PART_GC)-ZHe*delta_Epot(INNER_PART_GC);
        part2wave_power(INNER_PART_GC)=-delta_Ekin(INNER_PART_GC)/h;
        %         alphas_Ekin(INNER_PART_GC)=alphas_Ekin_th(INNER_PART_GC);
		if USE_VD==0
		    Pexchange_part_record(INNER_PART_GC)=Pexchange_part_record(INNER_PART_GC)+(part2wave_power(INNER_PART_GC));
        end
		if CALCULATE_VD_POWER==0
            part2wave_vD_power(INNER_PART_GC)=part2wave_power(INNER_PART_GC);
        end
	end
	
    if (f_counter<2)
        % no drive initially
        P_VD_PART_TAE=0;
        P_PART_TAE=0;
    else
        P_VD_PART_TAE=(NB_PART_RESCALE*eV)*sum(part2wave_vD_power(INNER_PART_GC));
        P_PART_TAE=(NB_PART_RESCALE*eV)*sum(part2wave_power(INNER_PART_GC));
   end

    % TAE growth rate from power transfer and mode energy
	
    gamma_TAE_prev=gamma_TAE;
    if (TAE_amplitude>=MINIMUM_TAE_AMPLITUDE) && (f_counter>=2)
        % only the vD.E part is transfered to the wave ?
        if USE_VD==0
            gamma_TAE=0.5*(P_PART_TAE/W_TAE);
        else
            gamma_TAE=0.5*(P_VD_PART_TAE/W_TAE);
        end       
        gamma_TAE=max(gamma_TAE,gamma_TAE_min);
        gamma_TAE=min(gamma_TAE,gamma_TAE_max);
    else
        gamma_TAE=0;
    end
%     gamme_TAE_half=max(0.5*(gamme_TAE+gamme_TAE_prev),0);
    gamma_TAE_half=0.5*(gamma_TAE+gamma_TAE_prev);
    if (gamma_TAE_prev~=0)
        if EVOLVE_TAE_AMPLITUDE==1
            if (TAE_amplitude>=MINIMUM_TAE_AMPLITUDE)
                TAE_amplitude=max(TAE_amplitude+TAE_amplitude*gamma_TAE_half*h,0);
            else
                TAE_amplitude=max(TAE_amplitude+TAE_amplitude*gamma_TAE_half*h,MINIMUM_TAE_AMPLITUDE);
            end
        end
        Epart_tot_vD_rel_ini=Epart_tot_vD_rel_ini-P_VD_PART_TAE*h;
        Epart_tot_rel_ini=Epart_tot_rel_ini-P_PART_TAE*h;
        % Ekin evolution follows gc equation
        if USE_VD==0
            alphas_Ekin(INNER_PART_GC)=max(alphas_Ekin(INNER_PART_GC)-h*part2wave_power(INNER_PART_GC),0);
        else
            alphas_Ekin(INNER_PART_GC)=max(alphas_Ekin(INNER_PART_GC)-h*part2wave_vD_power(INNER_PART_GC),0);
        end
    end
    W_TAE=(TAE_amplitude)^2;    

    
    % Canonical angular momentum evolution
    if (CALCULATE_TRUE_PPHI==1)        
        alphas_vphi_grad_psi_star(INNER_PART)=(v_phi_step(INNER_PART).*alphas_grad_psi_star(INNER_PART));
		alphas_grad_Phi(INNER_PART)=-(omega_TAE/kTAE)*alphas_grad_psi_star(INNER_PART);
		alphas_pphi0(INNER_PART)=alphas_pphi0(INNER_PART)-...
			(h*ZHe)*((alphas_grad_Phi(INNER_PART))+alphas_vphi_grad_psi_star(INNER_PART));
    else
        alphas_psi_value=interp2_XZ(interp_x,interp_z,psi_XZsmall_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
        if ~isempty(OUTER_PART)
            alphas_psi_value_corr(OUTER_PART)=alphas_psi_value(OUTER_PART);
        end
        alphas_psi_value_corr(INNER_PART)=alphas_psi_value(INNER_PART)+alphas_psi_star(INNER_PART);
        v_phi_step=0.5*(v_phi+v_phi_prev);
        alphas_pphi0=mHe*(alphas_pos_x+R0).*v_phi_step-ZHe*alphas_psi_value_corr;
	end   
	
	
	gamma_TAE_record=gamma_TAE_record+gamma_TAE;
%     WTAE_record(round(0.5*record_step))=W_TAE;
    Ekin_record=Ekin_record+NB_PART_RESCALE*eV*sum(alphas_Ekin);
    Etot_record=Etot_record+NB_PART_RESCALE*eV*sum(alphas_Etot);
    pphi_record=pphi_record+alphas_pphi0;