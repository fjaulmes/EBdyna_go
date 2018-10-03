    
    
    % -------------------------------------------------------------------
    % ----------------------- next time step ----------------------------
    % -------------------------------------------------------------------
    v_X=v_X_next;
    v_Z=v_Z_next;
    v_phi=v_phi_next;
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
    % sorting out which particles are inside mixing region and which are out
    interpolate_theta_psi_fromXZ;
    
    % getting field values for this time step and current position	
	Btot_map_phi_prev=Btot_map_phi;
    update_E_b_fields;
	
	% the potentials can actually be used to improve the knowledge of the toroidal Electric field
    calculate_potentials;
    alphas_psi_star=-bphi.*alphas_Apll_value.*alphas_Rpos;
    alphas_psi_value=interp2_XZ(interp_x,interp_z,psi_XZsmall_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
    if ~isempty(OUTER_PART)
        alphas_psi_value_corr(OUTER_PART)=alphas_psi_value(OUTER_PART);		
    end
	if ~isempty(EJECTED_POP)
		alphas_psi_value_corr(EJECTED_POP)=alphas_psi_value(EJECTED_POP);
	end
    alphas_psi_value_corr(INNER_PART)=alphas_psi_value(INNER_PART)+alphas_psi_star(INNER_PART);

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
% 	v_phi_step=0.5*(v_phi+v_phi_next);
    v_phi_step=calc_vstep_value(v_phi_prev_prev,v_phi,v_phi_next);


	
    if CALCULATE_VD_POWER==1
        calculate_gc_pos;
        if ~isempty(OUTER_PART_GC)
            alphas_Bfield_gc_prev(OUTER_PART_GC)=interp2_XZ(interp_x_gc,interp_z_gc,Btot_XZ_map,INDEX_LIST_1_GC,INDEX_LIST_2_GC,INDEX_LIST_3_GC,INDEX_LIST_4_GC,OUTER_PART_GC);
        end
        alphas_Bfield_gc_prev(INNER_PART_GC)=lininterp3(Btot_map_phi_prev,IL3D_1_GC,IL3D_2_GC,IL3D_3_GC,IL3D_4_GC,IL3D_5_GC,IL3D_6_GC,IL3D_7_GC,IL3D_8_GC, slopex_gc,slopey_gc,slopez_gc);
    end
    


    
   if (CALCULATE_DELTAE_POWER==1) %&& (quasi_lin_flag==1)
		if ~isempty(OUTER_PART)
           alphas_dt_Epot(OUTER_PART)=zeros(length(OUTER_PART),1);
           delta_ps(OUTER_PART)=zeros(length(OUTER_PART),1);          
		end
		if ~isempty(OUTER_PART)
			alphas_dt_Epot(EJECTED_POP)=zeros(length(EJECTED_POP),1);
			delta_ps(EJECTED_POP)=zeros(length(EJECTED_POP),1);
			dE_prev(EJECTED_POP)=zeros(length(EJECTED_POP),1);
			delta_E_th(EJECTED_POP)=zeros(length(EJECTED_POP),1);
			delta_E(EJECTED_POP)=zeros(length(EJECTED_POP),1);
		end
        % the division by R is not necessary (we did not multiply)
%         delta_ps(INNER_PART)=ZHe*h*((alphas_dt_Epot(INNER_PART))-(v_phi_step(INNER_PART)./bphi(INNER_PART)).*alphas_dt_psi_star(INNER_PART));
        delta_ps(INNER_PART)=-(ZHe*h*omega_TAE)*(alphas_iEpot(INNER_PART)-(v_phi_step(INNER_PART)./bphi(INNER_PART)).*alphas_iA(INNER_PART));
		%dE=delta_ps;
		dE=0.5*(delta_ps+dE_prev);
		dE_prev=delta_ps;
        if USE_DELTAE_TH==1
%             alphas_Etot_th(INNER_PART_GC)=alphas_Etot_th(INNER_PART_GC)+delta_ps(INNER_PART_GC);
            delta_E_th(INNER_PART)=delta_E_th(INNER_PART)+dE(INNER_PART);          
            alphas_Etot_th(INNER_PART)=alphas_Etot_th_prev(INNER_PART)+delta_E_th(INNER_PART);
			
            %delta_E_th=delta_ps;
        else
			if USE_DELTAE_EVOL==1
                delta_E(INNER_PART)=delta_E(INNER_PART)+dE(INNER_PART);
				alphas_Etot(INNER_PART)=alphas_Etot_prev(INNER_PART)+delta_E(INNER_PART);
			end
			%alphas_Etot=double(alphas_Etot);
%             alphas_Etot(INNER_PART_GC)=alphas_Etot(INNER_PART_GC)+delta_ps(INNER_PART_GC);
%             delta_E(INNER_PART_GC)=delta_E(INNER_PART_GC)+delta_ps(INNER_PART_GC);
%             disp(sum(delta_ps(INNER_PART_GC)))
        end
    end

    % amplitude of the B field and toroidal gradients at local positions
    %if ~isempty(OUTER_PART)
    %    if CALCULATE_TRUE_PPHI==1
    %        alphas_grad_Phi(OUTER_PART)=zeros(length(OUTER_PART),1);
    %        alphas_grad_psi_star(OUTER_PART)=zeros(length(OUTER_PART),1);
	%    end
    %end
    %if ~isempty(EJECTED_POP)
    %    if CALCULATE_TRUE_PPHI==1
    %        alphas_grad_Phi(EJECTED_POP)=zeros(length(EJECTED_POP),1);
    %        alphas_grad_psi_star(EJECTED_POP)=zeros(length(EJECTED_POP),1);
	%    end
    %end   

    
    if CALCULATE_VD_POWER==1
        calculate_vD_E_gc_components;
    end
	




    % particles -> wave interaction

    part2wave_power(OUTER_PART)=0;
	if ~isempty(EJECTED_POP)
		part2wave_power(EJECTED_POP)=0;
		part2wave_vD_power(EJECTED_POP)=0;
	end
    if (CALCULATE_VD_POWER==1)

        delta_B(INNER_PART_GC)=(alphas_Bfield_gc(INNER_PART_GC)-alphas_Bfield_gc_prev(INNER_PART_GC));
        part2wave_vD_power(OUTER_PART_GC)=0;
        part2wave_vD_power(INNER_PART_GC)=-ZHe*vD0(INNER_PART_GC).*(vDX(INNER_PART_GC).*EX_gc(INNER_PART_GC)+vDZ(INNER_PART_GC).*EZ_gc(INNER_PART_GC)+vDphi(INNER_PART_GC).*Ephi_gc(INNER_PART_GC));
        part2wave_vD_power(INNER_PART_GC)=part2wave_vD_power(INNER_PART_GC)-(alphas_mm_gc(INNER_PART_GC).*delta_B(INNER_PART_GC))/h;

        delta_Ekin_vD(INNER_PART_GC)=delta_Ekin_vD(INNER_PART_GC)-part2wave_vD_power(INNER_PART_GC)*h;

        if CALCULATE_DELTAE_POWER==0
            part2wave_power(INNER_PART_GC)=part2wave_vD_power(INNER_PART_GC);
        end
    end



    
    % Canonical angular momentum evolution
	alphas_pphi_prev=alphas_pphi0;
 	if ~isempty(EJECTED_POP)
		delta_pphi_prev(EJECTED_POP)=0;
		delta_pphi_raw(EJECTED_POP)=0;
		delta_pphi(EJECTED_POP)=0;
		delta_pphi_th(EJECTED_POP)=0;
	end   
    if (CALCULATE_TRUE_PPHI==1)        
        % for single eigen omega
%         alphas_vphi_grad_psi_star(INNER_PART)=(v_phi_step(INNER_PART).*bphi(INNER_PART).*alphas_grad_psi_star(INNER_PART));
        
        delta_pphi_raw(INNER_PART)=-(h*ZHe*nTAE)*((alphas_iEpot(INNER_PART))-v_phi_step(INNER_PART).*bphi(INNER_PART).*alphas_iA(INNER_PART));
%         delta_pphi2(INNER_PART)=-(h*ZHe)*((alphas_grad_Phi(INNER_PART))+alphas_vphi_grad_psi_star(INNER_PART));
        delta_pphi=0.5*(delta_pphi_raw+delta_pphi_prev);
		delta_pphi_prev=delta_pphi_raw;

        % this is an half next time step pphi
		alphas_pphi0(INNER_PART)=alphas_pphi0(INNER_PART)+delta_pphi(INNER_PART);
		
        
        % compromise for precision
%         alphas_pphi0(INNER_PART)=0.99999*alphas_pphi0(INNER_PART)+0.00001*((mHe/eV)*(alphas_Rpos(INNER_PART)).*v_phi_step(INNER_PART)-ZHe*alphas_psi_value_corr(INNER_PART));
    else
        alphas_pphi0(INNER_PART)=((mHe/eV)*(alphas_Rpos(INNER_PART)).*v_phi_step(INNER_PART)-ZHe*alphas_psi_value_corr(INNER_PART));
    end
    delta_pphi_th(INNER_PART)=delta_pphi_th(INNER_PART)+delta_pphi(INNER_PART);

    % Total energy evolution
	if ~isempty(EJECTED_POP)
		delta_E(EJECTED_POP)=0;
		delta_E_th(EJECTED_POP)=0;
		delta_ps(EJECTED_POP)=0;
	end
    if (CALCULATE_DELTAE_POWER==1) % && (quasi_lin_flag==1)
        % pphi evolution
        delta_ps(INNER_PART)=alphas_iEpot(INNER_PART)-(v_phi_step(INNER_PART).*bphi(INNER_PART)).*alphas_iA(INNER_PART);
        delta_ps(OUTER_PART)=zeros(length(OUTER_PART),1);
%         
        INNER_TH_PART=find(NON_EJECTED_POP.*(abs(delta_ps)>MIN_DIV_DPPHI));
        INNER_BLOWUP_PART=find(NON_EJECTED_POP.*(abs(delta_ps)<=MIN_DIV_DPPHI).*(alphas_psi<=(pTAE_sup-1)).*(alphas_psi>=(pTAE_inf+1)));
%         if ~isempty(INNER_BLOWUP_PART)
%             disp('found some resonant crazy particle')
%         end

        % energy evolution
        delta_ps(INNER_TH_PART)=(omega_TAE/nTAE)*(alphas_iEpot(INNER_TH_PART)-(v_phi_step(INNER_TH_PART)./bphi(INNER_TH_PART)).*alphas_iA(INNER_TH_PART))./(alphas_iEpot(INNER_TH_PART)-(v_phi_step(INNER_TH_PART).*bphi(INNER_TH_PART)).*alphas_iA(INNER_TH_PART));
        delta_ps(INNER_BLOWUP_PART)=(omega_TAE/nTAE)*(alphas_iEpot(INNER_BLOWUP_PART)-(v_phi_step(INNER_BLOWUP_PART)./bphi(INNER_BLOWUP_PART)).*alphas_iA(INNER_BLOWUP_PART))./(MIN_DIV_DPPHI);
        delta_ps(INNER_PART)=delta_ps(INNER_PART).*delta_pphi(INNER_PART);
		delta_ps(OUTER_PART)=zeros(length(OUTER_PART),1);
		%deltaE_th=delta_ps;
        
%         delta_ps(INNER_TH_PART)=(omega_TAE/nTAE)*((1-(kTAE/omega_TAE)*v_phi_step(INNER_TH_PART)./bphi(INNER_TH_PART))./(1-(kTAE/omega_TAE)*v_phi_step(INNER_TH_PART).*bphi(INNER_TH_PART))).*delta_pphi(INNER_TH_PART);
%         % some particles just can't take this lovely division
%         % we give them a moderate kick from approximate expression
%         delta_ps(INNER_BLOWUP_PART)=(omega_TAE/nTAE)*((1-(kTAE/omega_TAE)*v_phi_step(INNER_BLOWUP_PART)./bphi(INNER_BLOWUP_PART))/MIN_DIV_DPPHI).*delta_pphi(INNER_BLOWUP_PART);
        if USE_DELTAE_TH==0
%             alphas_Etot_th(INNER_PART)=alphas_Etot_th(INNER_PART)+delta_ps(INNER_PART);
             delta_E_th(INNER_PART)=delta_E_th(INNER_PART)+delta_ps(INNER_PART);          
             alphas_Etot_th(INNER_PART)=alphas_Etot_th_prev(INNER_PART)+delta_E_th(INNER_PART);
			 %alphas_Etot_th=double(alphas_Etot_th);
        else
        % compromise for precision
%             alphas_Etot=0.5*(alphas_Etot)+0.5*(alphas_Etot_th);
%             alphas_Etot(INNER_PART)=0.5*(alphas_Etot(INNER_PART)+alphas_Etot_th(INNER_PART));
			if USE_DELTAE_EVOL==1
                delta_E(INNER_PART)=delta_E(INNER_PART)+delta_ps(INNER_PART);			
				alphas_Etot(INNER_PART)=alphas_Etot_prev(INNER_PART)+delta_E(INNER_PART);
			end
        end
        %         alphas_Ekin_th(INNER_PART_GC)=max(alphas_Etot(INNER_PART_GC)-ZHe*(alphas_Epot_gc(INNER_PART_GC)),0);
    end

    % finally adjust Ekin to make everybody happy
    % compromise between GC evol and instant particle picture
	if USE_DELTAE_EVOL==1
		if USE_VD==0
	%         alphas_Epot_gc(INNER_PART_GC)=lininterp3(E_potential_map_phi,IL3D_1_GC,IL3D_2_GC,IL3D_3_GC,IL3D_4_GC,IL3D_5_GC,IL3D_6_GC,IL3D_7_GC,IL3D_8_GC, slopex_gc,slopey_gc,slopez_gc);
	%         if ~isempty(OUTER_PART_GC)
	%             alphas_Epot_gc(OUTER_PART_GC)=zeros(length(OUTER_PART_GC),1);
	%         end
		   alphas_Ekin=max(alphas_Etot-ZHe*alphas_Epot,0);
	%        alphas_Ekin=0.5*(alphas_Ekin+max(alphas_Etot_th-ZHe*alphas_Epot_gc,0));
		else
			alphas_Ekin=alphas_Ekin-part2wave_vD_power*h;
		end
    else
	    alphas_Ekin=0.5*(mHe/eV)*(v_X_step.^2+v_Z_step.^2+v_phi_step.^2);
	    alphas_Etot=alphas_Ekin+ZHe*alphas_Epot;
		delta_E=alphas_Etot-alphas_Etot_prev;
	end
    if (quasi_lin_flag==0)
        %after change of amplitude, do not evolve total energy
        quasi_lin_flag=1;
    end
	
%	gamma_TAE_record=gamma_TAE_record+gamma_TAE;
%     WTAE_record(round(0.5*record_step))=W_TAE;
%     Etot_th_record=Etot_th_record+NB_PART_RESCALE*eV*sum(alphas_weight.*alphas_Etot_th);
%     Ekin_record=Ekin_record+NB_PART_RESCALE*eV*sum(alphas_weight.*alphas_Ekin);
%     Etot_record=Etot_record+NB_PART_RESCALE*eV*sum(alphas_weight.*alphas_Etot);
	

    % pphi_record=pphi_record+((mHe/eV)*(alphas_Rpos).*v_phi_step-ZHe*alphas_psi_value_corr);
    % pphi_record=pphi_record+alphas_pphi0;
    
    
v_phi_step_recalc=(alphas_pphi0+(ZHe)*alphas_psi_value_corr)./(alphas_Rpos);
v_phi_step_recalc=(eV/mHe)*v_phi_step_recalc;
	
delta_vphi=PPHI_CORR_FACTOR*(v_phi_step_recalc-v_phi_step);
v_phi_corr_integ=v_phi_corr_integ+PPHI_CORR_INTEG_FACTOR*delta_vphi;


    
adapt_speed_Ekin_G;
v_phi_next=v_phi_next+delta_vphi+v_phi_corr_integ;
% % 
%adapt_speed_Ekin_G;
%v_phi_next=v_phi_next+delta_vphi+v_phi_corr_integ;
% 
%adapt_speed_Ekin_G;
%v_phi_next=v_phi_next+delta_vphi+v_phi_corr_integ;
    