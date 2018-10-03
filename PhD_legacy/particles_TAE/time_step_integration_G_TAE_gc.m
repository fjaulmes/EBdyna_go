    

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
    
    %dunno why we cannot reuse previous iteration next speed values
%     update_G_3D_TAE;
    
	
    
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
    iEpot_map_phi_int=iEpot_map_phi;
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





    calculate_gc_pos;
	
    if CALCULATE_VD_POWER==1
        if ~isempty(OUTER_PART_GC)
            alphas_Bfield_gc_prev(OUTER_PART_GC)=interp2_XZ(interp_x_gc,interp_z_gc,Btot_XZ_map,INDEX_LIST_1_GC,INDEX_LIST_2_GC,INDEX_LIST_3_GC,INDEX_LIST_4_GC,OUTER_PART_GC);
        end
        alphas_Bfield_gc_prev(INNER_PART_GC)=lininterp3(Btot_map_phi_prev,IL3D_1_GC,IL3D_2_GC,IL3D_3_GC,IL3D_4_GC,IL3D_5_GC,IL3D_6_GC,IL3D_7_GC,IL3D_8_GC, slopex_gc,slopey_gc,slopez_gc);
    end

    
   if (CALCULATE_DELTAE_POWER==1) %&& (quasi_lin_flag==1)
        iEpot_map_phi_int=0.5*(iEpot_map_phi_int+iEpot_map_phi);
        %  Epot prev is partial dEpot/dt at the NEW location
        alphas_Epot_gc_prev(INNER_PART_GC)=-(h*omega_TAE/nTAE)*lininterp3(iEpot_map_phi_int,IL3D_1_GC,IL3D_2_GC,IL3D_3_GC,IL3D_4_GC,IL3D_5_GC,IL3D_6_GC,IL3D_7_GC,IL3D_8_GC, slopex_gc,slopey_gc,slopez_gc);
        delta_ps(INNER_PART_GC)=-(kTAE/omega_TAE)*bphi_gc(INNER_PART_GC).*alphas_Epot_gc_prev(INNER_PART_GC);
%        % this is Epot prev at the NEW location
%        alphas_Epot_gc_prev(INNER_PART_GC)=lininterp3( E_potential_map_phi_prev,IL3D_1_GC,IL3D_2_GC,IL3D_3_GC,IL3D_4_GC,IL3D_5_GC,IL3D_6_GC,IL3D_7_GC,IL3D_8_GC, slopex_gc,slopey_gc,slopez_gc);
       if ~isempty(OUTER_PART_GC)
           % %         alphas_psi_star_prev(OUTER_PART)=zeros(length(OUTER_PART),1);
           alphas_Epot_gc_prev(OUTER_PART_GC)=zeros(length(OUTER_PART_GC),1);
           delta_ps(OUTER_PART_GC)=zeros(length(OUTER_PART_GC),1);
       end
       % the multiplication by Rpos is not useful here
%        alphas_psi_star_prev=-(kTAE/omega_TAE)*bphi_gc.*alphas_Epot_gc_prev;
%        alphas_Epot_gc(INNER_PART_GC)=lininterp3(E_potential_map_phi,IL3D_1_GC,IL3D_2_GC,IL3D_3_GC,IL3D_4_GC,IL3D_5_GC,IL3D_6_GC,IL3D_7_GC,IL3D_8_GC, slopex_gc,slopey_gc,slopez_gc);
%        if ~isempty(OUTER_PART_GC)
%            alphas_Epot_gc(OUTER_PART_GC)=zeros(length(OUTER_PART_GC),1);
%        end
        % the multiplication by Rpos is not useful here
%        alphas_psi_star=-(kTAE/omega_TAE)*bphi_gc.*alphas_Epot_gc;
%    end
%    if (quasi_lin_flag==0)
%        %after change of amplitude, do not evolve total energy
%        alphas_Epot_gc(INNER_PART_GC)=lininterp3(E_potential_map_phi,IL3D_1_GC,IL3D_2_GC,IL3D_3_GC,IL3D_4_GC,IL3D_5_GC,IL3D_6_GC,IL3D_7_GC,IL3D_8_GC, slopex_gc,slopey_gc,slopez_gc);
%        if ~isempty(OUTER_PART_GC)
%            alphas_Epot_gc(OUTER_PART_GC)=zeros(length(OUTER_PART_GC),1);
%        end
%        alphas_Epot_gc_prev=alphas_Epot_gc;
%         % the multiplication by Rpos is not useful here
%        alphas_psi_star=-(kTAE/omega_TAE)*bphi_gc.*alphas_Epot_gc;
%        alphas_psi_star_prev=alphas_psi_star;
%    end
%     if (CALCULATE_DELTAE_POWER==1) % && (quasi_lin_flag==1)
        v_phi_gc=v_phi.*(bphi_gc./bphi);
        % psi star evolution
%         delta_ps(INNER_PART_GC)=(alphas_psi_star(INNER_PART_GC)-alphas_psi_star_prev(INNER_PART_GC));
%         delta_Epot(INNER_PART_GC)=delta_Epot(INNER_PART_GC)+(alphas_Epot_gc(INNER_PART_GC)-alphas_Epot_gc_part_prev(INNER_PART_GC));
       % Kinetic energy evolution
        alphas_Ekin_prev=alphas_Ekin;
        % the division by R is not necessary (we did not multiply)
        delta_Epot(INNER_PART_GC)=ZHe*((alphas_Epot_gc_prev(INNER_PART_GC))+(v_phi_gc(INNER_PART_GC)).*delta_ps(INNER_PART_GC));
        if USE_DELTAE_TH==1
            alphas_Etot_th(INNER_PART)=alphas_Etot_th(INNER_PART)+delta_Epot(INNER_PART);
        else
            alphas_Etot(INNER_PART)=alphas_Etot(INNER_PART)+delta_Epot(INNER_PART);
            delta_E(INNER_PART_GC)=delta_E(INNER_PART_GC)+delta_Epot(INNER_PART_GC);
        end
    end

    % amplitude of the B field and toroidal gradients at local positions
    if ~isempty(OUTER_PART)
        if CALCULATE_TRUE_PPHI==1
            alphas_grad_Phi(OUTER_PART)=zeros(length(OUTER_PART),1);
            alphas_grad_psi_star(OUTER_PART)=zeros(length(OUTER_PART),1);
	    end
    end
    if CALCULATE_TRUE_PPHI==1
        alphas_grad_Phi(INNER_PART)=lininterp3( iEpot_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
        % no adjustment with R necessary since we will divide by R this value
        % for pphi calculation
        %         alphas_grad_psi_star(INNER_PART)=alphas_grad_psi_star(INNER_PART).*alphas_Rpos(INNER_PART);
    end
    

    
    if CALCULATE_VD_POWER==1
        calculate_vD_E_gc_components;
    end
	




    % particles -> wave interaction

    part2wave_power(OUTER_PART_GC)=0;
    if (CALCULATE_VD_POWER==1)

        delta_B(INNER_PART_GC)=(alphas_Bfield_gc(INNER_PART_GC)-alphas_Bfield_gc_prev(INNER_PART_GC));
        part2wave_vD_power(OUTER_PART_GC)=0;
        part2wave_vD_power(INNER_PART_GC)=-ZHe*vD0(INNER_PART_GC).*(vDX(INNER_PART_GC).*EX_gc(INNER_PART_GC)+vDZ(INNER_PART_GC).*EZ_gc(INNER_PART_GC)+vDphi(INNER_PART_GC).*Ephi_gc(INNER_PART_GC));
        part2wave_vD_power(INNER_PART_GC)=part2wave_vD_power(INNER_PART_GC)-(alphas_mm_gc(INNER_PART_GC).*delta_B(INNER_PART_GC))/h;

        delta_Ekin_vD(INNER_PART_GC)=delta_Ekin_vD(INNER_PART_GC)-part2wave_vD_power(INNER_PART_GC)*h;

        %         disp('STOP')
        %         pause

        if CALCULATE_DELTAE_POWER==0
            part2wave_power(INNER_PART_GC)=part2wave_vD_power(INNER_PART_GC);
        end
    end



    
    % Canonical angular momentum evolution
	alphas_pphi0_prev=alphas_pphi0;
    if (CALCULATE_TRUE_PPHI==1)        
        % for single omega
		alphas_grad_psi_star(INNER_PART)=-(kTAE/omega_TAE)*alphas_grad_Phi(INNER_PART);
        alphas_vphi_grad_psi_star(INNER_PART)=(v_phi_step(INNER_PART).*bphi(INNER_PART).*alphas_grad_psi_star(INNER_PART));
% 		alphas_grad_Phi(INNER_PART)=-(omega_TAE/kTAE)*alphas_grad_psi_star(INNER_PART);
        % this is an half next time step pphi
		alphas_pphi0(INNER_PART)=alphas_pphi0(INNER_PART)-...
			(h*ZHe)*((alphas_grad_Phi(INNER_PART))+alphas_vphi_grad_psi_star(INNER_PART));
        alphas_pphi0=0.5*(alphas_pphi0+alphas_pphi0_prev);
%     else
        alphas_Epot(INNER_PART)=lininterp3(E_potential_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
        if ~isempty(OUTER_PART)
            alphas_Epot(OUTER_PART)=0;
        end
%        alphas_Epot_gc(INNER_PART_GC)=lininterp3(E_potential_map_phi,IL3D_1_GC,IL3D_2_GC,IL3D_3_GC,IL3D_4_GC,IL3D_5_GC,IL3D_6_GC,IL3D_7_GC,IL3D_8_GC, slopex_gc,slopey_gc,slopez_gc);
%        if ~isempty(OUTER_PART_GC)
%            alphas_Epot_gc(OUTER_PART_GC)=zeros(length(OUTER_PART_GC),1);
%        end
        alphas_psi_star=-(kTAE/omega_TAE)*bphi.*alphas_Epot.*alphas_Rpos;
        alphas_psi_value=interp2_XZ(interp_x,interp_z,psi_XZsmall_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
        if ~isempty(OUTER_PART)
            alphas_psi_value_corr(OUTER_PART)=alphas_psi_value(OUTER_PART);
        end
        alphas_psi_value_corr(INNER_PART)=alphas_psi_value(INNER_PART)+alphas_psi_star(INNER_PART);
        % compromise for precision
        alphas_pphi0(INNER_PART)=0.9999*alphas_pphi0(INNER_PART)+0.0001*((mHe/eV)*(alphas_Rpos(INNER_PART)).*v_phi_step(INNER_PART)-ZHe*alphas_psi_value_corr(INNER_PART));
    end

    % Total energy evolution
    if (CALCULATE_DELTAE_POWER==1) % && (quasi_lin_flag==1)
        %v_phi_gc=v_phi.*(bphi_gc./bphi);
        % psi star evolution
        delta_ps(INNER_PART)=1-(kTAE/omega_TAE)*v_phi(INNER_PART).*bphi(INNER_PART);
        delta_ps(OUTER_PART)=0;
%         INNER_TH_PART=find((abs(delta_ps)>MIN_DIV_DPPHI).*(alphas_psi<=(pTAE_sup-1)).*(alphas_psi>=(pTAE_inf+1)));
        INNER_TH_PART=find(abs(delta_ps)>MIN_DIV_DPPHI);
        INNER_BLOWUP_PART=find((abs(delta_ps)<=MIN_DIV_DPPHI).*(alphas_psi<=(pTAE_sup-1)).*(alphas_psi>=(pTAE_inf+1)));

        % energy evolution
        delta_ps=alphas_pphi0-alphas_pphi0_prev;
        delta_ps(INNER_TH_PART)=(omega_TAE/nTAE)*((1-(kTAE/omega_TAE)*v_phi(INNER_TH_PART)./bphi(INNER_TH_PART))./(1-(kTAE/omega_TAE)*v_phi(INNER_TH_PART).*bphi(INNER_TH_PART))).*delta_ps(INNER_TH_PART);
        % some particles just can't take this lovely division
        % we give them a moderate kick from approximate expression
        delta_ps(INNER_BLOWUP_PART)=(omega_TAE/nTAE)*delta_ps(INNER_BLOWUP_PART);
%         delta_ps(INNER_BLOWUP_PART)=delta_Epot(INNER_BLOWUP_PART);
        
%         delta_Epot(INNER_PART)=(omega_TAE/nTAE)*delta_ps(INNER_PART);
        if USE_DELTAE_TH==0
            alphas_Etot_th(INNER_PART)=alphas_Etot_th(INNER_PART)+delta_ps(INNER_PART);
        else
        % compromise for precision
            alphas_Etot(INNER_PART)=alphas_Etot(INNER_PART)+delta_ps(INNER_PART);
%             alphas_Etot(INNER_PART)=0.5*(alphas_Etot(INNER_PART)+alphas_Etot_th(INNER_PART));
            delta_E(INNER_PART_GC)=delta_E(INNER_PART_GC)+delta_ps(INNER_PART_GC);
        end
        %         alphas_Ekin_th(INNER_PART_GC)=max(alphas_Etot(INNER_PART_GC)-ZHe*(alphas_Epot_gc(INNER_PART_GC)),0);
    end
    % finally adjust Ekin to make everybody happy
    % compromise between GC evol and instant particle picture
    %     if USE_VD==0
    %         alphas_Epot(INNER_PART)=lininterp3(E_potential_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
    %         if ~isempty(OUTER_PART)
    %             alphas_Epot(OUTER_PART)=0;
    %         end
    alphas_Epot_gc(INNER_PART_GC)=lininterp3(E_potential_map_phi,IL3D_1_GC,IL3D_2_GC,IL3D_3_GC,IL3D_4_GC,IL3D_5_GC,IL3D_6_GC,IL3D_7_GC,IL3D_8_GC, slopex_gc,slopey_gc,slopez_gc);
    if ~isempty(OUTER_PART_GC)
        alphas_Epot_gc(OUTER_PART_GC)=zeros(length(OUTER_PART_GC),1);
    end
    
    alphas_vpll=v_X_step.*bX+v_Z_step.*bZ+v_phi_step.*bphi;
    alphas_vpll_gc=v_X_step.*bX_gc+v_Z_step.*bZ_gc+v_phi_step.*bphi_gc;
    vExB_X=(EZ_gc.*bphi_gc-Ephi_gc.*bZ_gc)./alphas_Bfield_gc;
    vExB_Z=(Ephi_gc.*bX_gc-EX_gc.*bphi_gc)./alphas_Bfield_gc;
    vExB_phi=(EX_gc.*bZ_gc-EZ_gc.*bX_gc)./alphas_Bfield_gc;
    alphas_vE_sq=vExB_X.^2+vExB_Z.^2+vExB_phi.^2;
    alphas_vD_sq=(vD0.^2).*(vDX.^2+vDZ.^2+vDphi.^2);
    alphas_mm_gc=max((alphas_Ekin-0.5*(mHe/eV)*(alphas_vE_sq+alphas_vD_sq+alphas_vpll_gc.^2))./alphas_Bfield_gc,0);
    alphas_Ekin_gc=alphas_mm_gc.*alphas_Bfield_gc+0.5*(mHe/eV)*(alphas_vE_sq+alphas_vD_sq+alphas_vpll_gc.^2);
% 

    alphas_Ekin(INNER_PART_GC)=alphas_Ekin(INNER_PART_GC)-part2wave_vD_power(INNER_PART_GC)*h;
%     if USE_VD==0
%         alphas_Ekin=max(alphas_Etot-ZHe*alphas_Epot,0);
%     else
%         %         alphas_Ekin=alphas_Ekin_gc;
%         alphas_Ekin(INNER_PART_GC)=alphas_Ekin(INNER_PART_GC)-part2wave_vD_power(INNER_PART_GC)*h;
% %         alphas_Ekin=alphas_Ekin_gc+ZHe*(alphas_Epot-alphas_Epot_gc);
%     end
    %         alphas_Ekin=alphas_Ekin-part2wave_vD_power*h;
    %         alphas_Ekin=0.5*(alphas_Ekin)+0.5*max(alphas_Etot-ZHe*alphas_Epot,0);
    %     end
    if (quasi_lin_flag==0)
        %after change of amplitude, do not evolve total energy
        quasi_lin_flag=1;
    end
	
%	gamma_TAE_record=gamma_TAE_record+gamma_TAE;
%     WTAE_record(round(0.5*record_step))=W_TAE;
    Etot_th_record=Etot_th_record+NB_PART_RESCALE*eV*sum(alphas_Etot_th);
    Ekin_record=Ekin_record+NB_PART_RESCALE*eV*sum(alphas_Ekin);
    Etot_record=Etot_record+NB_PART_RESCALE*eV*sum(alphas_Etot);
    pphi_record=pphi_record+alphas_pphi0;
    