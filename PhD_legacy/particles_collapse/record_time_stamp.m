    clear outcast recast
    % correct the particles that are out of simulation domain
    % by giving them a position randomly in the initial distribution
    outcast=find(alphas_psi>Nradial-2);
        
    if (~isempty(outcast))
		recast=randi(Nalphas_simulated,size(outcast,1),1);
		reposition_lost_particles_3DG;
		alphas_ejected(outcast)=1;
    end  
	
    if (mod(time_step,TIME_STAMP_PRECISION)==0)
        adapt_speed_Ekin_G;
%         v_X_step=(1.5*v_X-0.5*v_X_prev);
%         v_Z_step=(1.5*v_Z-0.5*v_Z_prev);
%         v_phi_step=(1.5*v_phi-0.5*v_phi_prev);
        v_X_prev=v_X;
        v_Z_prev=v_Z;
        v_phi_prev=v_phi;
		
        alphas_psi_value=interp2_XZ(interp_x,interp_z,psi_XZsmall_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
        alphas_psiH_value=interp2_XZ(interp_x,interp_z,psiH_XZsmall_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
       
        time_stamp=ceil((time_step)/TIME_STAMP_PRECISION);
        alphas_omega=wrap2pi(alphas_theta-alphas_pos_phi_wrapped);
        %alphas_psi_star(INNER_PART)=interp2(1:size_r,scale_theta,psi_star_omega_map',alphas_psi(INNER_PART), alphas_omega(INNER_PART),'*linear');
        if ~isempty(OUTER_PART)
            alphas_psi_value_corr(OUTER_PART)=alphas_psi_value(OUTER_PART);
        end
        alphas_psi_value_corr(INNER_PART)=alphas_psiH_value(INNER_PART)+alphas_psi_star(INNER_PART);
        
        % rescaling the value of pphi0 to our knowledge of the toroidal
        % speed ans psi value
        if CALCULATE_TRUE_PPHI==0
			alphas_pphi0=(mHe/eV)*alphas_Rpos.*v_phi-(ZHe)*alphas_psi_value_corr;
        end
		
        % amplitude of the B field and potential at time step positions
        update_Gfields_collapse;

        alphas_vpll=v_X_step.*bX+v_Z_step.*bZ+v_phi_step.*bphi;
        vExB_X=(EZ.*bphi-Ephi.*bZ)./alphas_Bfield;
        vExB_Z=(Ephi.*bX-EX.*bphi)./alphas_Bfield;
        vExB_phi=(EX.*bZ-EX.*bX)./alphas_Bfield;
        alphas_vE_sq=vExB_X.^2+vExB_Z.^2+vExB_phi.^2;
        alphas_Epll=0.5*(mHe/eV)*alphas_vpll.^2;
        alphas_Eperp=max(0.5*(mHe/eV)*(v_X_step.^2+v_Z_step.^2+v_phi_step.^2)-alphas_Epll,0);
		alphas_mm=alphas_Eperp./alphas_Bfield;
        outcast=find(alphas_psi>simulation_size_r+27);
        alphas_vperp=sqrt(2*alphas_Eperp*eV/mHe);
        alphas_rhoL=(mHe/eV)*(1/ZHe)*alphas_vperp./alphas_Bfield;
        dr_X=interp2(scale_theta,1:Nradial,dr_X_PR_map', alphas_omega,alphas_psi,'*linear');
        dr_Z=interp2(scale_theta,1:Nradial,dr_Z_PR_map', alphas_omega,alphas_psi,'*linear');
		alphas_vEradial=vExB_X.*dr_X+vExB_Z.*dr_Z;
        if (~isempty(outcast))
            recast=randi(Nalphas_simulated,size(outcast,1),1);
            reposition_lost_particles_3DG;
            alphas_ejected(outcast)=1;
            disp(strcat('number of ejected particles =  ',num2str(size(outcast,1))));
        end
    end
