
        TAE_amplitude
        record_step=0;
        calculate_gc_pos;
		calculate_vD_E_gc_components;

       
        
		alphas_psi_value=interp2_XZ(interp_x,interp_z,psi_XZsmall_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
        
        time_stamp=ceil((time_step)/TIME_STAMP_PRECISION);


		alphas_Epot_gc=alphas_Epot;
 
        %adapt_speed_Ekin_G;
 
        %alphas_Epot(INNER_PART)=lininterp3(E_potential_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
        %if ~isempty(OUTER_PART)
        %    alphas_Epot(OUTER_PART)=0;
        %end
        alphas_psi_star=-bphi.*alphas_Apll_value.*alphas_Rpos;
        
        if ~isempty(OUTER_PART)
            alphas_psi_value_corr(OUTER_PART)=alphas_psi_value(OUTER_PART);
        end
        alphas_psi_value_corr(INNER_PART)=alphas_psi_value(INNER_PART)+alphas_psi_star(INNER_PART);
        
        % rescaling the value of pphi0 to our knowledge of the toroidal
        % speed ans psi value
        if (CALCULATE_TRUE_PPHI==0)||(f_counter<=3)
            alphas_pphi0=(mHe/eV)*(alphas_pos_x+R0).*v_phi_step-(ZHe)*alphas_psi_value_corr;
            alphas_Etot=alphas_Ekin+ZHe*alphas_Epot_gc;
            alphas_Etot_th=alphas_Etot;
			Epart_tot_rel_ini=0;
            Epart_tot_vD_rel_ini=0;
        end
        % amplitude of the B field and potential at half time step local positions
%         update_Gfields_collapse;

        alphas_vpll=v_X_step.*bX+v_Z_step.*bZ+v_phi_step.*bphi;
        alphas_vpll_gc=v_X_step.*bX_gc+v_Z_step.*bZ_gc+v_phi_step.*bphi_gc;
        vExB_X=(EZ_gc.*bphi_gc-Ephi_gc.*bZ_gc)./alphas_Bfield_gc;
        vExB_Z=(Ephi_gc.*bX_gc-EX_gc.*bphi_gc)./alphas_Bfield_gc;
        vExB_phi=(EX_gc.*bZ_gc-EZ_gc.*bX_gc)./alphas_Bfield_gc;
        alphas_vE_sq=vExB_X.^2+vExB_Z.^2+vExB_phi.^2;
        alphas_vD_sq=(vD0.^2).*(vDX.^2+vDZ.^2+vDphi.^2);
        alphas_Epll=0.5*(mHe/eV)*alphas_vpll.^2;
        alphas_Eperp=max(0.5*(mHe/eV)*(v_X_step.^2+v_Z_step.^2+v_phi_step.^2)-alphas_Epll,0);
        alphas_mm_part=alphas_Eperp./alphas_Bfield;
        alphas_mm_gc=max((alphas_Ekin-0.5*(mHe/eV)*(alphas_vE_sq+alphas_vD_sq+alphas_vpll_gc.^2))./alphas_Bfield_gc,0);
        %alphas_Ekin_gc=alphas_mm_gc.*alphas_Bfield_gc+0.5*(mHe/eV)*(alphas_vE_sq+alphas_vD_sq+alphas_vpll_gc.^2);
        
        alphas_vperp=sqrt(2*alphas_Eperp*eV/mHe);
        alphas_rhoL=(mHe/eV)*(1/ZHe)*alphas_vperp./alphas_Bfield;
%         dr_X=interp2(scale_theta,1:Nradial,dr_X_PR_map', alphas_omega,alphas_psi,'*linear');
%         dr_Z=interp2(scale_theta,1:Nradial,dr_Z_PR_map', alphas_omega,alphas_psi,'*linear');
        outcast=find(alphas_psi>NB_PSI-2);
        if (~isempty(outcast))
            recast=randi(Nalphas_simulated,size(outcast,1),1);
            reposition_lost_particles_3DG;
            alphas_ejected(outcast)=1;
            disp(strcat('number of ejected particles =  ',num2str(size(outcast,1))));
        end
% 		 pos_X_gc=(mHe/eV)*(1/ZHe)*(v_Z_step.*bphi-v_phi_step.*bZ)./alphas_Bfield;
%        pos_Z_gc=(mHe/eV)*(1/ZHe)*(v_phi_step.*bX-v_X_step.*bphi)./alphas_Bfield;

        alphas_Ekin_gc=0.5*(mHe/eV)*(v_X_step.^2+v_Z_step.^2+v_phi_step.^2);
        
        Xpos_part_output(time_stamp,:)=alphas_pos_x;
        Xpos_output(time_stamp,:)=alphas_pos_x+pos_X_gc;
        Zpos_output(time_stamp,:)=alphas_pos_z+pos_Z_gc;
        phipos_output(time_stamp,:)=alphas_pos_phi_gc;
        psipos_output(time_stamp,:)=alphas_psi;
        vparallel_output(time_stamp,:)=alphas_vpll_gc;
        Eperp_output(time_stamp,:)=alphas_Eperp;
        vD_output(time_stamp,:)=sqrt(alphas_vD_sq);
        vEsq_output(time_stamp,:)=alphas_vE_sq;
        theta_output(time_stamp,:)=alphas_theta;
        Ekin_output(time_stamp,:)=alphas_Ekin;
        Ekin_gc_output(time_stamp,:)=alphas_Ekin_gc;
        vphi_output(time_stamp,:)=v_phi_step;
        pphi_output(time_stamp,:)=alphas_pphi0;
        psi_star_output(time_stamp,:)=alphas_psi_star;
        psi_value_output(time_stamp,:)=alphas_psi_value_corr;
        rhoL_output(time_stamp,:)=alphas_rhoL;
        Epot_output(time_stamp,:)=alphas_Epot_gc;
        Etot_output(time_stamp,:)=alphas_Etot;
        Bfield_output(time_stamp,:)=alphas_Bfield;
        mm_output(time_stamp,:)=alphas_mm_gc;
%         vEradial_output(time_stamp,:)=vExB_X.*dr_X+vExB_Z.*dr_Z;
        if CALCULATE_VD_POWER==1
            part2wave_power_output(time_stamp,:)=part2wave_power;
            part2wave_vD_power_output(time_stamp,:)=part2wave_vD_power;
            part2wave_vD_power_evol(time_stamp)=P_VD_PART_TAE;
        end
        power_exchange_evol(time_stamp,:)=Pexchange_part_record*TIME_GO_SIZE/TIME_STAMP_SIZE;
        pphi_evol(time_stamp,:)=pphi_record*TIME_GO_SIZE/TIME_STAMP_SIZE;
        Pexchange_part_record=zeros(Nalphas_simulated,1);
        pphi_record=zeros(Nalphas_simulated,1);