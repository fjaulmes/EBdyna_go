
        
		if (mod(time_step,RECORD_PRECISION)==0)
		    disp('verifying kinetic energy consistency:')
			max(alphas_Ekin)
			
			v_X_half=v_X;
			v_Z_half=v_Z;
			v_phi_half=v_phi;
			
			update_GT_3D_eq_centrifug;
			
			v_X_step=(0.5*v_X+0.5*v_X_prev);
			v_Z_step=(0.5*v_Z+0.5*v_Z_prev);
			v_phi_step=(0.5*v_phi+0.5*v_phi_prev);
	 
			v_X=v_X_half;
			v_Z=v_Z_half;
			v_phi=v_phi_half;
	 
			v_X_prev=v_X_half;
			v_Z_prev=v_Z_half;
			v_phi_prev=v_phi_half;
			
			time_stamp=ceil((time_step)/RECORD_PRECISION);
			alphas_omega=wrap2pi(alphas_theta-alphas_pos_phi_wrapped);
			alphas_psi_value_corr=alphas_psi_value;
			
			% rescaling the value of pphi0 to our knowledge of the toroidal
			% speed ans psi value
			% alphas_pphi0=(mHe/eV)*alphas_Rpos.*v_phi-(ZHe)*alphas_psi_value_corr;

			% amplitude of the B field and potential at time step local positions

			alphas_Bfield=interp2_XZ(interp_x,interp_z,Btot_XZ_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);

			% direction of the B field at local positions

			bX=interp2_XZ(interp_x,interp_z,bX_XZ_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
			bZ=interp2_XZ(interp_x,interp_z,bZ_XZ_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);

			bphi=TOROIDAL_FIELD_DIRECTION*sqrt(1-(bX.^2+bZ.^2));

			alphas_vpll=v_X_step.*bX+v_Z_step.*bZ+v_phi_step.*bphi;
			alphas_Epll=0.5*(mHe/eV)*alphas_vpll.^2;
			alphas_Eperp=max(alphas_Ekin-alphas_Epll,0);
			outcast=find(alphas_psi>SIMULATION_RADIAL_LIMIT);
			alphas_vperp=sqrt(2*alphas_Eperp*eV/mHe);
			alphas_rhoL=(mHe/eV)*(1/ZHe)*alphas_vperp./alphas_Bfield;
			if (~isempty(outcast))
				recast=randi(Nalphas_simulated,size(outcast,1),1);
				reposition_lost_particles_3DG;
				alphas_ejected(outcast)=1;
				disp(strcat('number of ejected particles =  ',num2str(size(outcast,1))));
			end
			
			pos_X_gc=(mHe/eV)*(1/ZHe)*(v_Z_step.*bphi-v_phi_step.*bZ)./alphas_Bfield;
			pos_Z_gc=(mHe/eV)*(1/ZHe)*(v_phi_step.*bX-v_X_step.*bphi)./alphas_Bfield;
	%         pos_phi_gc=(mHe/eV)*(1/ZHe)*(v_X.*bZ-v_Z.*bX)./alphas_Bfield;

			
			% recording for precession calculations
	%         phipos_output(time_stamp,:)=alphas_pos_phi;
	%         theta_output(time_stamp,:)=alphas_theta;
			vpll_output(time_stamp,:)=alphas_vpll;
			Xpos_gc_output(time_stamp,:)=alphas_pos_x+pos_X_gc;
			Zpos_gc_output(time_stamp,:)=alphas_pos_z+pos_Z_gc;
			phipos_output(time_stamp,:)=alphas_pos_phi;
        end