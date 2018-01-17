    DELTAE_ION=0;
    if (alphas_time_orbit(n)>0)
        alphas_Ekin_t=alphas_Ekin_pop(n)*ones(alphas_time_orbit(n)-TIMEINI+1,1);
        alphas_mm_t=alphas_mm_pop(n)*ones(alphas_time_orbit(n)-TIMEINI+1,1);
        alphas_vD_X=zeros(alphas_time_orbit(n)-TIMEINI+1,1);
        alphas_vD_Z=zeros(alphas_time_orbit(n)-TIMEINI+1,1);
        alphas_vD_phi=zeros(alphas_time_orbit(n)-TIMEINI+1,1);
        alphas_E_X=zeros(alphas_time_orbit(n)-TIMEINI+1,1);
        alphas_E_Z=zeros(alphas_time_orbit(n)-TIMEINI+1,1);
        alphas_E_phi=zeros(alphas_time_orbit(n)-TIMEINI+1,1);
        alphas_psi=zeros(alphas_time_orbit(n)-TIMEINI+1,1);
        alphas_theta=zeros(alphas_time_orbit(n)-TIMEINI+1,1);
        alphas_dB=zeros(alphas_time_orbit(n)-TIMEINI+1,1);
        %for time_stamp=2:alphas_time_orbit(n)
        %time=time_scale(time_stamp);
        alphas_pos_phi_t=squeeze(phipos_output(TIMEINI:alphas_time_orbit(n),PART_POP(n)));
        alphas_pos_phi_wrapped=wrap2pi(alphas_pos_phi_t);
        alphas_pos_x=squeeze(Xpos_gc_output(TIMEINI:alphas_time_orbit(n),PART_POP(n)));
        alphas_pos_z=squeeze(Zpos_gc_output(TIMEINI:alphas_time_orbit(n),PART_POP(n)));
        alphas_vpll_t=squeeze(vpll_output(TIMEINI:alphas_time_orbit(n),PART_POP(n)));
        interpolate_theta_psi_fromXZ;
        %alphas_theta=alphas_theta';
        %alphas_psi=alphas_psi';
        
        %     alphas_psi=min(alphas_psi,size_r);
        alphas_psi_avg=mean(alphas_psi);
        INNER_PART=find(alphas_psi<=(size_r-1));
        OUTER_PART=find(alphas_psi>(size_r-1));
        [IL3D_1 IL3D_2 IL3D_3 IL3D_4 IL3D_5 IL3D_6 IL3D_7 IL3D_8 slopex slopey slopez] = ...
            build_3Dinterp_indexarrays(scale_phi, scale_theta, scale_psi, DPHI,DTHETA,1,alphas_pos_phi_wrapped(INNER_PART), alphas_theta(INNER_PART), alphas_psi(INNER_PART));
        % alphas_vD_X(INNER_PART)=lininterp3(vD_X_PR_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
        % alphas_vD_X(OUTER_PART)=0;
        % alphas_vD_Z(INNER_PART)=lininterp3(vD_Z_PR_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
        % alphas_vD_Z(OUTER_PART)=0;
        % alphas_vD_phi(INNER_PART)=lininterp3(vD_phi_PR_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
        % alphas_vD_phi(OUTER_PART)=0;

        alphas_vD_X(INNER_PART)=interp2(theta_scale,1:size_r,vD_X_map_theta',alphas_theta(INNER_PART),alphas_psi(INNER_PART),'*linear');
        alphas_vD_X(OUTER_PART)=0;
        alphas_vD_Z(INNER_PART)=interp2(theta_scale,1:size_r,vD_Z_map_theta',alphas_theta(INNER_PART),alphas_psi(INNER_PART),'*linear');
        alphas_vD_Z(OUTER_PART)=0;
        alphas_vD_phi(INNER_PART)=interp2(theta_scale,1:size_r,vD_phi_map_theta',alphas_theta(INNER_PART),alphas_psi(INNER_PART),'*linear');
        alphas_vD_phi(OUTER_PART)=0;
        alphas_Eperp=max(alphas_Ekin_t-0.5*(mHe/eV)*alphas_vpll_t.^2,0);
        alphas_vD_coef=(1/ZHe)*(2*alphas_Ekin_t-alphas_Eperp);
        alphas_vD_X=alphas_vD_coef.*alphas_vD_X;
        alphas_vD_Z=alphas_vD_coef.*alphas_vD_Z;
        alphas_vD_phi=alphas_vD_coef.*alphas_vD_phi;

        alphas_E_X(INNER_PART)=interp2(theta_scale,1:size_r,Efield_X_map_theta',alphas_theta(INNER_PART),alphas_psi(INNER_PART),'*linear');
        alphas_E_X(OUTER_PART)=0;
        alphas_E_Z(INNER_PART)=interp2(theta_scale,1:size_r,Efield_Z_map_theta',alphas_theta(INNER_PART),alphas_psi(INNER_PART),'*linear');
        alphas_E_Z(OUTER_PART)=0;
        alphas_E_phi(INNER_PART)=interp2(theta_scale,1:size_r,Efield_phi_map_theta',alphas_theta(INNER_PART),alphas_psi(INNER_PART),'*linear');
        alphas_E_phi(OUTER_PART)=0;

%         alphas_E_X(INNER_PART)=lininterp3(Efield_X_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
%         alphas_E_X(OUTER_PART)=0;
%         alphas_E_Z(INNER_PART)=lininterp3(Efield_Z_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
%         alphas_E_Z(OUTER_PART)=0;
%         alphas_E_phi(INNER_PART)=lininterp3(Efield_phi_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
%         alphas_E_phi(OUTER_PART)=0;
        
         %alphas_dB(INNER_PART)=lininterp3(dBtot_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
         %alphas_dB(OUTER_PART)=0;
        
         alphas_dB(INNER_PART)=interp2(theta_scale,1:size_r,dBtot_map_theta',alphas_theta(INNER_PART),alphas_psi(INNER_PART),'*linear');
         alphas_dB(OUTER_PART)=0;
        
        % the result depends on the scenario for mapped exponential
        % should be rescaled to get gamma=1 and the results will be
        % the ratio of the damping over the growth
        dot_product_vD=alphas_E_X.*alphas_vD_X+alphas_E_Z.*alphas_vD_Z+alphas_E_phi.*alphas_vD_phi;
        POWER_ION_VD=(16/NBSPLITS)*transp_weight*ZHe*eV*mean(dot_product_vD);
        POWER_ION_MU=(16/NBSPLITS)*transp_weight*eV*mean(alphas_mm_t.*alphas_dB);
        POWER_ION=POWER_ION_VD+POWER_ION_MU;
%         POWER_ION=NBSPLITS*transp_weight*eV*mean(alphas_mm_t.*alphas_dB);
        DELTAE_ION=POWER_ION*alphas_delta_time(n);
        %delta_E_ions_evol(time_stamp)=DELTAE_IONS;
        % kinetic energy is so small......
        GAMMA_ION_VD=-0.5*POWER_ION_VD/(WMHD_KINK);
        GAMMA_ION_MU=-0.5*POWER_ION_MU/(WMHD_KINK);
		
		% the mu part should average to 0
        GAMMA_ION=-0.5*POWER_ION_VD/(WMHD_KINK);
        
        alphas_omega=wrap2pi(alphas_theta-alphas_pos_phi_wrapped);
        
        psi_bin_pos=interp1(bins_psi_bounds,1:length(bins_psi_bounds),alphas_psi(INNER_PART))-0.5;
        omega_bin_pos=interp1(bins_omega_bounds,1:length(bins_omega_bounds),alphas_omega(INNER_PART))-0.5;
        theta_bin_pos=interp1(bins_theta_bounds,1:length(bins_theta_bounds),alphas_theta(INNER_PART))-0.5;
        phi_bin_pos=interp1(0:(2*pi/128):2*pi,1:129,alphas_pos_phi_wrapped(INNER_PART))-0.5;
        
        for position=1:length(INNER_PART)
            gamma_vD_phi_psi_map(round(phi_bin_pos(position)),round(psi_bin_pos(position)))=...
                gamma_vD_phi_psi_map(round(phi_bin_pos(position)),round(psi_bin_pos(position)))+GAMMA_ION_VD;
            gamma_vD_omega_psi_map(round(omega_bin_pos(position)),round(psi_bin_pos(position)))=...
                gamma_vD_omega_psi_map(round(omega_bin_pos(position)),round(psi_bin_pos(position)))+GAMMA_ION_VD;
            gamma_vD_theta_psi_map(round(theta_bin_pos(position)),round(psi_bin_pos(position)))=...
                gamma_vD_theta_psi_map(round(theta_bin_pos(position)),round(psi_bin_pos(position)))+GAMMA_ION_VD;
            gamma_mu_omega_psi_map(round(omega_bin_pos(position)),round(psi_bin_pos(position)))=...
                gamma_mu_omega_psi_map(round(omega_bin_pos(position)),round(psi_bin_pos(position)))+GAMMA_ION_MU;
            gamma_mu_theta_psi_map(round(theta_bin_pos(position)),round(psi_bin_pos(position)))=...
                gamma_mu_theta_psi_map(round(theta_bin_pos(position)),round(psi_bin_pos(position)))+GAMMA_ION_MU;
        end
        
    else
        GAMMA_ION=0;
        GAMMA_ION_VD=0;
        GAMMA_ION_MU=0;
    end
    ions_gamma_vD(n)=GAMMA_ION_VD;
    ions_gamma_mu(n)=GAMMA_ION_MU;
    ions_gamma(n)=GAMMA_ION;
    ions_deltaE(n)=DELTAE_ION;
    

    
    if (abs(GAMMA_ION)>120) %&& (abs(GAMMA_ION)>90)
        n
        %
        figure(2)
        hold on
        POS_PART=find((-dot_product_vD>0).*(abs(dot_product_vD)>1e2));
        NEG_PART=find((-dot_product_vD<0).*(abs(dot_product_vD)>1e2));
        
        plot(alphas_theta(POS_PART),alphas_psi(POS_PART),'r.')
        plot(alphas_theta(NEG_PART),alphas_psi(NEG_PART),'b.')
        
        GAMMA_ION
    end
