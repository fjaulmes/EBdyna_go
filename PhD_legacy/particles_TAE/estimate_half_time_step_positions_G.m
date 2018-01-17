    % estimate half time step positions

    alphas_pos_x_int=alphas_pos_x+0.5*h*(v_X_step)+0.125*(h^2)*(ZHe*eV/mHe)*(EX+(v_Z.*bphi_corr-v_phi.*bZ_corr)./alphas_Bfield);
    alphas_pos_z_int=alphas_pos_z+0.5*h*(v_Z_step)+0.125*(h^2)*(ZHe*eV/mHe)*(EZ+(v_phi.*bX_corr-v_X.*bphi_corr)./alphas_Bfield);
    alphas_pos_phi_int=alphas_pos_phi+0.5*h*(v_phi_step)./(0.5*(alphas_pos_x_int+alphas_pos_x)+R0);
    alphas_pos_phi_int=alphas_pos_phi_int+0.125*(h^2)*(ZHe*eV/mHe)*(Ephi+(v_X.*bZ_corr-v_Z.*bX_corr)./alphas_Bfield);

    % limited map range
    alphas_pos_x_int=min(alphas_pos_x_int,scale_X(end-3));
    alphas_pos_z_int=min(alphas_pos_z_int,scale_Z(end-3));
    alphas_pos_x_int=max(alphas_pos_x_int,scale_X(4));
    alphas_pos_z_int=max(alphas_pos_z_int,scale_Z(4));

    alphas_Rpos_int=alphas_pos_x_int+R0;

    % estimating psi and theta values
    Xindex_precise=(alphas_pos_x_int/DX)+mid_Xzero;
    Zindex_precise=(alphas_pos_z_int/DX)+mid_Z;

    Xindex=floor(Xindex_precise);
    Zindex=floor(Zindex_precise);
    interp_x=Xindex_precise-Xindex;
    interp_z=Zindex_precise-Zindex;
    INDEX_LIST_1 = sub2ind(size(Btot_XZ_map), Xindex, Zindex );
    INDEX_LIST_2 = sub2ind(size(Btot_XZ_map), Xindex, Zindex+1 );
    INDEX_LIST_3 = sub2ind(size(Btot_XZ_map), Xindex+1, Zindex );
    INDEX_LIST_4 = sub2ind(size(Btot_XZ_map), Xindex+1, Zindex+1 );

    alphas_pos_phi_wrapped=wrap2pi(alphas_pos_phi_int);
    alphas_pos_phi_TAE=wrapTAEangle(alphas_pos_phi_wrapped,TAE_angle);

	alphas_psi_value=interp2_XZ(interp_x,interp_z,psi_XZsmall_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
	alphas_psiH_value=interp2_XZ(interp_x,interp_z,psiH_XZsmall_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
	if (USE_LAP_PSI==1)
	    alphas_lap_psi=interp2_XZ(interp_x,interp_z,lap_psi_XZsmall_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
	    alphas_psi_value=alphas_psi_value+0.25*(alphas_rhoL.^2).*alphas_lap_psi;
	    alphas_lap_psiH=interp2_XZ(interp_x,interp_z,lap_psiH_XZsmall_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
	    alphas_psiH_value=alphas_psiH_value+0.25*(alphas_rhoL.^2).*alphas_lap_psiH;
    else
%         alphas_psi=interp2_XZ(interp_x,interp_z,psi_norm_XZsmall_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);    % sorting out which particles are inside mixing region and which are out
        alphas_psi=interp1(psi_scale,1:257,alphas_psi_value);
    end

    INNER_PART=find(alphas_psi<=(size_r-1));
    OUTER_PART=find(alphas_psi>(size_r-1));

	INDEX_LIST = sub2ind(size(Btot_XZ_map), round(Xindex_precise), round(Zindex_precise) );
	alphas_theta=theta_XZsmall_map(INDEX_LIST);
	THETA_BOUNDARY=find((alphas_theta>=HQNB_THETA*DTHETA).*(alphas_theta<=(NB_THETA-HQNB_THETA-2)*DTHETA));
	THETA_UP=find(alphas_theta>(NB_THETA-HQNB_THETA-2)*DTHETA);
	THETA_LOW=find(alphas_theta<HQNB_THETA*DTHETA);
	alphas_theta(THETA_BOUNDARY)=interp2_XZ(interp_x,interp_z,theta_XZsmall_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4,THETA_BOUNDARY);
	alphas_theta(THETA_UP)=(1-interp_x(THETA_UP)).*theta_up_XZsmall_map(INDEX_LIST_1(THETA_UP)).*(1-interp_z(THETA_UP))+(interp_x(THETA_UP)).*theta_up_XZsmall_map(INDEX_LIST_3(THETA_UP)).*(1-interp_z(THETA_UP))+(1-interp_x(THETA_UP)).*theta_up_XZsmall_map(INDEX_LIST_2(THETA_UP)).*(interp_z(THETA_UP))+(interp_x(THETA_UP)).*theta_up_XZsmall_map(INDEX_LIST_4(THETA_UP)).*(interp_z(THETA_UP));
	alphas_theta(THETA_LOW)=(1-interp_x(THETA_LOW)).*theta_low_XZsmall_map(INDEX_LIST_1(THETA_LOW)).*(1-interp_z(THETA_LOW))+(interp_x(THETA_LOW)).*theta_low_XZsmall_map(INDEX_LIST_3(THETA_LOW)).*(1-interp_z(THETA_LOW))+(1-interp_x(THETA_LOW)).*theta_low_XZsmall_map(INDEX_LIST_2(THETA_LOW)).*(interp_z(THETA_LOW))+(interp_x(THETA_LOW)).*theta_low_XZsmall_map(INDEX_LIST_4(THETA_LOW)).*(interp_z(THETA_LOW));
	alphas_theta=wrap2pi(alphas_theta);

