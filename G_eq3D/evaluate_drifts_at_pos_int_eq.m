%     alphas_psi_value=interp2_XZ(interp_x,interp_z,psi_XZsmall_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
    alphas_psi_value=max(alphas_psi_value,-psi_global);
    alphas_psi=interp1(psi_scale,1:257,alphas_psi_value);

    INDEX_LIST = sub2ind(size(Btot_XZ_map), round(Xindex_precise), round(Zindex_precise) );
    alphas_theta=theta_XZsmall_map(INDEX_LIST);
    THETA_BOUNDARY=find((alphas_theta>=HQNB_THETA*DTHETA).*(alphas_theta<=(NB_THETA-HQNB_THETA-2)*DTHETA));
    THETA_UP=find(alphas_theta>(NB_THETA-HQNB_THETA-2)*DTHETA);
    THETA_LOW=find(alphas_theta<HQNB_THETA*DTHETA);
    alphas_theta(THETA_BOUNDARY)=interp2_XZ(interp_x,interp_z,theta_XZsmall_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4,THETA_BOUNDARY);
    alphas_theta(THETA_UP)=(1-interp_x(THETA_UP)).*theta_up_XZsmall_map(INDEX_LIST_1(THETA_UP)).*(1-interp_z(THETA_UP))+(interp_x(THETA_UP)).*theta_up_XZsmall_map(INDEX_LIST_3(THETA_UP)).*(1-interp_z(THETA_UP))+(1-interp_x(THETA_UP)).*theta_up_XZsmall_map(INDEX_LIST_2(THETA_UP)).*(interp_z(THETA_UP))+(interp_x(THETA_UP)).*theta_up_XZsmall_map(INDEX_LIST_4(THETA_UP)).*(interp_z(THETA_UP));
    alphas_theta(THETA_LOW)=(1-interp_x(THETA_LOW)).*theta_low_XZsmall_map(INDEX_LIST_1(THETA_LOW)).*(1-interp_z(THETA_LOW))+(interp_x(THETA_LOW)).*theta_low_XZsmall_map(INDEX_LIST_3(THETA_LOW)).*(1-interp_z(THETA_LOW))+(1-interp_x(THETA_LOW)).*theta_low_XZsmall_map(INDEX_LIST_2(THETA_LOW)).*(interp_z(THETA_LOW))+(interp_x(THETA_LOW)).*theta_low_XZsmall_map(INDEX_LIST_4(THETA_LOW)).*(interp_z(THETA_LOW));
    alphas_theta=wrap2pi(alphas_theta);


    % amplitude of the B field and potential at half time step local positions
    alphas_Bfield=interp2_XZ(interp_x,interp_z,Btot_XZ_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
    
    % direction of the B field at local positions
    bX=interp2_XZ(interp_x,interp_z,bX_XZ_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
    bZ=interp2_XZ(interp_x,interp_z,bZ_XZ_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
    bphi=sqrt(1-(bX.^2+bZ.^2));

    alphas_Eperp=alphas_mm_part.*alphas_Bfield;


    vD_coef=(2*alphas_Ekin-alphas_Eperp)/(ZHe);


    vD_X=vD_coef.*vDraw_X;
    vD_Z=vD_coef.*vDraw_Z;
%     vD_phi=vD_coef.*vDraw_phi;
    vD_phi=-(vD_X.*bX+vD_Z.*bZ)./bphi;

    alphas_vD_sq=(vD_X.^2+vD_Z.^2+vD_phi.^2);
