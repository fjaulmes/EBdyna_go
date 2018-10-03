% preparing interpolation vectors
%limited map range
alphas_pos_x=min(alphas_pos_x,scale_X(end-3));
alphas_pos_z=min(alphas_pos_z,scale_Z(end-3));
alphas_pos_x=max(alphas_pos_x,scale_X(4));
alphas_pos_z=max(alphas_pos_z,scale_Z(4));
    

Xindex_precise=(alphas_pos_x/DX)+mid_Xzero;
Zindex_precise=(alphas_pos_z/DX)+mid_Z;
%PHIindex=round(wrap(alphas_pos_phi,2)/DPHI);

alphas_Rpos=alphas_pos_x+R0;

Xindex=floor(Xindex_precise);
Zindex=floor(Zindex_precise);
interp_x=Xindex_precise-Xindex;
interp_z=Zindex_precise-Zindex;
INDEX_LIST_1 = sub2ind(size(Btot_XZ_map), Xindex, Zindex );
INDEX_LIST_2 = sub2ind(size(Btot_XZ_map), Xindex, Zindex+1 );
INDEX_LIST_3 = sub2ind(size(Btot_XZ_map), Xindex+1, Zindex );
INDEX_LIST_4 = sub2ind(size(Btot_XZ_map), Xindex+1, Zindex+1 );



alphas_psi_value=interp2_XZ(interp_x,interp_z,psi_XZsmall_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);

alphas_psi=interp1(psi_scale,1:NB_PSI,alphas_psi_value);



INDEX_LIST = sub2ind(size(Btot_XZ_map), round(Xindex_precise), round(Zindex_precise) );
alphas_theta=theta_XZsmall_map(INDEX_LIST);
THETA_BOUNDARY=find((alphas_theta>=HQNB_THETA*DTHETA).*(alphas_theta<=(NB_THETA-HQNB_THETA-2)*DTHETA));
THETA_UP=find(alphas_theta>(NB_THETA-HQNB_THETA-2)*DTHETA);
THETA_LOW=find(alphas_theta<HQNB_THETA*DTHETA);
alphas_theta(THETA_BOUNDARY)=interp2_XZ(interp_x,interp_z,theta_XZsmall_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4,THETA_BOUNDARY);
alphas_theta(THETA_UP)=(1-interp_x(THETA_UP)).*theta_up_XZsmall_map(INDEX_LIST_1(THETA_UP)).*(1-interp_z(THETA_UP))+(interp_x(THETA_UP)).*theta_up_XZsmall_map(INDEX_LIST_3(THETA_UP)).*(1-interp_z(THETA_UP))+(1-interp_x(THETA_UP)).*theta_up_XZsmall_map(INDEX_LIST_2(THETA_UP)).*(interp_z(THETA_UP))+(interp_x(THETA_UP)).*theta_up_XZsmall_map(INDEX_LIST_4(THETA_UP)).*(interp_z(THETA_UP));
alphas_theta(THETA_LOW)=(1-interp_x(THETA_LOW)).*theta_low_XZsmall_map(INDEX_LIST_1(THETA_LOW)).*(1-interp_z(THETA_LOW))+(interp_x(THETA_LOW)).*theta_low_XZsmall_map(INDEX_LIST_3(THETA_LOW)).*(1-interp_z(THETA_LOW))+(1-interp_x(THETA_LOW)).*theta_low_XZsmall_map(INDEX_LIST_2(THETA_LOW)).*(interp_z(THETA_LOW))+(interp_x(THETA_LOW)).*theta_low_XZsmall_map(INDEX_LIST_4(THETA_LOW)).*(interp_z(THETA_LOW));
alphas_theta=wrap2pi(alphas_theta);


