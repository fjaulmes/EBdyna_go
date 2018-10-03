%initialize Bfield properly

Xindex_precise=(alphas_pos_x/DX)+mid_Xzero;
Zindex_precise=(alphas_pos_z/DX)+mid_Z;
Xindex=floor(Xindex_precise);
Zindex=floor(Zindex_precise);
interp_x=Xindex_precise-Xindex;
interp_z=Zindex_precise-Zindex;
INDEX_LIST_1 = sub2ind(size(Btot_XZ_map), Xindex, Zindex );
INDEX_LIST_2 = sub2ind(size(Btot_XZ_map), Xindex, Zindex+1 );
INDEX_LIST_3 = sub2ind(size(Btot_XZ_map), Xindex+1, Zindex );
INDEX_LIST_4 = sub2ind(size(Btot_XZ_map), Xindex+1, Zindex+1 );
INDEX_LIST = sub2ind(size(Btot_XZ_map), round(Xindex_precise), round(Zindex_precise) );
alphas_psi=interp2_XZ(interp_x,interp_z,psi_norm_XZsmall_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
% alphas_psi=psi_XZsmall_map(INDEX_LIST);
% alphas_psi=interp1(psi_scale,1:257,alphas_psi_value);
alphas_theta=theta_XZsmall_map(INDEX_LIST);
THETA_BOUNDARY=find((alphas_theta>=HQNB_THETA*DTHETA).*(alphas_theta<=(NB_THETA-HQNB_THETA-2)*DTHETA));
THETA_UP=find(alphas_theta>(NB_THETA-HQNB_THETA-2)*DTHETA);
THETA_LOW=find(alphas_theta<HQNB_THETA*DTHETA);
alphas_theta(THETA_BOUNDARY)=interp2_XZ(interp_x,interp_z,theta_XZsmall_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4,THETA_BOUNDARY);
alphas_theta(THETA_UP)=(1-interp_x(THETA_UP)).*theta_up_XZsmall_map(INDEX_LIST_1(THETA_UP)).*(1-interp_z(THETA_UP))+(interp_x(THETA_UP)).*theta_up_XZsmall_map(INDEX_LIST_3(THETA_UP)).*(1-interp_z(THETA_UP))+(1-interp_x(THETA_UP)).*theta_up_XZsmall_map(INDEX_LIST_2(THETA_UP)).*(interp_z(THETA_UP))+(interp_x(THETA_UP)).*theta_up_XZsmall_map(INDEX_LIST_4(THETA_UP)).*(interp_z(THETA_UP));
alphas_theta(THETA_LOW)=(1-interp_x(THETA_LOW)).*theta_low_XZsmall_map(INDEX_LIST_1(THETA_LOW)).*(1-interp_z(THETA_LOW))+(interp_x(THETA_LOW)).*theta_low_XZsmall_map(INDEX_LIST_3(THETA_LOW)).*(1-interp_z(THETA_LOW))+(1-interp_x(THETA_LOW)).*theta_low_XZsmall_map(INDEX_LIST_2(THETA_LOW)).*(interp_z(THETA_LOW))+(interp_x(THETA_LOW)).*theta_low_XZsmall_map(INDEX_LIST_4(THETA_LOW)).*(interp_z(THETA_LOW));
alphas_theta=wrap2pi(alphas_theta);

alphas_pos_phi_wrapped=wrap2pi(alphas_pos_phi);
INNER_PART=find(alphas_psi<=(size_r-1));
OUTER_PART=find(alphas_psi>(size_r-1));

% amplitude of the B field at local positions
if ~isempty(OUTER_PART)
    alphas_Bfield(OUTER_PART)=interp2_XZ(interp_x,interp_z,Btot_XZ_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4,OUTER_PART);
end
alphas_Bfield(INNER_PART)=interp3_phi_theta_psi(scale_theta,scale_phi,scale_psi,Btot_map_phi,alphas_theta,alphas_pos_phi_wrapped,alphas_psi,INNER_PART);

% direction of the B field at local positions
if ~isempty(OUTER_PART)
    bX(OUTER_PART)=interp2_XZ(interp_x,interp_z,bX_XZ_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4,OUTER_PART);
    bZ(OUTER_PART)=interp2_XZ(interp_x,interp_z,bZ_XZ_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4,OUTER_PART);
end
bX(INNER_PART)=interp3_phi_theta_psi(scale_theta,scale_phi,scale_psi,bX_map_phi,alphas_theta,alphas_pos_phi_wrapped,alphas_psi,INNER_PART);
bZ(INNER_PART)=interp3_phi_theta_psi(scale_theta,scale_phi,scale_psi,bZ_map_phi,alphas_theta,alphas_pos_phi_wrapped,alphas_psi,INNER_PART);

bphi=sqrt(1-(bX.^2+bZ.^2));