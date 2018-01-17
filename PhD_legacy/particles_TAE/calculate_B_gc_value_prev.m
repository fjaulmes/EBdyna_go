%evaluating the time step value
v_X_step=(1.5*v_X-0.5*v_X_prev);
v_Z_step=(1.5*v_Z-0.5*v_Z_prev);
v_phi_step=(1.5*v_phi-0.5*v_phi_prev);


pos_X_gc=(mHe/eV)*(1/ZHe)*(v_Z_step.*bphi-v_phi_step.*bZ)./alphas_Bfield;
pos_Z_gc=(mHe/eV)*(1/ZHe)*(v_phi_step.*bX-v_X_step.*bphi)./alphas_Bfield;
    
Xindex_gc_precise=((alphas_pos_x+pos_X_gc)/DX)+mid_Xzero;
Zindex_gc_precise=((alphas_pos_z+pos_Z_gc)/DX)+mid_Z;

alphas_Rpos_gc=alphas_pos_x+pos_X_gc+R0;
% alphas_pos_phi_gc=alphas_pos_phi.*(alphas_Rpos./alphas_Rpos_gc)+(mHe/eV)*(1/ZHe)*(v_X_step.*bZ-v_Z_step.*bX)./(alphas_Rpos_gc.*alphas_Bfield);
alphas_pos_phi_gc=alphas_pos_phi+(mHe/eV)*(1/ZHe)*(v_X_step.*bZ-v_Z_step.*bX)./(alphas_Rpos_gc.*alphas_Bfield);
alphas_pos_phi_gc_wrapped=wrap2pi(alphas_pos_phi_gc);
alphas_pos_phi_gc_TAE=wrapTAEangle(alphas_pos_phi_gc_wrapped,TAE_ANGLE);


Xindex=floor(Xindex_gc_precise);
Zindex=floor(Zindex_gc_precise);
interp_x_gc=Xindex_gc_precise-Xindex;
interp_z_gc=Zindex_gc_precise-Zindex;
INDEX_LIST_1_GC = sub2ind(size(Btot_XZ_map), Xindex, Zindex );
INDEX_LIST_2_GC = sub2ind(size(Btot_XZ_map), Xindex, Zindex+1 );
INDEX_LIST_3_GC = sub2ind(size(Btot_XZ_map), Xindex+1, Zindex );
INDEX_LIST_4_GC = sub2ind(size(Btot_XZ_map), Xindex+1, Zindex+1 );


alphas_psi=interp2_XZ(interp_x_gc,interp_z_gc,psi_norm_XZsmall_map,INDEX_LIST_1_GC,INDEX_LIST_2_GC,INDEX_LIST_3_GC,INDEX_LIST_4_GC);


INDEX_LIST = sub2ind(size(Btot_XZ_map), round(Xindex_gc_precise), round(Zindex_gc_precise) );
alphas_theta=theta_XZsmall_map(INDEX_LIST);
THETA_BOUNDARY=find((alphas_theta>=HQNB_THETA*DTHETA).*(alphas_theta<=(NB_THETA-HQNB_THETA-2)*DTHETA));
THETA_UP=find(alphas_theta>(NB_THETA-HQNB_THETA-2)*DTHETA);
THETA_LOW=find(alphas_theta<HQNB_THETA*DTHETA);
alphas_theta(THETA_BOUNDARY)=interp2_XZ(interp_x_gc,interp_z_gc,theta_XZsmall_map,INDEX_LIST_1_GC,INDEX_LIST_2_GC,INDEX_LIST_3_GC,INDEX_LIST_4_GC,THETA_BOUNDARY);
alphas_theta(THETA_UP)=(1-interp_x_gc(THETA_UP)).*theta_up_XZsmall_map(INDEX_LIST_1_GC(THETA_UP)).*(1-interp_z_gc(THETA_UP))+(interp_x_gc(THETA_UP)).*theta_up_XZsmall_map(INDEX_LIST_3_GC(THETA_UP)).*(1-interp_z_gc(THETA_UP))+(1-interp_x_gc(THETA_UP)).*theta_up_XZsmall_map(INDEX_LIST_2_GC(THETA_UP)).*(interp_z_gc(THETA_UP))+(interp_x_gc(THETA_UP)).*theta_up_XZsmall_map(INDEX_LIST_4_GC(THETA_UP)).*(interp_z_gc(THETA_UP));
alphas_theta(THETA_LOW)=(1-interp_x_gc(THETA_LOW)).*theta_low_XZsmall_map(INDEX_LIST_1_GC(THETA_LOW)).*(1-interp_z_gc(THETA_LOW))+(interp_x_gc(THETA_LOW)).*theta_low_XZsmall_map(INDEX_LIST_3_GC(THETA_LOW)).*(1-interp_z_gc(THETA_LOW))+(1-interp_x_gc(THETA_LOW)).*theta_low_XZsmall_map(INDEX_LIST_2_GC(THETA_LOW)).*(interp_z_gc(THETA_LOW))+(interp_x_gc(THETA_LOW)).*theta_low_XZsmall_map(INDEX_LIST_4_GC(THETA_LOW)).*(interp_z_gc(THETA_LOW));
alphas_theta=wrap2pi(alphas_theta);

% sorting out which particles are inside mixing region and which are out
INNER_PART_GC=find((alphas_psi<=(pTAE_sup-1)).*(alphas_psi>=(pTAE_inf+1)));
OUTER_PART_GC=find((alphas_psi>(pTAE_sup-1))+(alphas_psi<(pTAE_inf+1))>0);
[IL3D_1_GC IL3D_2_GC IL3D_3_GC IL3D_4_GC IL3D_5_GC IL3D_6_GC IL3D_7_GC IL3D_8_GC slopex_gc slopey_gc slopez_gc] = ...
    build_3Dinterp_indexarrays(scale_phi, scale_theta, scale_psi, DPHI,DTHETA,1,alphas_pos_phi_gc_TAE(INNER_PART_GC), alphas_theta(INNER_PART_GC), alphas_psi(INNER_PART_GC));

if ~isempty(OUTER_PART_GC)
   alphas_Bfield_gc_prev(OUTER_PART_GC)=interp2_XZ(interp_x_gc,interp_z_gc,Btot_XZ_map,INDEX_LIST_1_GC,INDEX_LIST_2_GC,INDEX_LIST_3_GC,INDEX_LIST_4_GC,OUTER_PART_GC);
end
alphas_Bfield_gc_prev(INNER_PART_GC)=lininterp3(Btot_map_phi_prev,IL3D_1_GC,IL3D_2_GC,IL3D_3_GC,IL3D_4_GC,IL3D_5_GC,IL3D_6_GC,IL3D_7_GC,IL3D_8_GC, slopex_gc,slopey_gc,slopez_gc);


