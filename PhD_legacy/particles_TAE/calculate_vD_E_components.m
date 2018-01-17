alphas_vpll=v_X_step.*bX+v_Z_step.*bZ+v_phi_step.*bphi;
alphas_Epll=0.5*(mHe/eV)*alphas_vpll.^2;
alphas_Eperp=max(0.5*(mHe/eV)*(v_X_step.^2+v_Z_step.^2+v_phi_step.^2)-alphas_Epll,0);


% pos_X_gc=(mHe/eV)*(1/ZHe)*(v_Z_step.*bphi-v_phi_step.*bZ)./alphas_Bfield;
% pos_Z_gc=(mHe/eV)*(1/ZHe)*(v_phi_step.*bX-v_X_step.*bphi)./alphas_Bfield;
%     
% Xindex_precise=((alphas_pos_x+pos_X_gc)/DX)+mid_Xzero;
% Zindex_precise=((alphas_pos_z+pos_Z_gc)/DX)+mid_Z;
% 
% alphas_pos_phi_gc=alphas_pos_phi+(mHe/eV)*(1/ZHe)*(v_X_step.*bZ-v_Z_step.*bX)./(alphas_Rpos.*alphas_Bfield);
% alphas_pos_phi_gc_wrapped=wrap2pi(alphas_pos_phi);
% alphas_pos_phi_gc_TAE=wrapTAEangle(alphas_pos_phi_gc_wrapped,TAE_ANGLE);
% 
% 
% Xindex=floor(Xindex_precise);
% Zindex=floor(Zindex_precise);
% interp_x_gc=Xindex_precise-Xindex;
% interp_z_gc=Zindex_precise-Zindex;
% INDEX_LIST_1 = sub2ind(size(Btot_XZ_map), Xindex, Zindex );
% INDEX_LIST_2 = sub2ind(size(Btot_XZ_map), Xindex, Zindex+1 );
% INDEX_LIST_3 = sub2ind(size(Btot_XZ_map), Xindex+1, Zindex );
% INDEX_LIST_4 = sub2ind(size(Btot_XZ_map), Xindex+1, Zindex+1 );
% 
% 
% alphas_psi=interp2_XZ(interp_x_gc,interp_z_gc,psi_norm_XZsmall_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
% alphas_Bfield_gc=alphas_Bfield;
% 
% 
% 
% INDEX_LIST = sub2ind(size(Btot_XZ_map), round(Xindex_precise), round(Zindex_precise) );
% alphas_theta=theta_XZsmall_map(INDEX_LIST);
% THETA_BOUNDARY=find((alphas_theta>=HQNB_THETA*DTHETA).*(alphas_theta<=(NB_THETA-HQNB_THETA-2)*DTHETA));
% THETA_UP=find(alphas_theta>(NB_THETA-HQNB_THETA-2)*DTHETA);
% THETA_LOW=find(alphas_theta<HQNB_THETA*DTHETA);
% alphas_theta(THETA_BOUNDARY)=interp2_XZ(interp_x_gc,interp_z_gc,theta_XZsmall_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4,THETA_BOUNDARY);
% alphas_theta(THETA_UP)=(1-interp_x_gc(THETA_UP)).*theta_up_XZsmall_map(INDEX_LIST_1(THETA_UP)).*(1-interp_z_gc(THETA_UP))+(interp_x_gc(THETA_UP)).*theta_up_XZsmall_map(INDEX_LIST_3(THETA_UP)).*(1-interp_z_gc(THETA_UP))+(1-interp_x_gc(THETA_UP)).*theta_up_XZsmall_map(INDEX_LIST_2(THETA_UP)).*(interp_z_gc(THETA_UP))+(interp_x_gc(THETA_UP)).*theta_up_XZsmall_map(INDEX_LIST_4(THETA_UP)).*(interp_z_gc(THETA_UP));
% alphas_theta(THETA_LOW)=(1-interp_x_gc(THETA_LOW)).*theta_low_XZsmall_map(INDEX_LIST_1(THETA_LOW)).*(1-interp_z_gc(THETA_LOW))+(interp_x_gc(THETA_LOW)).*theta_low_XZsmall_map(INDEX_LIST_3(THETA_LOW)).*(1-interp_z_gc(THETA_LOW))+(1-interp_x_gc(THETA_LOW)).*theta_low_XZsmall_map(INDEX_LIST_2(THETA_LOW)).*(interp_z_gc(THETA_LOW))+(interp_x_gc(THETA_LOW)).*theta_low_XZsmall_map(INDEX_LIST_4(THETA_LOW)).*(interp_z_gc(THETA_LOW));
% alphas_theta=wrap2pi(alphas_theta);
% 
% % sorting out which particles are inside mixing region and which are out
% INNER_PART=find((alphas_psi<=(pTAE_sup-1)).*(alphas_psi>=(pTAE_inf+1)));
% OUTER_PART=find((alphas_psi>(pTAE_sup-1))+(alphas_psi<(pTAE_inf+1))>0);
% [IL3D_1 IL3D_2 IL3D_3 IL3D_4 IL3D_5 IL3D_6 IL3D_7 IL3D_8 slopex_gc slopey_gc slopez_gc] = ...
%     build_3Dinterp_indexarrays(scale_phi, scale_theta, scale_psi, DPHI,DTHETA,1,alphas_pos_phi_gc_TAE(INNER_PART), alphas_theta(INNER_PART), alphas_psi(INNER_PART));
% 
if ~isempty(OUTER_PART_GC)
   alphas_Bfield_gc(OUTER_PART_GC)=interp2_XZ(interp_x_gc,interp_z_gc,Btot_XZ_map,INDEX_LIST_1_GC,INDEX_LIST_2_GC,INDEX_LIST_3_GC,INDEX_LIST_4_GC,OUTER_PART_GC);
end
alphas_Bfield_gc(INNER_PART_GC)=lininterp3(Btot_map_phi,IL3D_1_GC,IL3D_2_GC,IL3D_3_GC,IL3D_4_GC,IL3D_5_GC,IL3D_6_GC,IL3D_7_GC,IL3D_8_GC, slopex_gc,slopey_gc,slopez_gc);


vD0=(((2*alphas_Ekin)-alphas_Eperp)/ZHe)./alphas_Bfield.^3;
%vDX(INNER_PART)=interp2_XZ(interp_x_gc,interp_z_gc,vD_X_XZ_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4,INNER_PART);
%vDZ(INNER_PART)=interp2_XZ(interp_x_gc,interp_z_gc,vD_Z_XZ_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4,INNER_PART);
%vDphi(INNER_PART)=interp2_XZ(interp_x_gc,interp_z_gc,vD_phi_XZ_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4,INNER_PART);

% now we improve the calculation to get vD from equilibrium gradient
gB_X=interp2_XZ(interp_x_gc,interp_z_gc,gB0_X_XZ_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
gB_Z=interp2_XZ(interp_x_gc,interp_z_gc,gB0_Z_XZ_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);



bX_gc=interp2_XZ(interp_x_gc,interp_z_gc,BpolX_initial_XZsmall_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
bZ_gc=interp2_XZ(interp_x_gc,interp_z_gc,BpolZ_initial_XZsmall_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
BX_gc=bX_gc;
BZ_gc=bZ_gc;
% bphi_gc=interp2_XZ(interp_x_gc,interp_z_gc,Bphi_XZsmall_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
% bphi_gc(OUTER_PART)=bphi_gc(OUTER_PART)./alphas_Bfield_gc(OUTER_PART);


EX_gc=zeros(Nalphas_simulated,1);
EZ_gc=zeros(Nalphas_simulated,1);
Ephi_gc=zeros(Nalphas_simulated,1);

EX_gc(INNER_PART)=lininterp3( Efield_X_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
EZ_gc(INNER_PART)=lininterp3( Efield_Z_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
BX_gc(INNER_PART)=BX_gc(INNER_PART)+lininterp3( BsX_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
BZ_gc(INNER_PART)=BZ_gc(INNER_PART)+lininterp3( BsZ_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
% Btot(INNER_PART)=interp2_XZ(interp_x_gc,interp_z_gc,Btot_map_phi,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4,INNER_PART);
bX_gc=BX_gc./alphas_Bfield;
bZ_gc=BZ_gc./alphas_Bfield;
bphi_gc=sqrt(1-(bX_gc.^2+bZ_gc.^2));

Bphi_gc=alphas_Bfield.*bphi_gc;

vDX=-Bphi_gc.*gB_Z;
vDZ=Bphi_gc.*gB_X;
vDphi=BX_gc.*gB_Z-BZ_gc.*gB_X;



Ephi_gc(INNER_PART)=-(EX_gc(INNER_PART).*bX_gc(INNER_PART)+EZ_gc(INNER_PART).*bZ_gc(INNER_PART))./(bphi_gc(INNER_PART));


