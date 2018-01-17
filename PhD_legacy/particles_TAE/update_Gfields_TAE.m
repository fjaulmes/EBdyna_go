
% amplitude of the B field at local positions
% if NO_PERTURB_ORBITS==0
%    if ~isempty(OUTER_PART)
%        alphas_Bfield(OUTER_PART)=interp2_XZ(interp_x,interp_z,Btot_XZ_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4,OUTER_PART);
%    end
%  
%    alphas_Bfield(INNER_PART)=lininterp3(Btot_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
% else
%    alphas_Bfield=interp2_XZ(interp_x,interp_z,Btot_XZ_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
% end
% alphas_Bfield=interp2_XZ(interp_x,interp_z,Btot_XZ_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);

% direction of the B field at local positions
% if ~isempty(OUTER_PART)
    bX=interp2_XZ(interp_x,interp_z,BpolX_initial_XZsmall_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
    bZ=interp2_XZ(interp_x,interp_z,BpolZ_initial_XZsmall_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
	bphi=interp2_XZ(interp_x,interp_z,Bphi_XZsmall_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
% end

%adding TAE perturbation on top of equilibrium value
if NO_PERTURB_ORBITS==0
   bX(INNER_PART)=bX(INNER_PART)+lininterp3( BsX_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
   bZ(INNER_PART)=bZ(INNER_PART)+lininterp3( BsZ_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
   bphi(INNER_PART)=bphi(INNER_PART)+lininterp3( Bsphi_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
end
alphas_Bfield=sqrt(bX.^2+bZ.^2+bphi.^2);
bX=bX./alphas_Bfield;
bZ=bZ./alphas_Bfield;
bphi=bphi./alphas_Bfield;
% bphi=TOROIDAL_FIELD_DIRECTION*sqrt(1-(bX.^2+bZ.^2));


% % E field at local positions
if NO_PERTURB_ORBITS==0
    if ~isempty(OUTER_PART)
        EX(OUTER_PART)=zeros(size(OUTER_PART,1),1);
        EZ(OUTER_PART)=zeros(size(OUTER_PART,1),1);
		Ephi(OUTER_PART)=zeros(size(OUTER_PART,1),1);
%     alphas_Epot(OUTER_PART)=zeros(size(OUTER_PART,1),1);
    end
    EX(INNER_PART)=lininterp3( Efield_X_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
    EZ(INNER_PART)=lininterp3( Efield_Z_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
	
	Ephi(INNER_PART)=omega_TAE*bphi(INNER_PART).*alphas_iA(INNER_PART)-nTAE*alphas_iEpot(INNER_PART)./alphas_Rpos(INNER_PART);

%     Ephi=-(EX.*bX+EZ.*bZ)./bphi;
	
end


