
% amplitude of the B field at local positions
if ~isempty(OUTER_PART)
   alphas_Bfield(OUTER_PART)=interp2_XZ(interp_x,interp_z,Btot_XZ_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4,OUTER_PART);
end
 
alphas_Bfield(INNER_PART)=lininterp3(Btot_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);

% direction of the B field at local positions
if ~isempty(OUTER_PART)
    bX(OUTER_PART)=interp2_XZ(interp_x,interp_z,bX_XZ_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4,OUTER_PART);
    bZ(OUTER_PART)=interp2_XZ(interp_x,interp_z,bZ_XZ_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4,OUTER_PART);
end

bX(INNER_PART)=lininterp3( bX_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
bZ(INNER_PART)=lininterp3( bZ_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
bphi=TOROIDAL_FIELD_DIRECTION*sqrt(1-(bX.^2+bZ.^2));


% % E field at local positions
if ~isempty(OUTER_PART)
    EX(OUTER_PART)=zeros(size(OUTER_PART,1),1);
    EZ(OUTER_PART)=zeros(size(OUTER_PART,1),1);
%     alphas_Epot(OUTER_PART)=zeros(size(OUTER_PART,1),1);
end
EX(INNER_PART)=lininterp3( Efield_X_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
EZ(INNER_PART)=lininterp3( Efield_Z_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);

Ephi=-(EX.*bX+EZ.*bZ)./bphi;

