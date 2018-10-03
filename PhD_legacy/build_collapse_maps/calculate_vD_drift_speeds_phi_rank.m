

calculate_gradB_phi_rank;

vD_X_map=-Bphi_XZsmall_map.*gradB_Z;
vD_Z_map=Bphi_XZsmall_map.*gradB_X;
vD_phi_map=BpolX_XZ_map.*gradB_Z-BpolZ_XZ_map.*gradB_X;
%         vD_X_map=vD_coef*vD_X_map./(Btot_XZ_map.^3);
%         vD_Z_map=vD_coef*vD_Z_map./(Btot_XZ_map.^3);
%         vD_phi_map=vD_coef*vD_phi_map./(Btot_XZ_map.^3);
vD_X_map=vD_X_map./(Btot_XZ_map.^3);
vD_Z_map=vD_Z_map./(Btot_XZ_map.^3);
vD_phi_map=vD_phi_map./(Btot_XZ_map.^3);

%         vD_X_map=smooth_small_map(vD_X_map);
%         vD_Z_map=smooth_small_map(vD_Z_map);
%         vD_phi_map=smooth_small_map(vD_phi_map);



%             B_data=reshape(Btot_XZ_map(:,:)',sizeX*sizeZ,1);
%             Btot_PR_map_interp=dtinterp(XZ_mesh,XZ_mesh_dtri,B_data,RR,Z_PR_map(:,1:size_r),PR_DT_INTERPOLATION_METHOD);

vD_X_PR_map=-Btor_PR_map(:,1:size_r).*gradB_Z_PR_map(:,1:size_r);
vD_Z_PR_map=Btor_PR_map(:,1:size_r).*gradB_X_PR_map(:,1:size_r);
vD_phi_PR_map=BpolX_PR_map.*gradB_Z_PR_map(:,1:size_r)-BpolZ_PR_map.*gradB_X_PR_map(:,1:size_r);

vD_X_PR_map=vD_X_PR_map./(Btot_PR_map.^3);
vD_Z_PR_map=vD_Z_PR_map./(Btot_PR_map.^3);
vD_phi_PR_map=vD_phi_PR_map./(Btot_PR_map.^3);

%             v_data=reshape(vD_X_map(:,:)',sizeX*sizeZ,1);
%             vD_X_PR_map=dtinterp(XZ_mesh,XZ_mesh_dtri,v_data,RR,Z_PR_map(:,1:size_r),PR_DT_INTERPOLATION_METHOD);
%
%             v_data=reshape(vD_Z_map(:,:)',sizeX*sizeZ,1);
%             vD_Z_PR_map=dtinterp(XZ_mesh,XZ_mesh_dtri,v_data,RR,Z_PR_map(:,1:size_r),PR_DT_INTERPOLATION_METHOD);
%
%             v_data=reshape(vD_phi_map(:,:)',sizeX*sizeZ,1);
%             vD_phi_PR_map=dtinterp(XZ_mesh,XZ_mesh_dtri,v_data,RR,Z_PR_map(:,1:size_r),PR_DT_INTERPOLATION_METHOD);
