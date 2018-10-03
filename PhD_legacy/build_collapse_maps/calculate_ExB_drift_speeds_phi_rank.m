
    filename='../B_maps/B0';
    filename=strcat(filename,frame_name,'.mat');
    load(filename)

    bX_PR(:,:)=bX_map_phi(phi_rank,:,:);
    bZ_PR(:,:)=bZ_map_phi(phi_rank,:,:);
    Btot_PR_map(:,:)=Btot_map_phi(phi_rank,:,:);
    BpolX_PR_map=bX_PR.*Btot_PR_map;
    BpolZ_PR_map=bZ_PR.*Btot_PR_map;    
    
E_data=reshape(Ephi_PR_map(:,1:size_r),NP*size_r,1);
Ephi_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,E_data(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
Ephi_XZ_map=Ephi_XZ_map';
Ephi_XZ_map(isnan(Ephi_XZ_map))=0;


% BpolX_PR_map=BHpol_X_PR_map_size_r+BstarX_PR_map;
% BpolZ_PR_map=BHpol_Z_PR_map_size_r+BstarZ_PR_map;

B_data=reshape(BpolX_PR_map(:,1:size_r),NP*size_r,1);
BpolX_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,B_data(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
BpolX_XZ_map=BpolX_XZ_map';
BpolX_XZ_map(isnan(BpolX_XZ_map))=0;

B_data=reshape(BpolZ_PR_map(:,1:size_r),NP*size_r,1);
BpolZ_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,B_data(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
BpolZ_XZ_map=BpolZ_XZ_map';
BpolZ_XZ_map(isnan(BpolZ_XZ_map))=0;

Btot_XZ_map=sqrt(BpolX_XZ_map.^2+BpolZ_XZ_map.^2+Bphi_XZsmall_map.^2);


vExB_X_map=EZ_XZ_map.*Bphi_XZsmall_map-Ephi_XZ_map.*BpolZ_XZ_map;
vExB_Z_map=Ephi_XZ_map.*BpolX_XZ_map-Bphi_XZsmall_map.*EX_XZ_map;
vExB_phi_map=EX_XZ_map.*BpolZ_XZ_map-BpolX_XZ_map.*EZ_XZ_map;
vExB_X_map=vExB_X_map./(Btot_XZ_map.^2);
vExB_Z_map=vExB_Z_map./(Btot_XZ_map.^2);
vExB_phi_map=vExB_phi_map./(Btot_XZ_map.^2);

%         vExB_X_map=smooth_small_map(vExB_X_map);
%         vExB_Z_map=smooth_small_map(vExB_Z_map);
%         vExB_phi_map=smooth_small_map(vExB_phi_map);



Btot_PR_map=sqrt(BpolX_PR_map.^2+BpolZ_PR_map.^2+Btor_PR_map(:,1:size_r).^2);

vExB_X_PR_map=EZ_PR_map.*Btor_PR_map(:,1:size_r)-Ephi_PR_map.*BpolZ_PR_map;
vExB_Z_PR_map=Ephi_PR_map.*BpolX_PR_map-Btor_PR_map(:,1:size_r).*EX_PR_map;
vExB_phi_PR_map=EX_PR_map.*BpolZ_PR_map-BpolX_PR_map.*EZ_PR_map;
vExB_X_PR_map=vExB_X_PR_map./(Btot_PR_map(:,1:size_r).^2);
vExB_Z_PR_map=vExB_Z_PR_map./(Btot_PR_map(:,1:size_r).^2);
vExB_phi_PR_map=vExB_phi_PR_map./(Btot_PR_map(:,1:size_r).^2);

%             v_data=reshape(vExB_X_map(:,:)',sizeX*sizeZ,1);
%             vExB_X_PR_map=dtinterp(XZ_mesh,XZ_mesh_dtri,v_data,RR,Z_PR_map(:,1:size_r),PR_DT_INTERPOLATION_METHOD);
%
%             v_data=reshape(vExB_Z_map(:,:)',sizeX*sizeZ,1);
%             vExB_Z_PR_map=dtinterp(XZ_mesh,XZ_mesh_dtri,v_data,RR,Z_PR_map(:,1:size_r),PR_DT_INTERPOLATION_METHOD);
%
%             v_data=reshape(vExB_phi_map(:,:)',sizeX*sizeZ,1);
%             vExB_phi_PR_map=dtinterp(XZ_mesh,XZ_mesh_dtri,v_data,RR,Z_PR_map(:,1:size_r),PR_DT_INTERPOLATION_METHOD);

