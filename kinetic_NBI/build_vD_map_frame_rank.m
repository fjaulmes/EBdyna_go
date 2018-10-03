ikink_kinetic_energy_phi=zeros(NB_PHI,1);

% big arrays, be careful with them
vD_X_PR_map_phi=zeros(NB_PHI,NB_THETA,size_r);
vD_Z_PR_map_phi=zeros(NB_PHI,NB_THETA,size_r);
vD_phi_PR_map_phi=zeros(NB_PHI,NB_THETA,size_r);

% bX_XZ_map(isnan(bX_XZ_map))=0;
% bZ_XZ_map(isnan(bZ_XZ_map))=0;
% bphi_XZ_map(isnan(bphi_XZ_map))=0;


for phi_rank=1:NB_PHI
    phi_rank
    B_PR_map=squeeze(Btot_map_phi(phi_rank,:,:));
    bX_PR_map=squeeze(bX_map_phi(phi_rank,:,:));
    bZ_PR_map=squeeze(bZ_map_phi(phi_rank,:,:));
    BpolX_PR_map=bX_PR_map.*B_PR_map;
    BpolZ_PR_map=bZ_PR_map.*B_PR_map;

%     B_data=reshape(B_PR_map(:,1:NR),NB_THETA*NR,1);
%     B_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,B_PR_map(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
%     B_XZ_map(isnan(B_XZ_map))=0;
%     B_XZ_map=B_XZ_map';
    Bdata=reshape(B_PR_map,NP*size_r,1);
    B_XZ_map=gridfit(finesse_data_X,finesse_data_Z,Bdata,scale_X,scale_Z,'smoothness',0.5);
    B_XZ_map=B_XZ_map';

    BX_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,BpolX_PR_map(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
    BX_XZ_map(isnan(BX_XZ_map))=0;
    BX_XZ_map=BX_XZ_map';
    BZ_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,BpolZ_PR_map(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
    BZ_XZ_map(isnan(BZ_XZ_map))=0;
    BZ_XZ_map=BZ_XZ_map';
    
    calculate_kink_kinetic_energy;
    ikink_kinetic_energy_phi(phi_rank)=phi_rank_kinetic_energy;

    gradB_X=zeros(sizeX,sizeZ);
    gradB_Z=zeros(sizeX,sizeZ);
    
    for (x=3:sizeX-2)
%         for (z=)
%             if radial_XZsmall_map(x,z)<Nradial-1
                gradB_X(x,3:sizeZ-2)=(1/12)*(-B_XZ_map(x+2,3:sizeZ-2)+B_XZ_map(x-2,3:sizeZ-2))+(2/3)*(B_XZ_map(x+1,3:sizeZ-2)-B_XZ_map(x-1,3:sizeZ-2));
                gradB_Z(x,3:sizeZ-2)=(1/12)*(-B_XZ_map(x,(3:sizeZ-2)+2)+B_XZ_map(x,(3:sizeZ-2)-2))+(2/3)*(B_XZ_map(x,(3:sizeZ-2)+1)-B_XZ_map(x,(3:sizeZ-2)-1));
%             end
%         end
    end
    %remove outer derivative
    gradB_X=gradB_X.*radial_mask;
    gradB_Z=gradB_Z.*radial_mask;
    
    gradB_X(isnan(gradB_X))=0;
    gradB_Z(isnan(gradB_Z))=0;
    
    gradB_X=gradB_X/DX;
    gradB_Z=gradB_Z/DX;
    
    vD_X_XZ_map=-Bphi_XZsmall_map.*gradB_Z;
    vD_Z_XZ_map=Bphi_XZsmall_map.*gradB_X;
    vD_phi_XZ_map=BX_XZ_map.*gradB_Z-BZ_XZ_map.*gradB_X;
    
    vD_X_XZ_map=vD_X_XZ_map./(B_XZ_map.^3);
    vD_Z_XZ_map=vD_Z_XZ_map./(B_XZ_map.^3);
    vD_phi_XZ_map=vD_phi_XZ_map./(B_XZ_map.^3);
    
    vD_X_XZ_map(isnan(vD_X_XZ_map))=0;
    vD_Z_XZ_map(isnan(vD_Z_XZ_map))=0;
    vD_phi_XZ_map(isnan(vD_phi_XZ_map))=0;

    
    vD_data=reshape(vD_X_XZ_map(:,:)',sizeX*sizeZ,1);
    vD_X_PR_map=griddata(X_scale_data,Z_scale_data,vD_data,X_PR_map(:,1:size_r),Z_PR_map(:,1:size_r),'cubic');
    vD_data=reshape(vD_Z_XZ_map(:,:)',sizeX*sizeZ,1);
    vD_Z_PR_map=griddata(X_scale_data,Z_scale_data,vD_data,X_PR_map(:,1:size_r),Z_PR_map(:,1:size_r),'cubic');
    vD_data=reshape(vD_phi_XZ_map(:,:)',sizeX*sizeZ,1);
    vD_phi_PR_map=griddata(X_scale_data,Z_scale_data,vD_data,X_PR_map(:,1:size_r),Z_PR_map(:,1:size_r),'cubic');
    vD_X_PR_map(isnan(vD_X_PR_map))=0;
    vD_Z_PR_map(isnan(vD_Z_PR_map))=0;
    vD_phi_PR_map(isnan(vD_phi_PR_map))=0;
    vD_X_PR_map(:,size_r-5:end)=0;
    vD_Z_PR_map(:,size_r-5:end)=0;
    vD_phi_PR_map(:,size_r-5:end)=0;

    vD_X_PR_map_phi(phi_rank,:,:)=vD_X_PR_map;
    vD_Z_PR_map_phi(phi_rank,:,:)=vD_Z_PR_map;
    vD_phi_PR_map_phi(phi_rank,:,:)=vD_phi_PR_map;
    
end