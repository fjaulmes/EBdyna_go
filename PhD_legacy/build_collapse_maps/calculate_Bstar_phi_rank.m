

BstarX_XZ_map=zeros(sizeX,sizeZ);
BstarZ_XZ_map=zeros(sizeX,sizeZ);

psi_PR_map=zeros(NB_THETA,NR);
psi_PR_map(:,1:size_r-4)=psi_star_PR_map(:,1:size_r-4);
psi_data=reshape(psi_PR_map(:,1:NR),NB_THETA*NR,1);
% psi_XZ_map=griddata(finesse_data_X_extended,finesse_data_Z_extended,psi_data,XX_small,ZZ_small,'cubic');
psi_XZ_map=dtinterp(finesse_mesh_extended,finesse_mesh_extended_dtri,psi_PR_map(IDTRI_EXT),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
psi_XZ_map(isnan(psi_XZ_map))=0;
psi_XZ_map=psi_XZ_map';
psi2D=psi_XZ_map;

psi_data=reshape(psi_star_PR_map(:,1:size_r),NB_THETA*size_r,1);
psi_star_XZ_map=gridfit(finesse_data_X,finesse_data_Z,psi_data,scale_X,scale_Z,'smoothness',0.5);
psi_star_XZ_map=psi_star_XZ_map';
psi2D=psi_star_XZ_map;


%%
gpsi_R=zeros(sizeX,sizeZ);
gpsi_Z=zeros(sizeX,sizeZ);
VZ=(3:sizeZ-2);
for (x=3:sizeX-2)
    
    gpsi_R(x,VZ)=(1/12)*(-psi2D(x+2,VZ)+psi2D(x-2,VZ))+(2/3)*(psi2D(x+1,VZ)-psi2D(x-1,VZ));
    gpsi_Z(x,VZ)=(1/12)*(-psi2D(x,VZ+2)+psi2D(x,VZ-2))+(2/3)*(psi2D(x,VZ+1)-psi2D(x,VZ-1));
    
end
% for (x=3:sizeX-2)
%     for (z=3:sizeZ-2)
%         gpsi_R(x,z)=(1/12)*(-psi2D(x+2,z)+psi2D(x-2,z))+(2/3)*(psi2D(x+1,z)-psi2D(x-1,z));
%         gpsi_Z(x,z)=(1/12)*(-psi2D(x,z+2)+psi2D(x,z-2))+(2/3)*(psi2D(x,z+1)-psi2D(x,z-1));
%     end
% end
gpsi_R=gpsi_R/DX;
gpsi_Z=gpsi_Z/DX;

BstarX_XZ_map=-gpsi_Z./Rpos_XZsmall_map;
BstarZ_XZ_map=gpsi_R./Rpos_XZsmall_map;

Bstar_data=reshape(BstarX_XZ_map(:,:)',sizeX*sizeZ,1);
BstarX_PR_map=griddata(X_scale_data,Z_scale_data,Bstar_data,RR,Z_PR_map(:,1:size_r),'cubic');

Bstar_data=reshape(BstarZ_XZ_map(:,:)',sizeX*sizeZ,1);
BstarZ_PR_map=griddata(X_scale_data,Z_scale_data,Bstar_data,RR,Z_PR_map(:,1:size_r),'cubic');
        


