


gPhi_X=zeros(sizeX,sizeZ);
gPhi_Z=zeros(sizeX,sizeZ);
EX_XZ_map=zeros(sizeX,sizeZ);
EZ_XZ_map=zeros(sizeX,sizeZ);

% smooth_potential_map=smooth_small_map(Epot_XZ_map);

VZ=(3:sizeZ-2);
for (x=3:sizeX-2)
    
    gPhi_X(x,VZ)=(1/12)*(-Epot_XZ_map(x+2,VZ)+Epot_XZ_map(x-2,VZ))+(2/3)*(Epot_XZ_map(x+1,VZ)-Epot_XZ_map(x-1,VZ));
    gPhi_Z(x,VZ)=(1/12)*(-Epot_XZ_map(x,VZ+2)+Epot_XZ_map(x,VZ-2))+(2/3)*(Epot_XZ_map(x,VZ+1)-Epot_XZ_map(x,VZ-1));
    
end
% for (x=3:sizeX-2)
%     for (z=3:sizeZ-2)
%         gPhi_X(x,z)=(1/12)*(-Epot_XZ_map(x+2,z)+Epot_XZ_map(x-2,z))+(2/3)*(Epot_XZ_map(x+1,z)-Epot_XZ_map(x-1,z));
%         gPhi_Z(x,z)=(1/12)*(-Epot_XZ_map(x,z+2)+Epot_XZ_map(x,z-2))+(2/3)*(Epot_XZ_map(x,z+1)-Epot_XZ_map(x,z-1));
% %         gpsi_Z(x,z)=0.5*(psi2D(x,z+1)-psi2D(x,z-1));
%     end
% end
% for (x=3:sizeX-2)
%     for (z=3:sizeZ-2)
%         gPhi_X(x,z)=(1/12)*(-smooth_potential_map(x+2,z)+smooth_potential_map(x-2,z))+(2/3)*(smooth_potential_map(x+1,z)-smooth_potential_map(x-1,z));
%         gPhi_Z(x,z)=(1/12)*(-smooth_potential_map(x,z+2)+smooth_potential_map(x,z-2))+(2/3)*(smooth_potential_map(x,z+1)-smooth_potential_map(x,z-1));
%     end
% end
gPhi_X=gPhi_X/DX;
gPhi_Z=gPhi_Z/DX;

EX_XZ_map=-gPhi_X;
EZ_XZ_map=-gPhi_Z;




E_data=reshape(EX_XZ_map(:,:)',sizeX*sizeZ,1);
EX_PR_map=dtinterp(XZ_mesh,XZ_mesh_dtri,E_data,RR,Z_PR_map(:,1:size_r),PR_DT_INTERPOLATION_METHOD);
EX_PR_map(isnan(EX_PR_map))=0;

E_data=reshape(EZ_XZ_map(:,:)',sizeX*sizeZ,1);
EZ_PR_map=dtinterp(XZ_mesh,XZ_mesh_dtri,E_data,RR,Z_PR_map(:,1:size_r),PR_DT_INTERPOLATION_METHOD);
EZ_PR_map(isnan(EZ_PR_map))=0;
