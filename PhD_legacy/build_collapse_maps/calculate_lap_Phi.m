


gPhi_X=zeros(sizeX,sizeZ);
gPhi_Z=zeros(sizeX,sizeZ);
ggPhi_X=zeros(sizeX,sizeZ);
ggPhi_Z=zeros(sizeX,sizeZ);

% smooth_potential_map=smooth_small_map(Epot_XZ_map);

VZ=(3:sizeZ-2);
for (x=3:sizeX-2)
    
    gPhi_X(x,VZ)=(1/12)*(-Epot_XZ_map(x+2,VZ)+Epot_XZ_map(x-2,VZ))+(2/3)*(Epot_XZ_map(x+1,VZ)-Epot_XZ_map(x-1,VZ));
    gPhi_Z(x,VZ)=(1/12)*(-Epot_XZ_map(x,VZ+2)+Epot_XZ_map(x,VZ-2))+(2/3)*(Epot_XZ_map(x,VZ+1)-Epot_XZ_map(x,VZ-1));
    
end

gPhi_X=Rpos_XZsmall_map.*gPhi_X/DX;
gPhi_Z=gPhi_Z/DX;

for (x=3:sizeX-2)
    
    ggPhi_X(x,VZ)=(1/12)*(-gPhi_X(x+2,VZ)+gPhi_X(x-2,VZ))+(2/3)*(gPhi_X(x+1,VZ)-gPhi_X(x-1,VZ));
    ggPhi_Z(x,VZ)=(1/12)*(-gPhi_Z(x,VZ+2)+gPhi_Z(x,VZ-2))+(2/3)*(gPhi_Z(x,VZ+1)-gPhi_Z(x,VZ-1));
    
end

ggPhi_X=ggPhi_X./Rpos_XZsmall_map/DX;
ggPhi_Z=ggPhi_Z/DX;


E_data=reshape(ggPhi_X(:,:)',sizeX*sizeZ,1);
ggPhiX_PR_map=dtinterp(XZ_mesh,XZ_mesh_dtri,E_data,RR,Z_PR_map(:,1:size_r),PR_DT_INTERPOLATION_METHOD);
ggPhiX_PR_map(isnan(ggPhiX_PR_map))=0;

E_data=reshape(ggPhi_Z(:,:)',sizeX*sizeZ,1);
ggPhiZ_PR_map=dtinterp(XZ_mesh,XZ_mesh_dtri,E_data,RR,Z_PR_map(:,1:size_r),PR_DT_INTERPOLATION_METHOD);
ggPhiZ_PR_map(isnan(ggPhiZ_PR_map))=0;

ggPhitheta_PR_map=E_potential_PR_map_next_next-2*E_potential_PR_map+E_potential_PR_map_prev_prev;
ggPhitheta_PR_map=ggPhitheta_PR_map/(2*DOMEGA)/(2*DOMEGA);
ggPhitheta_PR_map=ggPhitheta_PR_map./(Rpos_PR_map(:,1:size_r).^2);

lap_Phi_PR_map=ggPhiX_PR_map+ggPhiZ_PR_map+ggPhitheta_PR_map;