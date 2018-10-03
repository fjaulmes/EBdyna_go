
gradB_X=zeros(sizeX,sizeZ);
gradB_Z=zeros(sizeX,sizeZ);

% smooth_B_map=smooth_small_map(Btot_XZ_map);



for (x=3:sizeX-2)
    for (z=3:sizeZ-2)
        gradB_X(x,z)=(1/12)*(-Btot_XZ_map(x+2,z)+Btot_XZ_map(x-2,z))+(2/3)*(Btot_XZ_map(x+1,z)-Btot_XZ_map(x-1,z));
        gradB_Z(x,z)=(1/12)*(-Btot_XZ_map(x,z+2)+Btot_XZ_map(x,z-2))+(2/3)*(Btot_XZ_map(x,z+1)-Btot_XZ_map(x,z-1));
    end
end
gradB_X=gradB_X/DX;
gradB_Z=gradB_Z/DX;



% B_data=reshape(gradB_X(:,:)',sizeX*sizeZ,1);
% gradB_X_PR_map=dtinterp(XZ_mesh,XZ_mesh_dtri,B_data,RR,Z_PR_map(:,1:size_r),PR_DT_INTERPOLATION_METHOD);
% gradB_X_PR_map(isnan(gradB_X_PR_map))=0;
% 
% B_data=reshape(gradB_Z(:,:)',sizeX*sizeZ,1);
% gradB_Z_PR_map=dtinterp(XZ_mesh,XZ_mesh_dtri,B_data,RR,Z_PR_map(:,1:size_r),PR_DT_INTERPOLATION_METHOD);
% gradB_Z_PR_map(isnan(gradB_Z_PR_map))=0;
% 
%         
% 
        