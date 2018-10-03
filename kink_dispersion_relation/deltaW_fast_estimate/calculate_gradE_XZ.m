
gradE_X=zeros(sizeX,sizeZ);
gradE_Z=zeros(sizeX,sizeZ);

% smooth_B_map=smooth_small_map(Btot_XZ_map);



for (x=3:sizeX-2)
    for (z=3:sizeZ-2)
        gradE_X(x,z)=(1/12)*(-Etot_XZ_map(x+2,z)+Etot_XZ_map(x-2,z))+(2/3)*(Etot_XZ_map(x+1,z)-Etot_XZ_map(x-1,z));
        gradE_Z(x,z)=(1/12)*(-Etot_XZ_map(x,z+2)+Etot_XZ_map(x,z-2))+(2/3)*(Etot_XZ_map(x,z+1)-Etot_XZ_map(x,z-1));
    end
end
gradE_X=gradE_X/DX;
gradE_Z=gradE_Z/DX;



% B_data=reshape(gradE_X(:,:)',sizeX*sizeZ,1);
% gradE_X_PR_map=dtinterp(XZ_mesh,XZ_mesh_dtri,B_data,RR,Z_PR_map(:,1:size_r),PR_DT_INTERPOLATION_METHOD);
% gradE_X_PR_map(isnan(gradE_X_PR_map))=0;
% 
% B_data=reshape(gradE_Z(:,:)',sizeX*sizeZ,1);
% gradE_Z_PR_map=dtinterp(XZ_mesh,XZ_mesh_dtri,B_data,RR,Z_PR_map(:,1:size_r),PR_DT_INTERPOLATION_METHOD);
% gradE_Z_PR_map(isnan(gradE_Z_PR_map))=0;
% 
%         
% 
        