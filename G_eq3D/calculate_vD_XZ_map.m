gradB_X=zeros(size_X,size_Z);
gradB_Z=zeros(size_X,size_Z);

for (x=3:size_X-2)
    for (z=3:size_Z-2)
        if radial_XZsmall_map(x,z)<Nradial-1
            gradB_X(x,z)=(1/12)*(-Btot_XZ_map(x+2,z)+Btot_XZ_map(x-2,z))+(2/3)*(Btot_XZ_map(x+1,z)-Btot_XZ_map(x-1,z));
            gradB_Z(x,z)=(1/12)*(-Btot_XZ_map(x,z+2)+Btot_XZ_map(x,z-2))+(2/3)*(Btot_XZ_map(x,z+1)-Btot_XZ_map(x,z-1));
        end
    end
end
gradB_X(isnan(gradB_X))=0;
gradB_Z(isnan(gradB_Z))=0;

gradB_X=gradB_X/DX;
gradB_Z=gradB_Z/DX;

vD_X_XZ_map=-Bphi_XZsmall_map.*gradB_Z;
vD_Z_XZ_map=Bphi_XZsmall_map.*gradB_X;
vD_phi_XZ_map=BpolX_initial_XZsmall_map.*gradB_Z-BpolZ_initial_XZsmall_map.*gradB_X;

vD_X_XZ_map=vD_X_XZ_map./(Btot_XZ_map).^3;
vD_Z_XZ_map=vD_Z_XZ_map./(Btot_XZ_map).^3;
vD_phi_XZ_map=vD_phi_XZ_map./(Btot_XZ_map).^3;