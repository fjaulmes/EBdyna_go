


Btot_XZ_map=sqrt(BpolX_XZ_map.^2+BpolZ_XZ_map.^2+Bphi_XZsmall_map.^2);

RBphi_XZsmall_map=Rpos_XZsmall_map.*Bphi_XZsmall_map;

%gradBX_X=zeros(size_X,size_Z);
gradBX_Z=zeros(sizeX,sizeZ);
gradBZ_X=zeros(sizeX,sizeZ);
%gradBZ_Z=zeros(size_X,size_Z);
gradRBphi_X=zeros(sizeX,sizeZ);
gradBphi_Z=zeros(sizeX,sizeZ);

for (x=3:sizeX-2)
    for (z=3:sizeZ-2)
        if radial_XZsmall_map(x,z)<Nradial-1
            %gradBX_X(x,z)=(1/12)*(-BpolX_initial_XZsmall_map(x+2,z)+BpolX_initial_XZsmall_map(x-2,z))+(2/3)*(BpolX_initial_XZsmall_map(x+1,z)-BpolX_initial_XZsmall_map(x-1,z));
            gradBX_Z(x,z)=(1/12)*(-BpolX_XZ_map(x,z+2)+BpolX_XZ_map(x,z-2))+(2/3)*(BpolX_XZ_map(x,z+1)-BpolX_XZ_map(x,z-1));
            gradBZ_X(x,z)=(1/12)*(-BpolZ_XZ_map(x+2,z)+BpolZ_XZ_map(x-2,z))+(2/3)*(BpolZ_XZ_map(x+1,z)-BpolZ_XZ_map(x-1,z));
            %gradBZ_Z(x,z)=(1/12)*(-BpolZ_initial_XZsmall_map(x,z+2)+BpolZ_initial_XZsmall_map(x,z-2))+(2/3)*(BpolZ_initial_XZsmall_map(x,z+1)-BpolZ_initial_XZsmall_map(x,z-1));
            gradRBphi_X(x,z)=(1/12)*(-RBphi_XZsmall_map(x+2,z)+RBphi_XZsmall_map(x-2,z))+(2/3)*(RBphi_XZsmall_map(x+1,z)-RBphi_XZsmall_map(x-1,z));
            gradBphi_Z(x,z)=(1/12)*(-Bphi_XZsmall_map(x,z+2)+Bphi_XZsmall_map(x,z-2))+(2/3)*(Bphi_XZsmall_map(x,z+1)-Bphi_XZsmall_map(x,z-1));
        end
    end
end
%gradBX_X(isnan(gradBX_X))=0;
gradBX_Z(isnan(gradBX_Z))=0;
gradBZ_X(isnan(gradBZ_X))=0;
%gradBZ_Z(isnan(gradBZ_Z))=0;
gradRBphi_X(isnan(gradRBphi_X))=0;
gradBphi_Z(isnan(gradBphi_Z))=0;

%gradBX_X=gradBX_X/DX;
gradBX_Z=gradBX_Z/DX;
gradBZ_X=gradBZ_X/DX;
%gradBZ_Z=gradBZ_Z/DX;
gradRBphi_X=(gradRBphi_X./Rpos_XZsmall_map)/DX;
gradBphi_Z=gradBphi_Z/DX;

gradRBphi_X(isnan(gradRBphi_X))=0;


rotB_X=gradBphi_Z;
rotB_Z=-gradRBphi_X;
rotB_phi=gradBZ_X-gradBX_Z;



