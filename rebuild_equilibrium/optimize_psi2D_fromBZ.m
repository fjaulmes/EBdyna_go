% load('tokamak_map_dimensions.mat')
% load('B_fields.mat')
% BZ_XZ_map=BpolZ_initial_map;

% BZ_XZ_map=BZ_XZ_map.*mask_XZ_map;
psi2D=psi_XZ_map;

BZ=zeros(NZ,NZ);
gpsi_R=zeros(NZ,NZ);
gpsi_Z=zeros(NZ,NZ);

for (x=3:NZ-2)
    for (z=3:NZ-2)
        gpsi_R(x,z)=(1/12)*(-psi2D(x+2,z)+psi2D(x-2,z))+(2/3)*(psi2D(x+1,z)-psi2D(x-1,z));
        gpsi_Z(x,z)=(1/12)*(-psi2D(x,z+2)+psi2D(x,z-2))+(2/3)*(psi2D(x,z+1)-psi2D(x,z-1));
%         gpsi_Z(x,z)=0.5*(psi2D(x,z+1)-psi2D(x,z-1));
    end
end
gpsi_R=gpsi_R/DX;
gpsi_Z=gpsi_Z/DX;

BZ=gpsi_R./Rpos_map;

psi_fac=fminsearch(@(x) psi_fun(x,BZ,BZ_XZ_map),1)


figure(4);
imagesc((BZ-BZ_XZ_map)',[-0.05 0.05]);
axis xy;
colorbar;
