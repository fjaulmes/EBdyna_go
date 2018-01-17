
% clear all

load('..\data_tokamak\physics_constants.mat')
load('..\data_tokamak\pressure_profile.mat')
load('..\data_tokamak\B_fields.mat')
load('..\data_tokamak\flux_geometry.mat')
load('..\data_tokamak\tokamak_map_dimensions')
load('..\data_tokamak\tokamak_PR_map')


%%
[XX,ZZ] = meshgrid(X_scale,Z_scale);

F2_PR_map=(Btor_PR_map.*Rpos_PR_map).^2;
F2_profile=mean(F2_PR_map(1:NP-1,:),1);
P_prime=gradient(P_initial_profile,psi_scale);
F2_prime=gradient(F2_profile,psi_scale);

for(p=1:NP)
    for(r=1:Nradial)
        Pprime_PR_map(p,r)=P_prime(r);
        FFprime_PR_map(p,r)=0.5*F2_prime(r);
    end
end

NX=500;
NZ=4*NX;
finesse_data_X=reshape((Rpos_PR_map(:,:)-R0),NP*Nradial,1);
finesse_data_Z=reshape(Z_PR_map(:,:),NP*Nradial,1);

for (x=1:NZ)
    for (z=1:NZ)
        Rpos_map(x,z)=R0+(x-mid_X-1)*DX;
    end
end


BZ_XZ_map=BpolZ_initial_map;

lap_star_XZ_map_psi=zeros(NZ,NZ);
lap_star_XZ_map=zeros(NZ,NZ);

gpsi_X=zeros(NZ,NZ);
gpsi_Z=zeros(NZ,NZ);
RdR_gpsi_R_R=zeros(NZ,NZ);
dZ_gpsi_Z=zeros(NZ,NZ);

for (x=3:NZ-2)
    for (z=3:NZ-2)
        gpsi_X(x,z)=(1/12)*(-psi_XZ_map(x+2,z)+psi_XZ_map(x-2,z))+(2/3)*(psi_XZ_map(x+1,z)-psi_XZ_map(x-1,z));
        gpsi_Z(x,z)=(1/12)*(-psi_XZ_map(x,z+2)+psi_XZ_map(x,z-2))+(2/3)*(psi_XZ_map(x,z+1)-psi_XZ_map(x,z-1));
    end
end
gpsi_X=gpsi_X/DX;
gpsi_Z=gpsi_Z/DX;

% here you can verify the good consistency
% of the psi XZ map and the BZ map
BZ_recalc=gpsi_X./Rpos_map;
BX_recalc=-gpsi_Z./Rpos_map;

for (x=3:NZ-2)
    for (z=3:NZ-2)
        RdR_gpsi_R_R(x,z)=(1/12)*(-BZ_XZ_map(x+2,z)+BZ_XZ_map(x-2,z))+(2/3)*(BZ_XZ_map(x+1,z)-BZ_XZ_map(x-1,z));
        dZ_gpsi_Z(x,z)=(1/12)*(-gpsi_Z(x,z+2)+gpsi_Z(x,z-2))+(2/3)*(gpsi_Z(x,z+1)-gpsi_Z(x,z-1));
    end
end
RdR_gpsi_R_R=RdR_gpsi_R_R/DX;
dZ_gpsi_Z=dZ_gpsi_Z/DX;

RdR_gpsi_R_R=Rpos_map.*RdR_gpsi_R_R;
lap_star_XZ_map_psi=RdR_gpsi_R_R+dZ_gpsi_Z;


lap_star_PR_map=-mu0*(Rpos_PR_map.^2).*Pprime_PR_map-FFprime_PR_map;

lap_star_XZ_map_psi=max(lap_star_XZ_map_psi,min(min(lap_star_PR_map))-1);
lap_star_XZ_map_psi=min(lap_star_XZ_map_psi,max(max(lap_star_PR_map))+1);

lap_data=reshape((lap_star_PR_map(:,:)),NP*Nradial,1);
lap_star_XZ_map=griddata(finesse_data_X,finesse_data_Z,lap_data,XX,ZZ,'cubic');
lap_star_XZ_map(isnan(lap_star_XZ_map)) = 0; 
lap_star_XZ_map=lap_star_XZ_map';


figure(6);
imagesc(X_scale,Z_scale,((BpolX_initial_map-(BX_recalc))/mean(mean(abs(BpolX_initial_map))))',[-0.01 0.01]);
axis xy;
colorbar;
title('B_X FINESSE error')

figure(5);
imagesc(X_scale,Z_scale,((BZ_XZ_map-BZ_recalc)/mean(mean(abs(BZ_XZ_map))))',[-0.01 0.01]);
axis xy;
colorbar;
title('B_Z FINESSE error')

figure(4);
imagesc(X_scale,Z_scale,((lap_star_XZ_map_psi-lap_star_XZ_map)./lap_star_XZ_map_psi)',[-1 1]);
axis xy;
colorbar;
title('GS FINESSE error')
