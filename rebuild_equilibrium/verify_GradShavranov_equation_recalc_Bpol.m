
% clear all

% BpolX_initial_map=BR_XZ_map;
% BpolZ_initial_map=BZ_XZ_map;



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

%NX=400;
%NZ=4*NX;
finesse_data_X=reshape((Rpos_PR_map(:,:)-R0),NP*Nradial,1);
finesse_data_Z=reshape(Z_PR_map(:,:),NP*Nradial,1);

for (x=1:NZ)
    for (z=1:NZ)
        Rpos_map(x,z)=R0+(x-mid_X-1)*DX;
    end
end



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


BR_XZ_map_recalc=BX_recalc;
BZ_XZ_map_recalc=BZ_recalc;
Bpol_XZ_map_recalc=sqrt(BR_XZ_map_recalc.^2+BZ_XZ_map_recalc.^2);

figure(7);
imagesc((Bpol_XZ_map-Bpol_XZ_map_recalc)',[-0.001 0.001])
axis xy

