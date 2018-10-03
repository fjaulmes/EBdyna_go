[XX,ZZ] = meshgrid(X_scale,Z_scale);

F2_PR_map=(Btor_PR_map.*Rpos_PR_map).^2;
F2_profile=mean(F2_PR_map(1:NP-1,:),1);
P_prime=gradient(P_initial_profile,psi_scale);
F2_prime=gradient(F2_profile,psi_scale);
% for(n=2:Nradial-1)
%     P_prime(n)=(P_initial_profile(n+1)-P_initial_profile(n-1))/(psi_scale(n+1)-psi_scale(n-1));
%     F2_prime(n)=(F2_profile(n+1)-F2_profile(n-1))/(psi_scale(n+1)-psi_scale(n-1));
% end

% P_prime(1)=P_prime(2);
% F2_prime(1)=F2_prime(2);
% P_prime(Nradial)=P_prime(Nradial-1);
% F2_prime(Nradial)=F2_prime(Nradial-1);

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

% lap_data=reshape((BZ_PR_map(:,:)),NP*Nradial,1);
% BZ_XZ_map=griddata(finesse_data_X,finesse_data_Z,lap_data,XX,ZZ,'cubic');
% BZ_XZ_map(isnan(BZ_XZ_map)) = 0; 
% BZ_XZ_map=BZ_XZ_map';
% lap_data=reshape((psi_PR_map(:,:)),NP*Nradial,1);
% psi_XZ_map=griddata(finesse_data_X,finesse_data_Z,lap_data,XX,ZZ,'cubic');
% psi_XZ_map(isnan(psi_XZ_map)) = 0; 
% psi_XZ_map=psi_XZ_map';


BZ_XZ_map=BpolZ_initial_map;

lap_star_XZ_map_psi=zeros(NZ,NZ);
lap_star_XZ_map=zeros(NZ,NZ);

gpsi_Z=zeros(NZ,NZ);
RdR_gpsi_R_R=zeros(NZ,NZ);
dZ_gpsi_Z=zeros(NZ,NZ);

for (x=3:NZ-2)
    for (z=3:NZ-2)
        gpsi_Z(x,z)=(1/12)*(-psi_XZ_map(x,z+2)+psi_XZ_map(x,z-2))+(2/3)*(psi_XZ_map(x,z+1)-psi_XZ_map(x,z-1));
    end
end
gpsi_Z=gpsi_Z/DX;

for (x=3:NZ-2)
    for (z=3:NZ-2)
        RdR_gpsi_R_R(x,z)=(1/12)*(-BZ_XZ_map(x+2,z)+BZ_XZ_map(x-2,z))+(2/3)*(BZ_XZ_map(x+1,z)-BZ_XZ_map(x-1,z));
        dZ_gpsi_Z(x,z)=(1/12)*(-gpsi_Z(x,z+2)+gpsi_Z(x,z-2))+(2/3)*(gpsi_Z(x,z+1)-gpsi_Z(x,z-1));
%         gpsi_Z(x,z)=0.5*(psi2D(x,z+1)-psi2D(x,z-1));
    end
end
RdR_gpsi_R_R=RdR_gpsi_R_R/DX;
dZ_gpsi_Z=dZ_gpsi_Z/DX;

RdR_gpsi_R_R=Rpos_map.*RdR_gpsi_R_R;
lap_star_XZ_map_psi=RdR_gpsi_R_R+dZ_gpsi_Z;

lap_star_PR_map=-mu0*(Rpos_PR_map.^2).*Pprime_PR_map-FFprime_PR_map;

lap_data=reshape((lap_star_PR_map(:,:)),NP*Nradial,1);
lap_star_XZ_map=griddata(finesse_data_X,finesse_data_Z,lap_data,XX,ZZ,'cubic');
lap_star_XZ_map(isnan(lap_star_XZ_map)) = 0; 
lap_star_XZ_map=lap_star_XZ_map';


% Pprime=zeros(Nradial,1);
% FFprime=zeros(Nradial,1);
% 
% for (p=2:Nradial-1)
%     Pprime(p)=(P_initial_profile(p+1)-P_initial_profile(p-1))/(psi_scale(p+1)-psi_scale(p-1));  
%     FFprime(p)=(F2_profile(p+1)-F2_profile(p-1))/(psi_scale(p+1)-psi_scale(p-1));
% end
% 
% Pprime=Pprime/psi_global;
% FFprime=FFprime/psi_global;

figure(4);
imagesc((lap_star_XZ_map_psi-lap_star_XZ_map)',[-2 2]);
axis xy;
colorbar;
