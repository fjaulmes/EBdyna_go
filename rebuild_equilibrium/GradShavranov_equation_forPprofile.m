% approximate scaling
Pmax=Pmax_finesse*(2*alpha_finesse^2)/((a/R0))/mu0; 

[XX,ZZ] = meshgrid(X_scale,Z_scale);

F2_PR_map=(Btor_PR_map.*Rpos_PR_map).^2;
F2_profile=mean(F2_PR_map(1:NP-1,:),1);
P_prime=gradient(Pmax*P_initial_profile,psi_scale);
F2_prime=gradient(F2_profile,psi_scale);



for(p=1:NP)
    for(r=1:Nradial)
%         Pprime_PR_map(p,r)=P_prime(r);
        FFprime_PR_map(p,r)=0.5*F2_prime(r);
    end
end


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
gpsi_X=gpsi_X./Rpos_map;
for (x=3:NZ-2)
    for (z=3:NZ-2)
        RdR_gpsi_R_R(x,z)=(1/12)*(-gpsi_X(x+2,z)+gpsi_X(x-2,z))+(2/3)*(gpsi_X(x+1,z)-gpsi_X(x-1,z));
        dZ_gpsi_Z(x,z)=(1/12)*(-gpsi_Z(x,z+2)+gpsi_Z(x,z-2))+(2/3)*(gpsi_Z(x,z+1)-gpsi_Z(x,z-1));
%         gpsi_Z(x,z)=0.5*(psi2D(x,z+1)-psi2D(x,z-1));
    end
end
RdR_gpsi_R_R=RdR_gpsi_R_R/DX;
dZ_gpsi_Z=dZ_gpsi_Z/DX;

RdR_gpsi_R_R=Rpos_map.*RdR_gpsi_R_R;
lap_star_XZ_map_psi=RdR_gpsi_R_R+dZ_gpsi_Z;

Rlap_star_XZ_map_psi=(lap_star_XZ_map_psi./Rpos_map.^2)/mu0;

RFFprime_PR_map=(FFprime_PR_map./(Rpos_PR_map.^2))/mu0;
FF_data=reshape((RFFprime_PR_map(:,:)),NP*Nradial,1);
RFFprime_XZ_map=griddata(finesse_data_X,finesse_data_Z,FF_data,XX,ZZ,'cubic');
RFFprime_XZ_map(isnan(RFFprime_XZ_map)) = 0; 
RFFprime_XZ_map=RFFprime_XZ_map';
FFprime_XZ_map=mu0*RFFprime_XZ_map.*(Rpos_map.^2);

Rlap_star_FF_XZ_map_psi=Rlap_star_XZ_map_psi+RFFprime_XZ_map;

lap_data=reshape((Rlap_star_FF_XZ_map_psi(:,:)'),NZ*NZ,1);

%gridddata is not working here, so we need terrible dtinterp
[XX,ZZ] = meshgrid(X_scale,Z_scale);
X_scale_data=reshape(XX(:,:),NZ*NZ,1);
Z_scale_data=reshape(ZZ(:,:),NZ*NZ,1);

finesse_mesh=[finesse_data_X  finesse_data_Z];
[finesse_mesh, IDTRI, J] = unique(finesse_mesh,'last','rows');
finesse_mesh_dtri=DelaunayTri(finesse_mesh);
XZ_mesh=[X_scale_data Z_scale_data];
XZ_mesh_dtri=DelaunayTri(XZ_mesh);
lap_star_FF_PR_map_psi=dtinterp(XZ_mesh,XZ_mesh_dtri,lap_data,X_PR_map,Z_PR_map,'quadratic');

lap_star_FF_PR_map_psi(isnan(lap_star_FF_PR_map_psi)) = 0; 

lap_star_FF_profile=mean(lap_star_FF_PR_map_psi(1:NP-1,:),1);
lap_star_FF_profile(end)=0.5*lap_star_FF_profile(end-1);


P_prime=-(lap_star_FF_profile);
%if (lap_star_FF_profile>0)
%P_prime=max(P_prime,0);
%else
%P_prime=min(P_prime,0);
%end

P_profile=P_initial_profile*Pmax;
for (p=Nradial-2:-1:1)
    P_profile(p)=P_profile(p+1)+(P_prime(p+1))*(psi_scale(p)-psi_scale(p+1));
end
for(p=1:NP)
    for(r=1:Nradial)
        Pprime_PR_map(p,r)=P_prime(r);
    end
end

Pp_data=reshape((mu0*(Rpos_PR_map.^2).*Pprime_PR_map(:,:)),NP*Nradial,1);
RPprime_XZ_map=griddata(finesse_data_X,finesse_data_Z,Pp_data,XX,ZZ,'cubic');
RPprime_XZ_map(isnan(RPprime_XZ_map)) = 0; 
RPprime_XZ_map=RPprime_XZ_map';
Pprime_XZ_map= RPprime_XZ_map./Rpos_map.^2/mu0;

Dsize=round(NX);
P_fac=fminsearch(@(x) lap_psi_fun(x,Rlap_star_XZ_map_psi(mid_X-Dsize:mid_X+Dsize,mid_Z-Dsize:mid_Z+Dsize),...
    RFFprime_XZ_map(mid_X-Dsize:mid_X+Dsize,mid_Z-Dsize:mid_Z+Dsize),Pprime_XZ_map(mid_X-Dsize:mid_X+Dsize,mid_Z-Dsize:mid_Z+Dsize)),1)

disp('Value of pressure on axis according to solving of GS equation:')
P0=P_profile(1)*P_fac
P_profile=P_profile*P0/P_profile(1);

