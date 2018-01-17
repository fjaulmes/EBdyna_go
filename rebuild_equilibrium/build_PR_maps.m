

theta_PR_map=zeros(NP,Nradial);
X_PR_map=zeros(NP,Nradial);
Z_PR_map=zeros(NP,Nradial);
inv_Rpos_PR_map=zeros(NP,Nradial);
q_PR_map=zeros(NP,Nradial);
psi_PR_map=zeros(NP,Nradial);
dl_PR_map=zeros(NP,Nradial);
dl_X_PR_map=zeros(NP,Nradial);
dl_Z_PR_map=zeros(NP,Nradial);
dr_PR_map=zeros(NP,Nradial);
dr_X_PR_map=zeros(NP,Nradial);
dr_Z_PR_map=zeros(NP,Nradial);
cos_ki_PR_map=zeros(NP,Nradial);
BX_PR_map=zeros(NP,Nradial);
BZ_PR_map=zeros(NP,Nradial);
dBXdX_PR_map=zeros(NP,Nradial);
dBXdZ_PR_map=zeros(NP,Nradial);
dBZdX_PR_map=zeros(NP,Nradial);
dBZdZ_PR_map=zeros(NP,Nradial);
Bpol_PR_map=zeros(NP,Nradial);
Btor_PR_map=zeros(NP,Nradial);
pressure_PR_map=zeros(NP,Nradial);
Pprime_PR_map=zeros(NP,Nradial);
FFprime_PR_map=zeros(NP,Nradial);
radial_PR_map=zeros(NP,Nradial);

Bpol_surf=zeros(1,Nradial);
psi_profile=zeros(1,Nradial);
r_scale=zeros(1,Nradial);
q_profile=zeros(1,Nradial);



%DP=(2*pi)/(NP-1);
%for(p=1:NP)
%    theta_scale(p)=(p-1)*DP;
%end

n=1;
for(p=1:NP)
    for(r=1:Nradial)
        X_PR_map(p,r)=finesse_data(n,1);
        Z_PR_map(p,r)=finesse_data(n,2);
        q_PR_map(p,r)=abs(finesse_data(n,end));
		% this relation should always yield psi=0 at the separatrix
        psi_PR_map(p,r)=SIGN_CO_CURRENT_FIELD*finesse_data(n,end-1)+0.5*(SIGN_TOROIDAL_FIELD-1);
        theta_PR_map(p,r)=theta_scale(p);
        BX_PR_map(p,r)=finesse_data(n,end-25);
        BZ_PR_map(p,r)=SIGN_CO_CURRENT_FIELD*finesse_data(n,end-19);
        dBXdX_PR_map(p,r)=finesse_data(n,end-24);
		dBXdZ_PR_map(p,r)=finesse_data(n,end-23);
        dBZdX_PR_map(p,r)=SIGN_CO_CURRENT_FIELD*finesse_data(n,end-18);
        dBZdZ_PR_map(p,r)=SIGN_CO_CURRENT_FIELD*finesse_data(n,end-17);
        Bpol_PR_map(p,r)=sqrt(finesse_data(n,end-19)^2+finesse_data(n,end-25)^2);
        Btor_PR_map(p,r)=SIGN_CO_CURRENT_FIELD*finesse_data(n,end-13);
        pressure_PR_map(p,r)=finesse_data(n,end-31);
        if (theta_scale(p)<=pi)
            theta_data(n)=theta_scale(p);
        else
            theta_data(n)=2*pi-theta_scale(p);
        end
        cos_theta_data(n)=cos(theta_data(n));
        r_data(n)=sqrt((X_PR_map(p,r)-X_axis)^2+(Z_PR_map(p,r)-Z_axis)^2);
        radial_data(n)=r;
        n=n+1;
    end
end



Rpos_PR_map=R0+X_PR_map;


% map for contour integration


for(p=2:NP-1)
    for(r=2:Nradial-1)
        dl_avg_X=X_PR_map(p+1,r)-X_PR_map(p-1,r);
        dl_avg_Z=Z_PR_map(p+1,r)-Z_PR_map(p-1,r);
        dl_avg=0.5*sqrt(dl_avg_X^2+dl_avg_Z^2);
        
        dr_avg_X=X_PR_map(p,r+1)-X_PR_map(p,r-1);
        dr_avg_Z=Z_PR_map(p,r+1)-Z_PR_map(p,r-1);
        dr_avg=0.5*sqrt(dr_avg_X^2+dr_avg_Z^2);
        
        cos_ki_avg=0.25*(dl_avg_X*dr_avg_X+dl_avg_Z*dr_avg_Z);
        dl_PR_map(p,r)=dl_avg;
        dl_X_PR_map(p,r)=0.5*dl_avg_X/dl_avg;
        dl_Z_PR_map(p,r)=0.5*dl_avg_Z/dl_avg;
        dr_PR_map(p,r)=dr_avg;
		dr_X_PR_map(p,r)=0.5*dr_avg_X/dr_avg;
        dr_Z_PR_map(p,r)=0.5*dr_avg_Z/dr_avg;
        cos_ki_PR_map(p,r)=cos_ki_avg/(dl_avg*dr_avg);
        radial_PR_map(p,r)=sqrt((X_PR_map(p,r)-X_axis)^2+(Z_PR_map(p,r)-Z_axis)^2);
    end
end
p=1;
for(r=2:Nradial-1)
	dl_avg_X=X_PR_map(p+1,r)-X_PR_map(NP-1,r);
	dl_avg_Z=Z_PR_map(p+1,r)-Z_PR_map(NP-1,r);
	dl_avg=0.5*sqrt(dl_avg_X^2+dl_avg_Z^2);
	
	dr_avg_X=X_PR_map(p,r+1)-X_PR_map(p,r-1);
	dr_avg_Z=Z_PR_map(p,r+1)-Z_PR_map(p,r-1);
	dr_avg=0.5*sqrt(dr_avg_X^2+dr_avg_Z^2);
	
	cos_ki_avg=0.25*(dl_avg_X*dr_avg_X+dl_avg_Z*dr_avg_Z);
	dl_PR_map(p,r)=dl_avg;
	dl_X_PR_map(p,r)=0.5*dl_avg_X/dl_avg;
	dl_Z_PR_map(p,r)=0.5*dl_avg_Z/dl_avg;
	dr_PR_map(p,r)=dr_avg;
	dr_X_PR_map(p,r)=0.5*dr_avg_X/dr_avg;
	dr_Z_PR_map(p,r)=0.5*dr_avg_Z/dr_avg;
	cos_ki_PR_map(p,r)=cos_ki_avg/(dl_avg*dr_avg);
    radial_PR_map(p,r)=sqrt((X_PR_map(p,r)-X_axis)^2+(Z_PR_map(p,r)-Z_axis)^2);
end
length_contour_profile=sum(dl_PR_map,1);
p=NP;
for(r=2:Nradial-1)
    dl_avg_X=X_PR_map(2,r)-X_PR_map(p-1,r);
    dl_avg_Z=Z_PR_map(2,r)-Z_PR_map(p-1,r);
    dl_avg=0.5*sqrt(dl_avg_X^2+dl_avg_Z^2);
    
    dr_avg_X=X_PR_map(p,r+1)-X_PR_map(p,r-1);
    dr_avg_Z=Z_PR_map(p,r+1)-Z_PR_map(p,r-1);
    dr_avg=0.5*sqrt(dr_avg_X^2+dr_avg_Z^2);
    
    cos_ki_avg=0.25*(dl_avg_X*dr_avg_X+dl_avg_Z*dr_avg_Z);
    dl_X_PR_map(p,r)=0.5*dl_avg_X/dl_avg;
    dl_Z_PR_map(p,r)=0.5*dl_avg_Z/dl_avg;
    
    dl_PR_map(p,r)=dl_avg;
    dr_PR_map(p,r)=dr_avg;
 	dr_X_PR_map(p,r)=0.5*dr_avg_X/dr_avg;
    dr_Z_PR_map(p,r)=0.5*dr_avg_Z/dr_avg;
    cos_ki_PR_map(p,r)=cos_ki_avg/(dl_avg*dr_avg);
    DX_axis(r)=0.5*(X_PR_map(p,r+1)-X_PR_map(p,r-1));
    radial_PR_map(p,r)=sqrt((X_PR_map(p,r)-X_axis)^2+(Z_PR_map(p,r)-Z_axis)^2);
end


% code below seems to be obsolete: 
% calculate some psi value according to BZ flux....
% no need to!

% finesse_mesh=[finesse_data(:,1) finesse_data(:,2)];
% [finesse_mesh, IDTRI, J] = unique(finesse_mesh,'last','rows');
% finesse_mesh_dtri=DelaunayTri(finesse_mesh);
% data=SIGN_TOROIDAL_FIELD*finesse_data(:,end-19);
% BZ_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,data(IDTRI),XX,ZZ,'quadratic');
% BZ_XZ_map(isnan(BZ_XZ_map)) = 0; 
% BZ_XZ_map=BZ_XZ_map';
% BZ_XZ_map=Baxis*BZ_XZ_map;
% 
% BZ_axis=zeros(NZ,1);
% BZ_axis=BZ_XZ_map(:,Z_axis_pos);
% integ_BZ=zeros(NZ,1);
% for(x=2:NZ)
%     integ_BZ(x)=0.5*(BZ_axis(x)+BZ_axis(x-1))*((R0+(x-1)*DX-(mid_X-1)*DX)+0.5*DX);
% end
% psi_axis=zeros(NZ,1);
% for(x=2:NZ)
%     psi_axis(x)=psi_axis(x-1)+integ_BZ(x)*DX;
% end
% if (SIGN_TOROIDAL_FIELD>0)
%     psi0_recalc=min(psi_axis)
% else
%     psi0_recalc=max(psi_axis)
% end
% 
% data=SIGN_TOROIDAL_FIELD*finesse_data(:,end-1)-0.5*(SIGN_TOROIDAL_FIELD-1);
% % psi_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,data(IDTRI),XX,ZZ,'quadratic');
% psi_XZ_map=griddata(finesse_data(:,1),finesse_data(:,2),data,XX,ZZ,'cubic');
% if (SIGN_TOROIDAL_FIELD>0)
% psi_XZ_map(isnan(psi_XZ_map))=1;
% else
% psi_XZ_map(isnan(psi_XZ_map))=0;
% end
% if (SIGN_TOROIDAL_FIELD>0)
%     psi_XZ_map=abs(psi0_recalc)*(psi_XZ_map)';
% else
%     psi_XZ_map=abs(psi0_recalc)*(psi_XZ_map)';
% end
% optimize_psi2D_fromBZ;
% 
% 
% psi_XZ_map=psi_fac*psi_XZ_map;
% psi0_recalc=psi0_recalc*psi_fac;



psi_global=abs(psi0_input)


%rescaling the magnetic fields
BX_PR_map=Baxis*BX_PR_map;
BZ_PR_map=Baxis*BZ_PR_map;
dBXdX_PR_map=Baxis*dBXdX_PR_map;
dBXdZ_PR_map=Baxis*dBXdZ_PR_map;
dBZdX_PR_map=Baxis*dBZdX_PR_map;
dBZdZ_PR_map=Baxis*dBZdZ_PR_map;
Bpol_PR_map=Baxis*Bpol_PR_map;
% bug in sign of field in finesse output
Btor_PR_map=SIGN_TOROIDAL_FIELD*abs(Baxis)*abs(Btor_PR_map);

Bpol_surf=mean(Bpol_PR_map(1:NP-1,:),1);
% psi_profile=mean(psi_PR_map(1:NP-1,:),1);
q_profile=mean(q_PR_map(1:NP-1,:),1);


%rescaling the poloidal flux
if (SIGN_TOROIDAL_FIELD>0)
    % psi_profile=psi_global*(psi_profile-1);
    psi_PR_map = psi_global*(psi_PR_map-1);
else
    %psi_profile=psi_global*(psi_profile);
    psi_PR_map = psi_global*(psi_PR_map);
end

P_initial_profile=mean(pressure_PR_map(1:NP-1,:),1);


% rescaling quantities to physical units

psi_scale=mean(psi_PR_map(1:NP-1,:),1);
% the names should be unified but this an heritage from many years ago....
psi_profile=psi_scale;


finesse_data_X=reshape((Rpos_PR_map(:,:)-R0),NP*Nradial,1);
finesse_data_Z=reshape(Z_PR_map(:,:),NP*Nradial,1);
finesse_data_X_gf=reshape((Rpos_PR_map(:,1:NP-1)-R0),(NP-1)*Nradial,1);
finesse_data_Z_gf=reshape(Z_PR_map(:,1:NP-1),(NP-1)*Nradial,1);

calculate_gradients_r_theta;

n=1;
for(p=1:NP)
    for(r=1:Nradial)
        dl_data(n)=dl_PR_map(p,r);
        n=n+1;
    end
end


if (n==number_of_data_points+1)
    disp('data arrays correctly processed....');
end

%----------------------------------
% distance between flux surfaces in (P,R) coordinates
%----------------------------------

for(p=1:NP)
    for(r=2:Nradial)
        dr2=sqrt((X_PR_map(p,r)-X_PR_map(p,r-1))^2+(Z_PR_map(p,r)-Z_PR_map(p,r-1))^2);
        dist_surf_PR_map(p,r)=dr2;
    end
    dist_surf_PR_map(p,1)=0;
end

pause(0.1)

figure(2)
F_PR_map=(Btor_PR_map.*Rpos_PR_map);
F2_profile=mean(F_PR_map(1:end-1,:).^2,1);
BZ_profile=mean(BZ_PR_map(1:end-1,:),1);
subplot(2,1,1)
title('F2')
grid on;hold on
plot(psi_scale, F2_profile);
subplot(2,1,2)
title('BZ')
grid on;hold on
plot(psi_scale ,BZ_PR_map(round(0.5*(NP-1)),:));
plot( psi_scale ,BZ_PR_map(1,:),'r');

