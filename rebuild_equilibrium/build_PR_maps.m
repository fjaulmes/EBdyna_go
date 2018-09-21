

theta_PR_map=zeros(NP,Nradial);
theta_PR_map_gf=theta_PR_map;
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
		% Because we have flipped up down the Z axis, BX sign is reversed!!
		% we thus lose the direct theta angle (Bpol<0) but it is better like this
		% (poloidal flux pointing upward)
		% q should have been negative but we take positive convention.
        BX_PR_map(p,r)=SIGN_TOROIDAL_FIELD*finesse_data(n,end-25);
        BZ_PR_map(p,r)=finesse_data(n,end-19);
        dBXdX_PR_map(p,r)=SIGN_TOROIDAL_FIELD*finesse_data(n,end-24);
		dBXdZ_PR_map(p,r)=SIGN_TOROIDAL_FIELD*finesse_data(n,end-23);
        dBZdX_PR_map(p,r)=finesse_data(n,end-18);
        dBZdZ_PR_map(p,r)=finesse_data(n,end-17);
        Bpol_PR_map(p,r)=sqrt(finesse_data(n,end-19)^2+finesse_data(n,end-25)^2);
        Btor_PR_map(p,r)=SIGN_CO_CURRENT_FIELD*finesse_data(n,end-13);
        pressure_PR_map(p,r)=finesse_data(n,end-31);
        if (theta_scale(p)<=pi)
            theta_data(n)=theta_scale(p);
			theta_PR_map_gf(p,r)=theta_scale(p);
        else
            theta_data(n)=2*pi-theta_scale(p);
			theta_PR_map_gf(p,r)=2*pi-theta_scale(p);
        end
		theta_full_data(n)=theta_scale(p);
        cos_theta_data(n)=cos(theta_data(n));
        r_data(n)=sqrt((X_PR_map(p,r)-X_axis)^2+(Z_PR_map(p,r)-Z_axis)^2);
        radial_data(n)=r;
        n=n+1;
    end
end


% Building triangles for area calculation
% edges are named a,b,c,d and e for the diagonal
dl_a_X=zeros(NP-1,Nradial-1);
dl_a_Z=zeros(NP-1,Nradial-1);
dl_b_X=zeros(NP-1,Nradial-1);
dl_b_Z=zeros(NP-1,Nradial-1);
dl_c_X=zeros(NP-1,Nradial-1);
dl_c_Z=zeros(NP-1,Nradial-1);
dl_d_X=zeros(NP-1,Nradial-1);
dl_d_Z=zeros(NP-1,Nradial-1);
dl_e_X=zeros(NP-1,Nradial-1);
dl_e_Z=zeros(NP-1,Nradial-1);
Rpos1_values=zeros(NP-1,Nradial-1);
Zpos1_values=zeros(NP-1,Nradial-1);
Rpos2_values=zeros(NP-1,Nradial-1);
Zpos2_values=zeros(NP-1,Nradial-1);
r1_Xpos_values=zeros(NP-1,Nradial-1);
r1_Zpos_values=zeros(NP-1,Nradial-1);
r2_Xpos_values=zeros(NP-1,Nradial-1);
r2_Zpos_values=zeros(NP-1,Nradial-1);

for(p=1:NP-1)
    for(r=1:Nradial-1)
		dl_a_X(p,r)=X_PR_map(p,r)-X_PR_map(p+1,r);
		dl_a_Z(p,r)=Z_PR_map(p,r)-Z_PR_map(p+1,r);
		dl_b_X(p,r)=X_PR_map(p,r+1)-X_PR_map(p,r);
		dl_b_Z(p,r)=Z_PR_map(p,r+1)-Z_PR_map(p,r);
		dl_c_X(p,r)=X_PR_map(p+1,r+1)-X_PR_map(p,r+1);
		dl_c_Z(p,r)=Z_PR_map(p+1,r+1)-Z_PR_map(p,r+1);
		dl_d_X(p,r)=X_PR_map(p+1,r)-X_PR_map(p+1,r+1);
		dl_d_Z(p,r)=Z_PR_map(p+1,r)-Z_PR_map(p+1,r+1);
		dl_e_X(p,r)=X_PR_map(p,r+1)-X_PR_map(p+1,r);
		dl_e_Z(p,r)=Z_PR_map(p,r+1)-Z_PR_map(p+1,r);
		r1_Xpos_values(p,r)=0.5*(X_PR_map(p+1,r)+X_PR_map(p+1,r+1));
		r1_Zpos_values(p,r)=0.5*(Z_PR_map(p+1,r)+Z_PR_map(p+1,r+1));
		r2_Xpos_values(p,r)=0.5*(X_PR_map(p,r)+X_PR_map(p,r+1));
		r2_Zpos_values(p,r)=0.5*(Z_PR_map(p,r)+Z_PR_map(p,r+1));
		Rpos1_values(p,r)=R0+(X_PR_map(p,r)+X_PR_map(p,r+1)+X_PR_map(p+1,r))/3;
		Zpos1_values(p,r)=(Z_PR_map(p,r)+Z_PR_map(p,r+1)+Z_PR_map(p+1,r))/3;
		Rpos2_values(p,r)=R0+(X_PR_map(p,r+1)+X_PR_map(p+1,r+1)+X_PR_map(p+1,r))/3;
		Zpos2_values(p,r)=(Z_PR_map(p,r+1)+Z_PR_map(p+1,r+1)+Z_PR_map(p+1,r))/3;
	end	
end

dl_a=sqrt(dl_a_X.^2+dl_a_Z.^2);
dl_b=sqrt(dl_b_X.^2+dl_b_Z.^2);
dl_c=sqrt(dl_c_X.^2+dl_c_Z.^2);
dl_d=sqrt(dl_d_X.^2+dl_d_Z.^2);
dl_e=sqrt(dl_e_X.^2+dl_e_Z.^2);
% half perimeter of triangles
dl_t1=0.5*(dl_a+dl_b+dl_e);
dl_t2=0.5*(dl_c+dl_d+dl_e);
%Areas of triangles
dA1=sqrt(dl_t1.*(dl_t1-dl_a).*(dl_t1-dl_b).*(dl_t1-dl_e));
dA2=sqrt(dl_t2.*(dl_t2-dl_c).*(dl_t2-dl_d).*(dl_t2-dl_e));

dSurf=dA1+dA2;
%barycenter of quadrilateral
Rpos_values=(Rpos1_values.*dA1+Rpos2_values.*dA2)./dSurf;
Xpos_values=Rpos_values-R0;
Zpos_values=(Zpos1_values.*dA1+Zpos2_values.*dA2)./dSurf;
r_values=sqrt((Xpos_values-X_axis).^2+(Zpos_values-Z_axis).^2);
r1_values=sqrt((r1_Xpos_values-X_axis).^2+(r1_Zpos_values-Z_axis).^2);
r2_values=sqrt((r2_Xpos_values-X_axis).^2+(r2_Zpos_values-Z_axis).^2);
angle_values=acos((r1_values.^2+r2_values.^2-dl_c.^2)./(2*r1_values.*r2_values));

% the angle calculations are far too approximate but we assume they give the right ratios
% we normalize them so that one turn is 2*pi
for(r=1:Nradial-1)
    angle_values(:,r)=2*pi*angle_values(:,r)/sum(angle_values(:,r));
end
% recalculating the surface elements according to curved contour
% WARNING : this is still too approximate and will b erescaled ultimately after surface integration
dSurf=0.5*r_values.*angle_values.*(dl_b+dl_d);

%elementary volumes
dVol=2*pi*Rpos_values.*dSurf;

Rpos_PR_map=R0+X_PR_map;


% maps for contour integration
dl_PR_map(2:end,2:end)=0.5*(dl_a+dl_c);
length_contour_profile=sum(dl_PR_map,1);
dr_PR_map=0.5*(dl_b+dl_d);
radial_PR_map=sqrt((X_PR_map-X_axis)^2+(Z_PR_map-Z_axis)^2);

p=NP;
for(r=2:Nradial-1)
    DX_axis(r)=0.5*(X_PR_map(p,r+1)-X_PR_map(p,r-1));
end






psi_global=abs(psi0_input)


%rescaling the magnetic fields
BX_PR_map=Baxis*BX_PR_map;
BZ_PR_map=Baxis*BZ_PR_map;
dBXdX_PR_map=Baxis*dBXdX_PR_map;
dBXdZ_PR_map=Baxis*dBXdZ_PR_map;
dBZdX_PR_map=Baxis*dBZdX_PR_map;
dBZdZ_PR_map=Baxis*dBZdZ_PR_map;
Bpol_PR_map=Baxis*Bpol_PR_map;
% sign of field taken here
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

