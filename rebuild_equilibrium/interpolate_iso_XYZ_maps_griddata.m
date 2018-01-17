


radial_XZ_map=ones(NZ,NZ)*Nradial;

%----------------------------------
% Interpolate (X,Z) data from (P,R) data
%----------------------------------

% q_XZ_map=griddata(finesse_data(:,1),finesse_data(:,2),finesse_data(:,end),XX,ZZ,'cubic');
theta_XZ_map=griddata(finesse_data(:,1),finesse_data(:,2),theta_data,XX,ZZ,'cubic');
% psi_XZ_map=griddata(finesse_data(:,1),finesse_data(:,2),finesse_data(:,end-1),XX,ZZ,'cubic');
r_XZ_map=griddata(finesse_data(:,1),finesse_data(:,2),r_data,XX,ZZ,'cubic');
radial_XZ_map=griddata(finesse_data(:,1),finesse_data(:,2),radial_data,XX,ZZ,'cubic');
% Bphi_XZ_map=griddata(finesse_data(:,1),finesse_data(:,2),finesse_data(:,end-13),XX,ZZ,'cubic');
% BZ_XZ_map=griddata(finesse_data(:,1),finesse_data(:,2),finesse_data(:,end-19),XX,ZZ,'cubic');
% BR_XZ_map=griddata(finesse_data(:,1),finesse_data(:,2),finesse_data(:,end-25),XX,ZZ,'cubic');
P_XZ_map=griddata(finesse_data(:,1),finesse_data(:,2),finesse_data(:,end-31),XX,ZZ,'cubic');
dl_XZ_map=griddata(finesse_data(:,1),finesse_data(:,2),dl_data,XX,ZZ,'cubic');

dl_data=reshape((dl_X_PR_map(:,:)),NP*Nradial,1);
dl_X=griddata(finesse_data_X,finesse_data_Z,dl_data,XX,ZZ,'cubic');
dl_X(isnan(dl_X))=0;
dl_X=dl_X';

dl_data=reshape((dl_Z_PR_map(:,:)),NP*Nradial,1);
dl_Z=griddata(finesse_data_X,finesse_data_Z,dl_data,XX,ZZ,'cubic');
dl_Z(isnan(dl_Z))=0;
dl_Z=dl_Z';


% the transpose is done later in this file with the scaling
gf_data=reshape((q_PR_map(:,:)),NP*Nradial,1);
q_XZ_map=griddata(finesse_data_X,finesse_data_Z,gf_data,XX,ZZ,'cubic');
q_XZ_map(isnan(q_XZ_map))=0;
% q_XZ_map=q_XZ_map';
gf_data=reshape((psi_PR_map(:,:)),NP*Nradial,1);
psi_XZ_map=griddata(finesse_data_X,finesse_data_Z,gf_data,XX,ZZ,'cubic');
psi_XZ_map(isnan(psi_XZ_map))=0;
% psi_XZ_map=psi_XZ_map';
gf_data=reshape((Btor_PR_map(:,:)),NP*Nradial,1)/Baxis;
Bphi_XZ_map=griddata(finesse_data_X,finesse_data_Z,gf_data,XX,ZZ,'cubic');
Bphi_XZ_map(isnan(Bphi_XZ_map))=0;
% Bphi_XZ_map=Bphi_XZ_map';
gf_data=reshape((BX_PR_map(:,:)),NP*Nradial,1)/Baxis;
BR_XZ_map=griddata(finesse_data_X,finesse_data_Z,gf_data,XX,ZZ,'cubic');
BR_XZ_map(isnan(BR_XZ_map))=0;
% BR_XZ_map=BR_XZ_map';
gf_data=reshape((BZ_PR_map(:,:)),NP*Nradial,1)/Baxis;
BZ_XZ_map=griddata(finesse_data_X,finesse_data_Z,gf_data,XX,ZZ,'cubic');
BZ_XZ_map(isnan(BZ_XZ_map))=0;
% BZ_XZ_map=BZ_XZ_map';



% load gridfit_XZ_map.mat

q_data=reshape((q_PR_map(:,:)),NP*Nradial,1);

q_XZ_map_linear=griddata(finesse_data_X,finesse_data_Z,q_data,XX,ZZ,'linear');
q_XZ_map_linear(isnan(q_XZ_map_linear)) = 0; 
q_XZ_map_linear=q_XZ_map_linear';

mask_XZ_map=ones(NZ,NZ).*(q_XZ_map_linear~=0);





q_XZ_map(isnan(q_XZ_map)) = 0; 
theta_XZ_map(isnan(theta_XZ_map)) = 0; 
psi_XZ_map(isnan(psi_XZ_map)) = 1; 
r_XZ_map(isnan(r_XZ_map)) = 0; 
radial_XZ_map(isnan(radial_XZ_map)) = Nradial; 
radial_XZ_map=min(radial_XZ_map,Nradial);
Bphi_XZ_map(isnan(Bphi_XZ_map)) = 0; 
BZ_XZ_map(isnan(BZ_XZ_map)) = 0; 
BR_XZ_map(isnan(BR_XZ_map)) = 0; 
dl_XZ_map(isnan(dl_XZ_map)) = 0; 
P_XZ_map(isnan(P_XZ_map)) = 0; 

if SIGN_TOROIDAL_FIELD==1
    [value X2_out]=max(psi_XZ_map_linear(mid_X:end,Z_axis_pos));
else
    [value X2_out]=min(psi_XZ_map_linear(mid_X:end,Z_axis_pos));
end
X2_out=X2_out+mid_X-1


Bpol_XZ_map=sqrt(BR_XZ_map.^2+BZ_XZ_map.^2);


%----------------------------------
% Rescaling the maps
%----------------------------------
% no need for griddata
% psi_XZ_map = psi_global*(psi_XZ_map-1); 
% psi_XZ_map = psi_XZ_map+(-psi_global-min(min(psi_XZ_map)));
% psi_XZ_map=psi_XZ_map';
 


q_XZ_map = q_XZ_map'; 
psi_XZ_map = psi_XZ_map'; 
r_XZ_map = a*r_XZ_map'; 
radial_XZ_map = radial_XZ_map'; 
Bphi_XZ_map = Baxis*Bphi_XZ_map'; 
BZ_XZ_map = Baxis*BZ_XZ_map'; 
BR_XZ_map = Baxis*BR_XZ_map'; 
dl_XZ_map = a*dl_XZ_map'; 
Bpol_XZ_map = Baxis*Bpol_XZ_map'; 
theta_XZ_map=theta_XZ_map';

F_XZ_map=Bphi_XZ_map.*(Rpos_map);


%----------------------------------
% Cleaning the theta_map
%----------------------------------
theta_axis=zeros(1,NR);
for (x=1:X_axis_pos)
    [eps theta_axis(x)]=min(abs(pi-theta_XZ_map(x,:)));
    if (eps == pi)
        theta_axis(x)=0;
    end
end
X1_out=min(find(theta_axis))
Z_up_half=max(theta_axis);

theta_XZ_map=theta_XZ_map+(mask_XZ_map-1)*(-20);
for (x=X_axis_pos:X2_out)
    [eps theta_axis(x)]=min(abs(theta_XZ_map(x,:)));
end
theta_axis(X_axis_pos-1)=Z_axis_pos;
theta_axis(X_axis_pos)=Z_axis_pos;
theta_axis(X_axis_pos+1)=Z_axis_pos;
theta_XZ_map=theta_XZ_map+(mask_XZ_map-1)*(20);

[theta_axis_max X_theta_axis_max]=max(theta_axis);
% [theta_axis_max X_theta_axis_max]=max(theta_axis(1:X_theta_axis_max-1));
% for (x=X_theta_axis_max:X2_out)
%     if theta_axis(x+1)>theta_axis(x)
%         theta_axis(x+1)=theta_axis(x);
%     end
% end

for(x=X1_out:X2_out)
    Z_up_half=theta_axis(x);
    for (z=1:Z_up_half-1)
        if(theta_XZ_map(x,z)~=0)
            theta_XZ_map(x,z)=2*pi-theta_XZ_map(x,z);
        end
    end
end


if SIGN_TOROIDAL_FIELD<0
    theta_XZ_map=2*pi-theta_XZ_map;
end

    
theta_XZ_map=theta_XZ_map.*mask_XZ_map;

%rescaling the pressure
% Pmax_finesse=max(max(P_XZ_map));
% n=1;
% for(p=1:NP)
%     for(r=1:Nradial)
%         pressure_PR_map(p,r)=finesse_data(n,end-31);
%         n=n+1;
%     end
% end
P_initial_profile=mean(pressure_PR_map(1:NP-1,:),1);
Pmax_finesse=max(max(P_initial_profile));

GradShavranov_equation_forPprofile;
% lap_star_XZ_map_psi=RdR_gpsi_R_R+dZ_gpsi_Z;

Ne0=P0/(2*Te0)

% pressure_PR_map=P0*(pressure_PR_map/Pmax_finesse);
P_initial_profile=P_profile;
for(p=1:NP)
    for(r=1:Nradial)
        pressure_PR_map(p,r)=P_initial_profile(r);
    end
end
P_prime=gradient(P_initial_profile,psi_scale);
for(p=1:NP)
    for(r=1:Nradial)
        Pprime_PR_map(p,r)=P_prime(r);
    end
end

lap_star_PR_map=-mu0*(Rpos_PR_map.^2).*Pprime_PR_map-FFprime_PR_map;
lap_data=reshape((lap_star_PR_map(:,:)),NP*Nradial,1);
% lap_star_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,lap_data(IDTRI),XX,ZZ,'quadratic');
lap_star_XZ_map=griddata(finesse_data_X,finesse_data_Z,lap_data,XX,ZZ,'cubic');
lap_star_XZ_map(isnan(lap_star_XZ_map)) = 0; 
lap_star_XZ_map=lap_star_XZ_map';

close all

figure(5);
imagesc((lap_star_XZ_map_psi-lap_star_XZ_map)',[-2 2]);
axis xy;
colorbar;


pause(0.4)

P_XZ_map = (P0/Pmax_finesse)*P_XZ_map'; 
