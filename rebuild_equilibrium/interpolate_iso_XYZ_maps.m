
GRIDFIT_SMOOTHNESS=0.1

radial_XZ_map=ones(NZ,NZ)*Nradial;

%----------------------------------
% Interpolate (X,Z) data from (P,R) data
%----------------------------------

% q_XZ_map=griddata(finesse_data(:,1),finesse_data(:,2),finesse_data(:,end),XX,ZZ,'cubic');
theta_XZ_map=griddata(finesse_data(:,1),finesse_data(:,2),theta_data,XX,ZZ,'cubic');
theta_XZ_map(isnan(theta_XZ_map)) = 0; 
theta_XZ_map=theta_XZ_map';

theta_full_XZ_map=griddata(finesse_data(:,1),finesse_data(:,2),theta_full_data,XX,ZZ,'cubic');
theta_full_XZ_map(isnan(theta_full_XZ_map)) = 0; 
theta_full_XZ_map=theta_full_XZ_map';
if SIGN_CO_CURRENT_FIELD<0
    theta_full_XZ_map=2*pi-theta_full_XZ_map;
end
% psi_XZ_map=griddata(finesse_data(:,1),finesse_data(:,2),finesse_data(:,end-1),XX,ZZ,'cubic');
r_XZ_map=griddata(finesse_data(:,1),finesse_data(:,2),r_data,XX,ZZ,'cubic');
%radial_XZ_map=griddata(finesse_data(:,1),finesse_data(:,2),radial_data,XX,ZZ,'cubic');


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

%gf_data=reshape(theta_PR_map_gf(:,1:NP-1),Nradial*(NP-1),1);
%theta_XZ_map_gridfit=gridfit(finesse_data_X_gf,finesse_data_Z_gf,gf_data,X_scale,Z_scale,'smoothness',GRIDFIT_SMOOTHNESS);
%theta_XZ_map_gridfit(isnan(theta_XZ_map_gridfit)) = 0; 
%theta_XZ_map_gridfit=theta_XZ_map_gridfit';

%theta_XZ_map_gridfit=max(theta_XZ_map_gridfit,0);
%theta_XZ_map_gridfit=min(theta_XZ_map_gridfit,pi);



% gf_data=reshape(q_PR_map,Nradial*NP,1);
% q_XZ_map=gridfit(finesse_data_X,finesse_data_Z,gf_data,X_scale,Z_scale,'smoothness',0.5);

gf_data=reshape(psi_PR_map(:,1:NP-1),Nradial*(NP-1),1);

% re-normalizing psi for precision
if SIGN_CO_CURRENT_FIELD<0
    gf_data = (gf_data/max(gf_data)); 
else
    gf_data = -(gf_data/min(gf_data)); 
end


psi_XZ_map_gridfit=gridfit(finesse_data_X_gf,finesse_data_Z_gf,gf_data,X_scale,Z_scale,'smoothness',GRIDFIT_SMOOTHNESS);
psi_XZ_map_gridfit(isnan(psi_XZ_map_gridfit)) = 1; 
psi_XZ_map_gridfit = psi_XZ_map_gridfit';

% innacuracy induced by interpolation needs to be fixed
if SIGN_CO_CURRENT_FIELD>0 
	psi_XZ_map=psi_XZ_map_gridfit-min(min(psi_XZ_map_gridfit))-1;
else
	psi_XZ_map=psi_XZ_map_gridfit-max(max(psi_XZ_map_gridfit))+1;
end


% The variable SIGN_CO_CURRENT_FIELD is used
% in the case of I<0 and Bt>0, to reserve the positive Bt 
% in the case SIGN_CO_CURRENT_FIELD>0, it is not possible to derive 
% a configuration with Bt>0 without swapping the Z axis 
% which triggered more difficulties than it solved
% so now we just keep Bt<0 in the case of aligned field and current
psi_XZ_map = psi_global*(psi_XZ_map);


%if SIGN_CO_CURRENT_FIELD<0
%    psi_XZ_map = psi_global*(psi_XZ_map); 
%else
%    psi_XZ_map = -psi_global*(psi_XZ_map); 
%end

gf_data=reshape(Btor_PR_map,Nradial*NP,1)/Baxis;
Bphi_XZ_map=gridfit(finesse_data_X,finesse_data_Z,gf_data,X_scale,Z_scale,'smoothness',GRIDFIT_SMOOTHNESS);

derpsiX=psi_XZ_map*0;
derpsiZ=psi_XZ_map*0;

for z=3:length(Z_scale)-2
    for x=3:length(X_scale)-2
        derpsiX(x,z)=(1/12)*(-psi_XZ_map(x+2,z)+psi_XZ_map(x-2,z))+(2/3)*(psi_XZ_map(x+1,z)-psi_XZ_map(x-1,z));
        derpsiZ(x,z)=(1/12)*(-psi_XZ_map(x,z+2)+psi_XZ_map(x,z-2))+(2/3)*(psi_XZ_map(x,z+1)-psi_XZ_map(x,z-1));
        derpsiX(x,z)=derpsiX(x,z)/DX;
        derpsiZ(x,z)=derpsiZ(x,z)/DX;       
    end
end

BR_XZ_map=-derpsiZ./Rpos_map;
BZ_XZ_map=derpsiX./Rpos_map;


% consitency between XZ and PR maps
% Is this really relevant ?

[XX ZZ]=meshgrid(X_scale,Z_scale);
sizeX=size(BR_XZ_map,1);
sizeZ=size(BR_XZ_map,2);
X_scale_data=reshape(XX,sizeX*sizeZ,1);
Z_scale_data=reshape(ZZ,sizeX*sizeZ,1);

B_data=reshape(BR_XZ_map(:,:)',sizeX*sizeZ,1);
BR_PR_map_recalc=griddata(X_scale_data,Z_scale_data,B_data,X_PR_map,Z_PR_map,'cubic');
BR_PR_map_recalc(isnan(BR_PR_map_recalc))=0;
BX_PR_map_recalc=BR_PR_map_recalc;

BX_PR_map=0.5*(BX_PR_map_recalc+BX_PR_map);

B_data=reshape(BZ_XZ_map(:,:)',sizeX*sizeZ,1);
BZ_PR_map_recalc=griddata(X_scale_data,Z_scale_data,B_data,X_PR_map,Z_PR_map,'cubic');
BZ_PR_map_recalc(isnan(BZ_PR_map_recalc))=0;

BZ_PR_map=0.5*(BZ_PR_map_recalc+BZ_PR_map);

%gf_data=reshape(BX_PR_map,Nradial*NP,1)/Baxis;
%BR_XZ_map=gridfit(finesse_data_X,finesse_data_Z,gf_data,X_scale,Z_scale,'smoothness',0.1);
%gf_data=reshape(BZ_PR_map,Nradial*NP,1)/Baxis;
%BZ_XZ_map=gridfit(finesse_data_X,finesse_data_Z,gf_data,X_scale,Z_scale,'smoothness',0.1);


q_data=reshape((q_PR_map(:,:)),NP*Nradial,1);

q_XZ_map_linear=griddata(finesse_data_X,finesse_data_Z,q_data,XX,ZZ,'cubic');
q_XZ_map_linear(isnan(q_XZ_map_linear)) = 0; 
q_XZ_map_linear=q_XZ_map_linear';

mask_XZ_map=ones(NZ,NZ).*(q_XZ_map_linear~=0);

%if SIGN_CO_CURRENT_FIELD<0
%    mask_XZ_map = ones(NZ,NZ).*(psi_XZ_map>=0);
%else
%    mask_XZ_map = ones(NZ,NZ).*(psi_XZ_map<=0);
%end






r_XZ_map(isnan(r_XZ_map)) = 0; 

Bphi_XZ_map(isnan(Bphi_XZ_map)) = 0; 
BZ_XZ_map(isnan(BZ_XZ_map)) = 0; 
BR_XZ_map(isnan(BR_XZ_map)) = 0; 
dl_XZ_map(isnan(dl_XZ_map)) = 0; 
P_XZ_map(isnan(P_XZ_map)) = 0; 

if SIGN_CO_CURRENT_FIELD==1
    [value X2_out]=max(psi_XZ_map_linear(mid_X:end,Z_axis_pos));
else
    [value X2_out]=min(psi_XZ_map_linear(mid_X:end,Z_axis_pos));
end
X2_out=X2_out+mid_X-1


Bpol_XZ_map=sqrt(BR_XZ_map.^2+BZ_XZ_map.^2);


%----------------------------------
% Rescaling the maps
%----------------------------------

% psi_XZ_map = psi_XZ_map+(-psi_global-min(min(psi_XZ_map)));

 
max_psi_scl=max(psi_scale);
min_psi_scl=min(psi_scale);
q_XZ_map=interp1(psi_scale,q_profile,max(min(psi_XZ_map,max_psi_scl),min_psi_scl));
q_XZ_map(isnan(q_XZ_map)) = 0; 

% q_XZ_map = q_XZ_map'; 

r_XZ_map = a*r_XZ_map'; 
% radial_XZ_map = radial_XZ_map'; 
Bphi_XZ_map = Baxis*Bphi_XZ_map'; 
%BZ_XZ_map = Baxis*BZ_XZ_map'; 
%BR_XZ_map = Baxis*BR_XZ_map'; 
dl_XZ_map = a*dl_XZ_map'; 
%Bpol_XZ_map = Baxis*Bpol_XZ_map'; 




%theta_XZ_map_gridfit=theta_XZ_map_gridfit.*mask_XZ_map;
%if SIGN_CO_CURRENT_FIELD<0
%    theta_XZ_map_gridfit=2*pi-theta_XZ_map_gridfit;
%end

%theta_XZ_map=0.5*(theta_XZ_map+theta_XZ_map_gridfit);

theta_XZ_map=max(theta_XZ_map,0);
theta_XZ_map=min(theta_XZ_map,pi);
theta_XZ_map=theta_XZ_map.*mask_XZ_map;

F_XZ_map=Bphi_XZ_map.*(Rpos_map);


%----------------------------------
% Cleaning the theta_map
%----------------------------------
theta_axis=zeros(1,NR);
theta_axis2=zeros(1,NR);
DELTA_THETA=0.001*pi;
for (x=1:X_axis_pos)
    [eps theta_axis(x)]=min(abs(pi-DELTA_THETA-theta_XZ_map(x,:)));
	theta_axis(x)=theta_axis(x);
    if (eps == pi-DELTA_THETA)
        theta_axis(x)=0;
    end
	[eps theta_axis2(x)]=min(abs(pi-theta_XZ_map(x,:)));
	%theta_axis(x)=Z_axis_pos-DELTA_AXIS+theta_axis(x);
    if (eps == pi)
        theta_axis2(x)=0;
    end
    if (theta_axis(x) ~= theta_axis2(x))
        theta_axis(x)=max(theta_axis(x),theta_axis2(x));
    end
end
X1_out=min(find(theta_axis>0))


%theta_XZ_map=theta_XZ_map+(mask_XZ_map-1)*(-20);
for (x=X_axis_pos:X2_out)
    %[eps theta_axis(x)]=min(abs(theta_XZ_map(x,:)));
	theta_axis(x)=Z_axis_pos;
end
theta_axis(X_axis_pos-1)=Z_axis_pos;
theta_axis(X_axis_pos)=Z_axis_pos;
theta_axis(X_axis_pos+1)=Z_axis_pos;
%theta_XZ_map=theta_XZ_map+(mask_XZ_map-1)*(+20);

if SIGN_CO_CURRENT_FIELD<0
[ theta_axis_max  X_theta_axis_max ]=max(theta_axis(X1_out+3:X2_out));
X_theta_axis_max=X_theta_axis_max+X1_out+2
else
[ theta_axis_max  X_theta_axis_max ]=min(theta_axis(X1_out+3:X2_out));
X_theta_axis_max=X_theta_axis_max+X1_out+2
end
% [theta_axis_max X_theta_axis_max]=max(theta_axis(1:X_theta_axis_max-1));
% for (x=X_theta_axis_max:X2_out)
%     if theta_axis(x+1)>theta_axis(x)
%         theta_axis(x+1)=theta_axis(x);
%     end
% end

for(x=X1_out:X_axis_pos-1)
    for (z=1:NZ)
        if(theta_full_XZ_map(x,z)~=0)
            theta_XZ_map(x,z)=theta_full_XZ_map(x,z);
        end
    end
end

for(x=X_axis_pos:X2_out)
    Z_up_half=theta_axis(x);
	if SIGN_CO_CURRENT_FIELD<0
    for (z=1:Z_up_half-1)
        if(theta_XZ_map(x,z)~=0)
            theta_XZ_map(x,z)=2*pi-theta_XZ_map(x,z);
        end
    end
	else
    for (z=Z_up_half+1:NZ)
        if(theta_XZ_map(x,z)~=0)
            theta_XZ_map(x,z)=2*pi-theta_XZ_map(x,z);
        end
    end
	end
end


if SIGN_CO_CURRENT_FIELD<0
    theta_XZ_map=2*pi-theta_XZ_map;
end

    
theta_XZ_map=theta_XZ_map.*mask_XZ_map;


save thetamap2 theta_XZ_map theta_full_XZ_map theta_axis theta_axis2

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
lap_star_XZ_map=griddata(finesse_data_X,finesse_data_Z,lap_data,XX,ZZ,'cubic');
lap_star_XZ_map(isnan(lap_star_XZ_map)) = 0; 
lap_star_XZ_map=lap_star_XZ_map';

P_XZ_map = (P0/Pmax_finesse)*P_XZ_map'; 

% the radial_XZ_map is used solely for the Epot mask
% and for display purpose

radial_XZ_map=griddata(finesse_data(:,1),finesse_data(:,2),radial_data,XX,ZZ,'cubic');

radial_XZ_map(isnan(radial_XZ_map)) = Nradial; 
radial_XZ_map=min(radial_XZ_map,Nradial);
radial_XZ_map = radial_XZ_map';

