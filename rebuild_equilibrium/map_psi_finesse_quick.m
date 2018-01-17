
% Evaluating the magnetic equilibrium configuration
% according to output from the FINESSE equilibrium code
% Thursday January 9th 2014

clc

clear all
close all
format compact

disp('... Initializing paramters...');
initialize_folder_names;
create_physics_constants;
Nradial=257;
NP=257;
DP=(2*pi)/(NP-1);
for(p=1:NP)
    theta_scale(p)=(p-1)*DP;    
end

create_tokamak_parameters;
initialize_XZ_maps;

for (x=1:NR)
    for (z=1:NZ)
        Rpos_map(x,z)=Rpos(x);
    end
end

disp('Building PR maps .... ');
build_PR_maps;

pause(0.4);



psi_data=reshape((psi_PR_map(:,:)),NP*Nradial,1);
psi_XZ_map_linear=griddata(finesse_data_X,finesse_data_Z,psi_data,XX,ZZ,'linear');
psi_XZ_map_linear(isnan(psi_XZ_map_linear)) = 1; 
psi_XZ_map_linear=psi_XZ_map_linear';

%%

disp('Interpolate XZ maps from PR data .... ');
interpolate_iso_XYZ_maps_griddata;

close all

%%
%----------------------------------
% maps for toroidal current density calculations
%----------------------------------
J_PHI_profile=(-P_prime.*R0^2-0.5*F2_prime/mu0)/R0;
J_PHI_PR_map=-Rpos_PR_map.*Pprime_PR_map-(FFprime_PR_map/mu0)./Rpos_PR_map;
J_PHI_flat_pressure_PR_map=-FFprime_PR_map/mu0/R0;
J_PHI_flat_profile=(-0.5*F2_prime/mu0)/R0;
JPHI_XZ_map=zeros(NZ,NZ);
JPHI_XZ_map=griddata(finesse_data_X,finesse_data_Z,reshape(J_PHI_PR_map,1,NP*Nradial),XX,ZZ,'cubic');
JPHI_XZ_map(find(isnan(JPHI_XZ_map)))=0;
JPHI_XZ_map=JPHI_XZ_map';
JPHI_flat_pressure_XZ_map=zeros(NZ,NZ);
JPHI_flat_pressure_XZ_map=griddata(finesse_data_X,finesse_data_Z,reshape(J_PHI_flat_pressure_PR_map,1,NP*Nradial),XX,ZZ,'cubic');
JPHI_flat_pressure_XZ_map=JPHI_flat_pressure_XZ_map';

q_data=reshape((q_PR_map(:,:)),NP*Nradial,1);
q_XZ_map_linear=griddata(finesse_data_X,finesse_data_Z,q_data,XX,ZZ,'linear');
q_XZ_map_linear(isnan(q_XZ_map_linear)) = 0; 
q_XZ_map_linear=q_XZ_map_linear';

mask_XZ_map=ones(NZ,NZ).*(q_XZ_map_linear~=0);



% ***************************************************************
% additional names
% ***************************************************************

psi2D=psi_XZ_map_linear;
% psi2D=0.9*psi_XZ_map+0.1*psi_XZ_map_linear;
% [Z_axis_pos val]=min(min(psi2D));

pos_axis=round(X_axis/DX)+mid_X;
Rpos_axis=Rpos(pos_axis);
Raxis_pos=pos_axis;

P_map=P_XZ_map*P0;
F_2_map=F_XZ_map.^2;

Btor=Bphi_XZ_map;

% Bpol=Bpol_XZ_map;
% bpol_X=BR_XZ_map;
% bpol_Z=BZ_XZ_map;




close all
create_tokamak_map_dimensions;



% *************************************************************
% defining curves approximating the shape of the flux surfaces
% *************************************************************

disp('find flux surface horizontal boundaries...');
find_flux_surface_horizontal_boundaries;

xi_psi_Nradial=zeros(1,Nradial);
xi_psi_Nradial(1)=X_axis;
disp('find displacement profile...');
fitting_displacement_profile

index=0;
%%
figure(2);
psi_XZ_map_lin=0.5*(psi_XZ_map+psi_XZ_map_linear);
Cont_psi2D=contour(X_scale,Z_scale,psi_XZ_map_lin',psi_profile(1:end));
Cont_psi2D_lin=Cont_psi2D;
%Cont_psi2D_lin=contour(X_scale,Z_scale,psi_XZ_map_linear',psi_profile(1:end));
Ncontours=find_number_of_contours(Cont_psi2D_lin);
for (n=1:Ncontours)
    index=index+1;
    psi_value_list_raw(n)=Cont_psi2D_lin(1,index);
    contour_length_raw(n)=Cont_psi2D_lin(2,index);
    for s=1:contour_length_raw(n)
        index=index+1;
        X_contour_raw(n,s)=Cont_psi2D_lin(1,index);
        Z_contour_raw(n,s)=Cont_psi2D_lin(2,index);
    end
end
[max_val LARGEST_CONTOUR]=max(contour_length_raw(1:end));
X_contour_raw_last=X_contour_raw(LARGEST_CONTOUR,:);
Z_contour_raw_last=Z_contour_raw(LARGEST_CONTOUR,:);

axis xy equal
xlabel('x')
ylabel('z')
title('psi2D');
hold on;


disp('Building Z curves of flux surface shapes ...');
Cont_psi2D=Cont_psi2D_lin;
extract_contour_information;

fitting_Zpsi_up_profiles;
fitting_Zpsi_down_profiles;


%%

evaluate_volume_flux_surfaces_Nradial;

% disp('Volume of vaccuum vessel = ');
% disp(volume_flux(Nradial));

Pmax=P0*volume_flux(Nradial)/(sum(P_initial_profile.*volume_flux));

radial_r_value_flux=sqrt(surf_flux_fit/pi);
for(p=1:NP)
    for(r=1:Nradial)
        r_PR_map(p,r)=radial_r_value_flux(r);
        r_data(n)=radial_r_value_flux(r);
    end
end

r_data=reshape(r_PR_map',NP*Nradial,1);
r_XZ_map=griddata(finesse_data(:,1),finesse_data(:,2),r_data,XX,ZZ,'cubic');
r_XZ_map(isnan(r_XZ_map)) = 0; 
r_XZ_map = r_XZ_map'; 


% *************************************************************
% defining the safety factor q=1 contour
% *************************************************************

build_q1_flux_surface;



% *************************************************************
% defining the helical field map and its normal vector
% *************************************************************


verify_GradShavranov_equation_recalc_Bpol;
Btot=sqrt(BR_XZ_map.^2+BZ_XZ_map.^2);

map_BH;



%TRIANGULATION_DISTORTION=1
%
Dpsi=0.5*(psi_scale(2)-psi_scale(1));
psi_scale_lin=(psi_scale(1):Dpsi:psi_scale(end)+10*Dpsi);

psiH_Nradial=zeros(1,Nradial);
psi_star_Nradial=zeros(NP,Nradial);
psi_star_Nradial_lin=zeros(NP,length(psi_scale_lin));
integ_q_profile_lin=zeros(NP,length(psi_scale_lin));


for (p=1:NP)
    integ_q_profile_lin(p,:)=interp1(psi_scale,(1-q_profile),psi_scale_lin);
%     integ_q_profile_lin(p,:)=interp1(psi_scale,(1-q_profile)./sqrt(1+(radial_PR_map(p,:)/Raxis).^2),psi_scale_lin);
%     integ_q_profile_lin(p,:)=interp1(psi_scale,(1-q_profile).*(1-0.5*(radial_PR_map(p,:)/Raxis).^2),psi_scale_lin);
end
% integ_q_profile_lin=interp1(psi_scale,(1-q_profile),psi_scale_lin);
% psi1_lin=interp1(integ_q_profile_lin,1:length(psi_scale_lin),1)
% inf_psi1=round(0.98*psi_rank_q1)

psi_star_Nradial_lin(:,1)=zeros(length(NP));

for (p=1:NP)
    for (x=2:length(integ_q_profile_lin))
        psi_star_Nradial_lin(p,x)=psi_star_Nradial_lin(p,x-1)+(integ_q_profile_lin(p,x)+integ_q_profile_lin(p,x-1))*0.5*Dpsi;
    end
    psi_star_Nradial(p,:)=interp1(psi_scale_lin,psi_star_Nradial_lin(p,:),psi_scale);
    psi_star_Nradial(p,end)=psi_star_Nradial(p,end-1)+0.5*(2-q_profile(end)-q_profile(end-1))*(psi_scale(end)-psi_scale(end-1));
    psi_star_Nradial(p,:)=psi_star_Nradial(p,:)-psi_star_Nradial(p,1);
end


psi_star_Nradial_avg=mean(psi_star_Nradial(1:NP-1,:),1);

% q_Nradial_recalc=1-gradient(psi_star_Nradial,psi_scale);


% psiH from psi star
psiH_Nradial=psi_scale-psi_star_Nradial_avg;



% mixing radius position
psi_mix_rank=interp1(psi_star_Nradial_avg(psi_rank_q1+1:end),(psi_rank_q1+1:Nradial),0)
psi_mix_value=interp1(psi_star_Nradial_avg(psi_rank_q1+1:end),psi_scale(psi_rank_q1+1:Nradial),0)


% for information
psiH_mix_value=interp1(psi_star_Nradial_avg(psi_rank_q1+1:end),psiH_Nradial(psi_rank_q1+1:Nradial),0);


% Bstar_PR_map_exact=Bpol_PR_map.*(1-q_PR_map).*(1-0.5*(radial_PR_map/R0).^2);
% Bstar_PR_map_exact=Bpol_PR_map.*(1-q_PR_map)./sqrt(1+(r_PR_map/Raxis).^2);
Bstar_PR_map_exact=Bpol_PR_map.*(1-q_PR_map);

psi_star_avg_PR_map=zeros(NP,1);
psiH_PR_map=zeros(NP,Nradial);

for (p=1:NP)
    for (r=2:Nradial)
        psi_star_avg_PR_map(p,r)=psi_star_Nradial_avg(r);
    end
end
psi_star_PR_map=psi_star_Nradial;
% psiH map from psi star map
psiH_PR_map=psi_PR_map-psi_star_PR_map;

% 
% psi_star_data=reshape((psi_star_PR_map(:,:)),NP*Nradial,1);
% psi_star_XZ_map=griddata(finesse_data_X,finesse_data_Z,psi_star_data,XX,ZZ,'cubic');
% psi_star_XZ_map(isnan(psi_star_XZ_map))=0;
% psi_star_XZ_map=psi_star_XZ_map';
%
% good interpolation for psi star XZ map
gf_data=reshape(psi_star_avg_PR_map(:,:),Nradial*NP,1);
if (SIGN_TOROIDAL_FIELD>0)
    gf_data=max(gf_data,-0.05);
else
    gf_data=min(gf_data,0.05);
end

psi_star_XZ_map=griddata(finesse_data_X,finesse_data_Z,gf_data,XX,ZZ,'cubic');
psi_star_XZ_map(isnan(psi_star_XZ_map))=0;
psi_star_XZ_map=psi_star_XZ_map';



% psiH map from psi star map
psiH_XZ_map=psi_XZ_map-psi_star_XZ_map;

calculate_Bstar_inital;
% recalculate_BH;

% the discrepancy between the two fields calculations
% is our error for helical flux estimate

% r_XZ_map=interp1(1:Nradial,radial_r_value_flux,radial_XZ_map);

% Bstar_exact_X=(1-q_XZ_map).*BR_XZ_map.*(1-(r_XZ_map/R0).^2);
% Bstar_exact_Z=(1-q_XZ_map).*BZ_XZ_map.*(1-(r_XZ_map/R0).^2);
% Bstar_exact_X=(1-q_XZ_map).*BR_XZ_map./sqrt(1+(r_XZ_map/Raxis).^2);
% Bstar_exact_Z=(1-q_XZ_map).*BZ_XZ_map./sqrt(1+(r_XZ_map/Raxis).^2);

Bstar_exact_X=(1-q_XZ_map).*BR_XZ_map;
Bstar_exact_Z=(1-q_XZ_map).*BZ_XZ_map;
Bstar_tot_exact=sqrt(Bstar_exact_X.^2+Bstar_exact_Z.^2);
Bstar_tot_recalc=sqrt(BstarX_XZ_map_ini.^2+BstarZ_XZ_map_ini.^2);


figure(8);
% imagesc((Bstar_PR_map_exact-Bstar_PR_map)',[-0.00001 0.00001])
imagesc(X_scale,Z_scale,abs(Bstar_tot_recalc-Bstar_tot_exact)',[-0.00005 0.00005])
xlim([-0.35 0.6])
ylim([-0.65 0.65])


% BH maps from Bstar maps
BHpolX_XZ_map_recalc=BR_XZ_map-BstarX_XZ_map_ini;
BHpolZ_XZ_map_recalc=BZ_XZ_map-BstarZ_XZ_map_ini;
BHpol_X_PR_map=BX_PR_map-BstarX_PR_map;
BHpol_Z_PR_map=BZ_PR_map-BstarZ_PR_map;

% BHpol_X_PR_map=BHpolX_PR_map_recalc;
% BHpol_Z_PR_map=BHpolZ_PR_map_recalc;
BHpol_PR_map=sqrt(BHpol_X_PR_map.^2+BHpol_Z_PR_map.^2);
BHpol_X_map=BHpolX_XZ_map_recalc.*mask_XZ_map;
BHpol_Z_map=BHpolZ_XZ_map_recalc.*mask_XZ_map;

% *************************************************************
% calculating metric coefficients in (P,R) maps
% *************************************************************

%
calculate_ki_angle_geometry;
% calculate_gradients_r_theta;

%clear psi_star_Nradial
%psi_star_Nradial=psi_star_Nradial_avg;
%
figure(3)
subplot(2,1,1)
plot(radial_r_value_flux,q_profile);
xlim([0 0.9*a]);
grid on
ylabel('q');

subplot(2,1,2)
plot(radial_r_value_flux,psi_star_Nradial_avg);
xlim([0 0.9*a]);
grid on
ylabel('\Psi _*');


figure(10);
imagesc((theta_XZ_map)');
axis xy;
colorbar;


calculate_density_temp_profiles

