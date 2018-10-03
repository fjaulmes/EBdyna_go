
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

% please update the polynomials to the actual FINESSE.inp values here
%pol_F2=[ 0 -0.14   -6.4*1e-3   0.19500 -0.11370 ];
%pol_P=[   1.00    -0.935700000  1.7319000  -4.602800   2.80980  ];

%pol_F2=[0.0  -0.132   -0.006405659016252    0.222   -0.122 -0.001]
%pol_P =[1.0     -0.65    1.212766253634022   -4.330277244534635    2.801249230231029]
%F2poly_coefs=pol_F2(end:-1:1);
%Ppoly_coefs=pol_P(end:-1:1);
%disp('Updating polynomials....')
%save('../METIS_finesse_data.mat','-append','F2poly_coefs','Ppoly_coefs')


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

%%
pause(0.4);

ASDEX_LIKE_EQUILIBRIUM=0;

psi_data=reshape((psi_PR_map(:,:)),NP*Nradial,1);
psi_XZ_map_linear=griddata(finesse_data_X,finesse_data_Z,psi_data,XX,ZZ,'linear');
psi_XZ_map_linear(isnan(psi_XZ_map_linear)) = SIGN_CO_CURRENT_FIELD*1; 
psi_XZ_map_linear=psi_XZ_map_linear';

%if ASDEX_LIKE_EQUILIBRIUM==0
	psi_data=reshape((psi_PR_map(:,:)),NP*Nradial,1);
	psi_XZ_map_cubic=griddata(finesse_data_X,finesse_data_Z,psi_data,XX,ZZ,'cubic');
	psi_XZ_map_cubic(isnan(psi_XZ_map_cubic)) = 1; 
	psi_XZ_map_cubic=psi_XZ_map_cubic';
%end
%%

disp('Interpolate XZ maps from PR data .... ');
interpolate_iso_XYZ_maps;

close all

%%


pos_axis=round(X_axis/DX)+mid_X;
Rpos_axis=Rpos(pos_axis);
Raxis_pos=pos_axis;


Btor=Bphi_XZ_map;


close all
create_tokamak_map_dimensions;

% *************************************************************
% defining curves approximating the shape of the flux surfaces
% *************************************************************

disp('find flux surface horizontal boundaries...');
if SIGN_CO_CURRENT_FIELD>0 
% 	psi2D=(psi_XZ_map.*(psi_XZ_map<=0));
    psi2D=(0.5*(psi_XZ_map+psi_XZ_map_linear).*(psi_XZ_map<=0));
	psi2D=psi2D-min(min(psi2D))-psi_global;
else
	%psi2D=(0.5*(psi_XZ_map+psi_XZ_map_cubic).*(psi_XZ_map<=0));
	%psi2D=0.5*(psi_XZ_map+psi_XZ_map_cubic).*(psi_XZ_map>=0);
	psi2D=psi_XZ_map.*(psi_XZ_map>=0);
	psi2D=psi2D-max(max(psi2D))+psi_global;
end
% psi2D=psi_XZ_map_cubic;
% psi2D=psi_XZ_map_linear;
find_flux_surface_horizontal_boundaries;

xi_psi_Nradial=zeros(1,Nradial);
xi_psi_Nradial(1)=X_axis;
disp('find displacement profile...');
fitting_displacement_profile

index=0;

figure(2);
% Cont_psi2D=contour(X_scale,Z_scale,psi2D',psi_scale(1:end));

% to be modified using psi2D if the linear map is not giving good enough results
if ASDEX_LIKE_EQUILIBRIUM==1
    Cont_psi2D_lin=contour(X_scale,Z_scale,psi_XZ_map_linear',psi_scale(1:end));
else
    %Cont_psi2D_lin=Cont_psi2D;
    Cont_psi2D_lin=contour(X_scale,Z_scale,psi_XZ_map_linear',psi_scale(1:end));
end
Cont_psi2D=Cont_psi2D_lin;

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
if LARGEST_CONTOUR>1
	X_contour_raw_prev_last=X_contour_raw(LARGEST_CONTOUR-1,:);
	Z_contour_raw_prev_last=Z_contour_raw(LARGEST_CONTOUR-1,:);
else
	X_contour_raw_prev_last=X_contour_raw(LARGEST_CONTOUR+1,:);
	Z_contour_raw_prev_last=Z_contour_raw(LARGEST_CONTOUR+1,:);
end

axis xy equal
xlabel('x')
ylabel('z')
title('psi2D');
hold on;


disp('Building Z curves of flux surface shapes ...');

%if ASDEX_LIKE_EQUILIBRIUM==1
%Cont_psi2D=Cont_psi2D_lin;
%end
extract_contour_information;

fitting_Zpsi_up_profiles;
fitting_Zpsi_down_profiles;
% further improvement : take Ncontours-1 values to be average between
% Ncontours-2 and Ncontours (to reflect innacuracy done before)


%%

%first get mixing radius value
NOQ1SURF=0;
if isempty(find(q_profile<1))
	NOQ1SURF=1;
	disp('No sawtooth possible since no q=1 surface is present!')
    psi_mix_rank=1
	psi_q1=psi_scale(1);
	psi_mix_value=psi_scale(1);
	r_value_q1_mean=0;
end

map_BH;

Dpsi=0.5*(psi_scale(2)-psi_scale(1));
psi_scale_lin=(psi_scale(1):Dpsi:psi_scale(end)+10*Dpsi);

psiH_Nradial=zeros(1,Nradial);
psi_star_Nradial=zeros(NP,Nradial);
psi_star_Nradial_lin=zeros(NP,length(psi_scale_lin));
integ_q_profile_lin=zeros(NP,length(psi_scale_lin));


for (p=1:NP)
    integ_q_profile_lin(p,:)=interp1(psi_scale,(1-q_profile),psi_scale_lin);
end


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

% psiH from psi star
psiH_Nradial=psi_scale-psi_star_Nradial_avg;

% mixing radius position
q1=ones(1,Nradial);
[eps Nq1]=min(abs(q1-q_profile));
psi_rank_q1=Nq1;
if NOQ1SURF==0
	psi_mix_rank=interp1(psi_star_Nradial_avg(psi_rank_q1+1:end),(psi_rank_q1+1:Nradial),0)
	psi_mix_value=interp1(psi_star_Nradial_avg(psi_rank_q1+1:end),psi_scale(psi_rank_q1+1:Nradial),0)
else
	psi_mix_rank=1
	psi_q1=psi_scale(1);
	psi_mix_value=psi_scale(1);
	r_value_q1_mean=0;
end




size_r=round(psi_mix_rank+4)
% build the radial coordinate in meters

evaluate_volume_flux_surfaces_Nradial;

surf_flux_diff=surf_flux(2:end)-surf_flux(1:end-1);
corr_surf=sum(dSurf,1)./surf_flux_diff;
[ mincorr mincorr_pos ]=min(corr_surf);
for r=1:mincorr_pos
	dSurf(:,r)=dSurf(:,r)/mincorr;
end
for r=mincorr_pos+1:Nradial-1
	dSurf(:,r)=dSurf(:,r)/corr_surf(r);
end
%elementary volumes
dVol=2*pi*Rpos_values.*dSurf;

% *************************************************************
% defining the safety factor q=1 contour
% *************************************************************
if NOQ1SURF==0
    build_q1_flux_surface;
end



for(p=1:NP)
    for(r=1:Nradial)
        r_PR_map(p,r)=radial_r_value_flux(r);
        r_data(n)=radial_r_value_flux(r);
    end
end




% verify_GradShavranov_equation_recalc_Bpol;
Btot=sqrt(BR_XZ_map.^2+BZ_XZ_map.^2);





% for information
if NOQ1SURF==0
psiH_mix_value=interp1(psi_star_Nradial_avg(psi_rank_q1+1:end),psiH_Nradial(psi_rank_q1+1:Nradial),0);
Bstar_PR_map_exact=Bpol_PR_map.*(1-q_PR_map);
else
psiH_mix_value=0;
Bstar_PR_map_exact=Bpol_PR_map.*0;
end

psi_star_avg_PR_map=zeros(NP,Nradial);
psiH_PR_map=zeros(NP,Nradial);
psi_star_PR_map=zeros(NP,Nradial);

if NOQ1SURF==0
	for (p=1:NP)
		for (r=2:Nradial)
			psi_star_avg_PR_map(p,r)=psi_star_Nradial_avg(r);
		end
	end
	psi_star_PR_map=psi_star_Nradial;
	% psiH map from psi star map
	psiH_PR_map=psi_PR_map-psi_star_PR_map;

	% good interpolation for psi star XZ map
	gf_data=reshape(psi_star_avg_PR_map(:,:),Nradial*NP,1);
	if (SIGN_TOROIDAL_FIELD>0)
		gf_data=max(gf_data,-0.05);
	else
		gf_data=min(gf_data,0.05);
	end
	if Nq1>1
		psi_star_XZ_map_unblurred=gridfit(finesse_data_X,finesse_data_Z,gf_data,X_scale,Z_scale,'smoothness' ,0.5);
		psi_star_XZ_map_unblurred(isnan(psi_star_XZ_map_unblurred))=0;
		psi_star_XZ_map=psi_star_XZ_map_unblurred';
	else
		psi_star_data=reshape((psi_star_PR_map(:,:)),NP*Nradial,1);
		psi_star_XZ_map=griddata(finesse_data_X,finesse_data_Z,psi_star_data,XX,ZZ,'cubic');
		psi_star_XZ_map(isnan(psi_star_XZ_map))=0;
		psi_star_XZ_map=psi_star_XZ_map';
	end
else
	psi_star_data=reshape((psi_star_PR_map(:,:)),NP*Nradial,1);
	psi_star_XZ_map=griddata(finesse_data_X,finesse_data_Z,psi_star_data,XX,ZZ,'cubic');
	psi_star_XZ_map(isnan(psi_star_XZ_map))=0;
	psi_star_XZ_map=psi_star_XZ_map';
end

% psiH map from psi star map
psiH_XZ_map=psi_XZ_map-psi_star_XZ_map;

calculate_Bstar_inital;
% recalculate_BH;

% the discrepancy between the two fields calculations
% is our error for helical flux estimate

Bstar_exact_X=(1-q_XZ_map).*BR_XZ_map;
Bstar_exact_Z=(1-q_XZ_map).*BZ_XZ_map;
Bstar_tot_exact=sqrt(Bstar_exact_X.^2+Bstar_exact_Z.^2);
Bstar_tot_recalc=sqrt(BstarX_XZ_map_ini.^2+BstarZ_XZ_map_ini.^2);


%figure(8);
%title('Bstar error')
%imagesc(X_scale,Z_scale,abs(Bstar_tot_recalc-Bstar_tot_exact)',[-0.00005 0.00005])
%xlim([-0.7 0.9]*a)
%ylim([-0.9 0.9]*a)


% BH maps from Bstar maps
BHpolX_XZ_map_recalc=BR_XZ_map-BstarX_XZ_map_ini;
BHpolZ_XZ_map_recalc=BZ_XZ_map-BstarZ_XZ_map_ini;
BHpol_X_PR_map=BX_PR_map-BstarX_PR_map;
BHpol_Z_PR_map=BZ_PR_map-BstarZ_PR_map;

BHpol_PR_map=sqrt(BHpol_X_PR_map.^2+BHpol_Z_PR_map.^2);
BHpol_X_map=BHpolX_XZ_map_recalc.*mask_XZ_map;
BHpol_Z_map=BHpolZ_XZ_map_recalc.*mask_XZ_map;

% *************************************************************
% calculating metric coefficients in (P,R) maps
% *************************************************************

%
calculate_ki_angle_geometry;

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


METIS_density_temp_profiles;

volume_flux=volume_profile;
disp('Volume of vaccuum vessel = ');
disp(volume_flux(end));

% 
P_prime=gradient(P_initial_profile,psi_scale);
P_prime(end-1)=P_prime(end-2);
P_prime(end)=P_prime(end-1);
for(p=1:NP)
    for(r=1:Nradial)
        Pprime_PR_map(p,r)=P_prime(r);
    end
end


%----------------------------------
% maps for toroidal current density calculations
%----------------------------------
SIGN_CURRENT=SIGN_CO_CURRENT_FIELD*SIGN_TOROIDAL_FIELD
J_PHI_profile=(-P_prime.*R0^2-0.5*F2_prime/mu0)/R0;
J_PHI_PR_map=-Rpos_PR_map.*Pprime_PR_map-(FFprime_PR_map/mu0)./Rpos_PR_map;
J_PHI_flat_pressure_PR_map=-FFprime_PR_map/mu0/R0;
J_PHI_flat_profile=(-0.5*F2_prime/mu0)/R0;

J_PHI_profile=SIGN_CURRENT*abs(J_PHI_profile);
J_PHI_PR_map=SIGN_CURRENT*abs(J_PHI_PR_map);
J_PHI_flat_pressure_PR_map=SIGN_CURRENT*abs(J_PHI_flat_pressure_PR_map);
J_PHI_flat_profile=SIGN_CURRENT*abs(J_PHI_flat_profile);

JPHI_XZ_map=zeros(NZ,NZ);
JPHI_XZ_map=griddata(finesse_data_X,finesse_data_Z,reshape(J_PHI_PR_map,1,NP*Nradial),XX,ZZ,'cubic');
JPHI_XZ_map(find(isnan(JPHI_XZ_map)))=0;
JPHI_XZ_map=JPHI_XZ_map';
JPHI_flat_pressure_XZ_map=zeros(NZ,NZ);
JPHI_flat_pressure_XZ_map=griddata(finesse_data_X,finesse_data_Z,reshape(J_PHI_flat_pressure_PR_map,1,NP*Nradial),XX,ZZ,'cubic');
JPHI_flat_pressure_XZ_map=JPHI_flat_pressure_XZ_map';

JPHI_XZ_map(find(isnan(JPHI_XZ_map)))=0;


I_flux_map=Z_psi_fit_up*0;
I_flux=zeros(1,Nradial);
I_flat_flux=zeros(1,Nradial);



%DS=(DX)^2;

for (n=2:Nradial)
    for(x=X1_Nradial(n):X2_Nradial(n))
        Z_max_up=Z_psi_fit_up(n,x);
        Z_max_down=Z_psi_fit_down(n,x);
        DS=(Z_max_up-Z_max_down)*DX;
        I_flux_map(n,x)=sum(JPHI_XZ_map(x,round(Z_max_down/DX)+mid_Z:round(Z_max_up/DX)+mid_Z))*DX*DX;
        I_flat_flux(n)=I_flat_flux(n)+JPHI_flat_pressure_XZ_map(x,mid_Z)*DS;

    end
end
for (n=2:Nradial)
    I_flux(n)=sum(I_flux_map(n,:));
end


for (n=2:Nradial)
    I_flux_diff(n)=I_flux(n-1)+I_flux(n);
    I_flat_flux_diff(n)=I_flat_flux(n-1)+I_flat_flux(n);
end

Iaxis=I_flux(end);
sign_Iaxis=sign(I_flux(round(0.5*Nradial)));

tor_flux_scale=q_profile*0;
for n=2:Nradial
    tor_flux_scale(n)=trapz(psi_scale(1:n)',q_profile(1:n));
end
rho_tor_scale=sqrt(tor_flux_scale/max(tor_flux_scale));


save_PR_maps
