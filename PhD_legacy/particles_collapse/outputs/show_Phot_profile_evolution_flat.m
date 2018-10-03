% close all;

psi_mix=size_r-4
delta_psi_q1=20;
psi_core=psi_rank_q1-delta_psi_q1

% load('map_dimensions.mat')
% load('geomtry_bins.mat')
% load('speed_bins.mat')
% load('XZsmall_fields_tokamak.mat')
% load('physics_constants.mat')
% 
% load('initial_alphas_distributions.mat');
load('initialG_1600keV_flat_pre_collapse.mat');

disp('Total number of particles simulates');
Nalphas_simulated=length(alphas_pos_x);
disp(Nalphas_simulated);

X0=alphas_pos_x;
Z0=alphas_pos_z;
vpll0=alphas_vpll;
alphas_mm0=alphas_mm;
Bfield0=interp2(scale_X,scale_Z,Btot_XZ_map',alphas_pos_x,alphas_pos_z,'*linear');
alphas_Eperp=alphas_mm0.*Bfield0;
alphas_Epll=0.5*(mHe/eV)*alphas_vpll.^2;
alphas_Ekin=alphas_Eperp+alphas_Epll;
Ekin0=alphas_Ekin;
psi_value0=interp2(scale_X,scale_Z,psi_XZsmall_map',alphas_pos_x,alphas_pos_z,'*linear');
% alphas_pos_psi=interp1(psi_scale,1:257,psi_value0);
alphas_pos_psi=interp2(scale_X,scale_Z,psi_norm_XZsmall_map',alphas_pos_x,alphas_pos_z,'*linear');
theta_value0=interp2(scale_X,scale_Z,theta_XZsmall_map',alphas_pos_x,alphas_pos_z,'*nearest');
alphas_lambda0=Bavg*alphas_mm0./Ekin0;
radial_pos=(a/257)*interp2(scale_X,scale_Z,radial_XZsmall_map',alphas_pos_x,alphas_pos_z,'*linear');
alphas_kappa=sqrt((alphas_Ekin*R0+Bavg*alphas_mm.*(radial_pos-R0))./(2*alphas_mm.*radial_pos*Bavg));

kappa_scale=[0 0.8 1.6 2.4 3.2 ];

% alphas_pos_psi=interp2(scale_X,scale_Z,psi_norm_XZsmall_map',alphas_pos_x,alphas_pos_z,'*linear');
psi_pos0=alphas_pos_psi;

clear ENERGY_POPULATION ENERGY_POPULATION_0

%load('alphas_ejected_count.mat');

% alphas_lambda0=alphas_mm0./Ekin0;
% ENERGY_POPULATION_0=find(alphas_Ekin>8*1e3);
load('flat2800keV_collapse_Glisa_fc1h2_G220713.mat');
alphas_ejected_half=alphas_ejected;
% load('alphas_collapse_Glisa_fc1h2_G280613.mat');
% alphas_ejected=alphas_ejected_half+alphas_ejected;
alphas_ejected(alphas_ejected>1)=1;
EKIN_INF=2
% ENERGY_POPULATION_0=find((~alphas_ejected).*(alphas_lambda0>1).*(Ekin0>240e4).*(Ekin0>160e4));
% ENERGY_POPULATION_0=find((~alphas_ejected).*(alphas_Ekin>1*1e6).*(psi_pos0<psi_mix).*(alphas_lambda0>0.8).*(alphas_lambda0<1.05));
% ENERGY_POPULATION_0=find((~alphas_ejected).*(alphas_Ekin>1*1e6));
% POPULATION_RANGE=(alphas_Ekin>EKIN_INF*1e3).*(alphas_lambda0<1).*(alphas_kappa>1).*(vpll0>0).*(psi_pos0>=psi_core).*(psi_pos0<=psi_rank_q1+delta_psi_q1);
% POPULATION_RANGE=(alphas_Ekin>EKIN_INF*1e3).*(alphas_lambda0>=0.8).*(alphas_kappa<=1).*(vpll0<0).*(psi_pos0>=psi_mix);
% POPULATION_RANGE=(alphas_Ekin>EKIN_INF*1e3).*(alphas_kappa>=1).*(vpll0>0);
POPULATION_RANGE=(alphas_Ekin>EKIN_INF*1e3).*(psi_pos0<=psi_core).*COUNTER_PASSING_POP;
ENERGY_POPULATION_0=find(POPULATION_RANGE);
ENERGY_POPULATION=find((~alphas_ejected).*POPULATION_RANGE);
% ENERGY_POPULATION=find((~alphas_ejected(ENERGY_POPULATION_0)));
Nkappa0=histc(alphas_kappa(ENERGY_POPULATION_0),kappa_scale+0.4);

post_dist_rescale=size(ENERGY_POPULATION_0,1)/size(ENERGY_POPULATION,1);

disp('EVALUATED_POPULATION')
EVALUATED_POPULATION=length(ENERGY_POPULATION);
disp(EVALUATED_POPULATION);
%LOW_ENERGY_POPULATION=find(alphas_Ekin(LOW_ENERGY_POPULATION)<120*1e3);

%load('initial_alphas_pre_collapse');
%load('initial_alphas_pre_collapse');

% time=4e-4;
L2=['after ' num2str(round(time*1e6)) ' \mus']

alphas_pos_psi=interp2(scale_X,scale_Z,psi_norm_XZsmall_map',alphas_pos_x,alphas_pos_z,'*linear');
alphas_value_psi=interp2(scale_X,scale_Z,psi_XZsmall_map',alphas_pos_x,alphas_pos_z,'*linear');
theta_value=interp2(scale_X,scale_Z,theta_XZsmall_map',alphas_pos_x,alphas_pos_z,'*nearest');






figure(3)
set(gca,'FontSize',22);
hold on;
grid on;

radial_bin_size_half=5;
radial_bin_size=radial_bin_size_half*2;
psi_bin_pos=(1:radial_bin_size:255-radial_bin_size);
nb_psi_bins=length(psi_bin_pos)


Npsi0=histc(psi_pos0(ENERGY_POPULATION_0),psi_bin_pos);
Npsi=histc(alphas_psi(ENERGY_POPULATION),psi_bin_pos);

for (psi=1:256)
    volume_radial(psi)=sum(volume_tor_diff(psi,:));
end
volume_radial_bin=zeros(nb_psi_bins,1);
for (psi=1:nb_psi_bins)
    volume_radial_bin(psi)=sum(volume_radial((psi-1)*radial_bin_size+1:psi*radial_bin_size));
end

density_radial_ini=Npsi0./volume_radial_bin;
density_radial_end=Npsi./volume_radial_bin;

plot(psi_bin_pos+0.5*radial_bin_size,density_radial_ini,'b--','LineWidth',2);
plot(psi_bin_pos+0.5*radial_bin_size,density_radial_end,'r','LineWidth',2);


legend('initial',L2);
title('flux surface position (\psi)');




lambda_bin_size=0.05;
lambda_bins=(0:lambda_bin_size:1.2);
alphas_lambda0=Bavg*alphas_mm0./Ekin0;
alphas_lambda=Bavg*alphas_mm./alphas_Ekin;
TRAPPED=(alphas_lambda0>=0.99);

Delta_lambda=alphas_lambda-alphas_lambda0;
Dlambda_scale=(-0.1:0.01:0.1);
lambda_scale=(0.05:0.05:1.25);

Nlambda=histc(Delta_lambda,Dlambda_scale-0.005);
% Nlambda=histc(alphas_lambda(ENERGY_POPULATION),lambda_scale-0.025);
figure(1)

plot(Dlambda_scale,Nlambda,'r','LineWidth',2)



Dpphi_bin_size=0.04;
Dp_phi_scale=(-0.8:Dpphi_bin_size:0.8);
p_phi_bin_size=0.25;
p_phi_scale=(-2:0.25:4);

alphas_pphi0=(mHe*(X0+R0).*(vpll0)-2*eV*(psi_value0))/eV;
alphas_pphi=(mHe*(alphas_pos_x+R0).*(alphas_vpll)-2*eV*(alphas_psi_value_corr))/eV;

p_phi0=mHe*(X0+R0).*(vpll0)-2*eV*(psi_value0);
p_phi=mHe*(alphas_pos_x+R0).*(alphas_vpll)-2*eV*(alphas_psi_value_corr);



figure(6);
set(gca,'FontSize',22);
hold on;
grid on;

hold on;
grid on;

Delta_pphi=(p_phi-p_phi0)/eV;
NDpphi=histc(Delta_pphi(ENERGY_POPULATION),Dp_phi_scale-0.5*Dpphi_bin_size);
plot(Dp_phi_scale,NDpphi,'b','LineWidth',2);
xlabel('\Delta p_\phi')

PPHI_IN_LIM=0.15
PPHI_OUT_LIM=-0.15

POP_FILTER=zeros(Nalphas_simulated,1);
POP_FILTER_0(ENERGY_POPULATION_0)=1;
POP_FILTER(ENERGY_POPULATION)=1;

REDIST_INPULSED=find((Delta_pphi>PPHI_IN_LIM).*POP_FILTER);
REDIST_EXPULSED=find((Delta_pphi<PPHI_OUT_LIM).*POP_FILTER);
INPULSED_PRECENTAGE=100*length(REDIST_INPULSED)/length(ENERGY_POPULATION_0)
EXPULSED_PRECENTAGE=100*length(REDIST_EXPULSED)/length(ENERGY_POPULATION)

MEAN_DPPHI_IN=mean(Delta_pphi(REDIST_INPULSED))
MEAN_DPPHI_OUT=mean(Delta_pphi(REDIST_EXPULSED))

AVG_DPPHI_IN=0.01*MEAN_DPPHI_IN*INPULSED_PRECENTAGE
AVG_DPPHI_OUT=0.01*MEAN_DPPHI_OUT*EXPULSED_PRECENTAGE

AVG_DPPHI=0.5*(AVG_DPPHI_IN+AVG_DPPHI_OUT)

% title(strcat('co-passing (Ekin>',num2str(EKIN_INF),'keV)'));


DELTA_2D_RES=DX*40;
XINF=-0.8;
XSUP=0.8;
ZINF=-0.8;
ZSUP=0.8;


[ posX_scale posZ_scale n_2D_map_ini temp_2D_map_ini ]=build_temperature_2D(X0(ENERGY_POPULATION_0),Z0(ENERGY_POPULATION_0),Ekin0(ENERGY_POPULATION_0)*eV,XINF,XSUP,ZINF,ZSUP,DELTA_2D_RES);
[ posX_scale posZ_scale n_2D_map_end temp_2D_map_end ]=build_temperature_2D(alphas_pos_x(ENERGY_POPULATION),alphas_pos_z(ENERGY_POPULATION),alphas_Ekin(ENERGY_POPULATION)*eV,XINF,XSUP,ZINF,ZSUP,DELTA_2D_RES);

volume_X=2*pi*(R0+posX_scale)*(DELTA_2D_RES)^2;


n_2D_map_ini=n_2D_map_ini*1e11;
n_2D_map_end=n_2D_map_end*1e11;

for (x=1:length(posX_scale))
    n_2D_map_ini(x,:)=n_2D_map_ini(x,:)/volume_X(x);
    n_2D_map_end(x,:)=n_2D_map_end(x,:)/volume_X(x);
end
% 
% figure(7)
% 
% subplot(2,1,1);
% 
% mHist = hist2d ([Z0(ENERGY_POPULATION_0), X0(ENERGY_POPULATION_0)], -0.7:0.1:0.7, -0.7:0.1:0.7);
% mHist=mHist*1e11;
% plot2Dhist(mHist, -0.7:0.1:0.7, -0.7:0.1:0.7,  'X', 'Z', 'Nb part ini'); 
% 
% subplot(2,1,2);
% 
% mHist = hist2d ([alphas_pos_z(ENERGY_POPULATION), alphas_pos_x(ENERGY_POPULATION)], -0.7:0.1:0.7, -0.7:0.1:0.7);
% mHist=mHist*1e11;
% plot2Dhist(mHist, -0.7:0.1:0.7, -0.7:0.1:0.7,  'X', 'Z', 'Nb part'); 



figure(8)

P_SCALE=[0 0.5]*1e5

subplot(2,1,1);
imagesc(posX_scale,posZ_scale,(n_2D_map_ini.*temp_2D_map_ini)',P_SCALE);
axis xy;colorbar;

subplot(2,1,2);
imagesc(posX_scale,posZ_scale,(n_2D_map_end.*temp_2D_map_end)',P_SCALE);
axis xy;colorbar;






