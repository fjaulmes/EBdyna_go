close all;
load('../../data_tokamak/q_profile.mat', 'psi_rank_q1')

psi_mix=size_r-4
delta_psi_q1=20
psi_core=round(psi_rank_q1-delta_psi_q1)

% load('map_dimensions.mat')
% load('geomtry_bins.mat')
% load('speed_bins.mat')
% load('XZsmall_fields_tokamak.mat')
% load('physics_constants.mat')
% 
VARNAME='Delta_pphi_trapped';
SAVENAME=strcat(VARNAME,'.mat')

load('initialG_2800keV_flat_pre_collapse.mat');
INPUTFILE='flat2800keV_collapse_Glisa_fc1h2_G220713.mat'
alphas_pphi0=(mHe*(X0+R0).*(vpll0)-ZHe*eV*(psi_value0))/eV;
alphas_lambda0=Bavg*alphas_mm0./Ekin0;
PPHI_CUTOFF=mean(alphas_pphi0(find(ALL_TRAPPED_POP.*(~alphas_ejected))))-0.1
ENERGY_POPULATION_0=find((~alphas_ejected).*(alphas_pphi0<=PPHI_CUTOFF).*ALL_TRAPPED_POP);

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
alphas_pos_psi=interp1(psi_scale,1:257,psi_value0);
theta_value0=interp2(scale_X,scale_Z,theta_XZsmall_map',alphas_pos_x,alphas_pos_z,'*nearest');
p_phi0=mHe*(X0+R0).*(vpll0)-ZHe*eV*(psi_value0);

radial_pos=(a/257)*interp2(scale_X,scale_Z,radial_XZsmall_map',alphas_pos_x,alphas_pos_z,'*linear');
alphas_kappa=sqrt((alphas_Ekin*R0+Bavg*alphas_mm.*(radial_pos-R0))./(2*alphas_mm.*radial_pos*Bavg));

EKIN_INF=2000
EKIN_SUP=2000


% alphas_pos_psi=interp2(scale_X,scale_Z,psi_norm_XZsmall_map',alphas_pos_x,alphas_pos_z,'*linear');
psi_pos0=alphas_pos_psi;

clear LOW_ENERGY_POPULATION

%load('alphas_ejected_count.mat');

% alphas_lambda0=alphas_mm0./Ekin0;
% ENERGY_POPULATION_0=find(alphas_Ekin>8*1e3);
load(INPUTFILE);
% load('flat8keV_collapse_Glisa_fc1h2_G030713.mat');
% ENERGY_POPULATION_0=find((~alphas_ejected).*(alphas_kappa>=1).*(alphas_kappa<=10));
% ENERGY_POPULATION_0=find((~alphas_ejected).*(alphas_kappa<=1).*(psi_pos0>=psi_mix).*(vpll0<0));
% ENERGY_POPULATION_0=find((~alphas_ejected).*(psi_pos0>psi_rank_q1).*(alphas_lambda0<1).*(vpll0>0));
% ENERGY_POPULATION=find(alphas_Ekin>1*1e3);
ENERGY_POPULATION=ENERGY_POPULATION_0;
post_dist_rescale=size(ENERGY_POPULATION_0,1)/size(ENERGY_POPULATION,1);

disp('EVALUATED_POPULATION')
disp(size(ENERGY_POPULATION,1));
%LOW_ENERGY_POPULATION=find(alphas_Ekin(LOW_ENERGY_POPULATION)<120*1e3);

%load('initial_alphas_pre_collapse');
%load('initial_alphas_pre_collapse');

% time=2e-5;
L2=['after ' num2str(round(time*1e6)) ' \mus']

alphas_pos_psi_recalc=interp1(psi_scale,1:257,alphas_psi_value_corr);

alphas_pos_psi=interp2(scale_X,scale_Z,psi_norm_XZsmall_map',alphas_pos_x,alphas_pos_z,'*linear');
% alphas_value_psi=interp2(scale_X,scale_Z,psi_XZsmall_map',alphas_pos_x,alphas_pos_z,'*linear');
theta_value=interp2(scale_X,scale_Z,theta_XZsmall_map',alphas_pos_x,alphas_pos_z,'*nearest');

figure(1)
set(gca,'FontSize',22);

hold on
grid on;

X_bin_size=0.06;
X_range=(-0.9:X_bin_size:0.9);


NXpos0=histc(X0(ENERGY_POPULATION_0),X_range);xlim([-0.9 0.9]);
NXpos=post_dist_rescale*histc(alphas_pos_x(ENERGY_POPULATION),X_range);xlim([-0.9 0.9]);

plot(X_range,NXpos0,'b');
plot(X_range,NXpos,'r');

% h = findobj(gca,'Type','patch');
% set(h(2),'FaceColor','b','EdgeColor','k');
% set(h(1),'FaceColor','r');

legend('initial',L2)
title('axial position (X)');


figure(2)
set(gca,'FontSize',22);

hold on;
grid on;

vpll_bin_size=0.05*1e7;
vpll_range=(-1.6:0.05:1.6)*1e7;
vpll_range=vpll_range';

Nvpll0=histc(vpll0(ENERGY_POPULATION_0),vpll_range-0.5*vpll_bin_size);
Nvpll=post_dist_rescale*histc(alphas_vpll(ENERGY_POPULATION),vpll_range-0.5*vpll_bin_size);

% h = findobj(gca,'Type','patch');
% set(h(2),'FaceColor','b','EdgeColor','k');
% set(h(1),'FaceColor','r');

plot(vpll_range,Nvpll0,'b');
plot(vpll_range,Nvpll,'r');

legend('initial',L2);
title('v_{||}');


figure(3)

radial_bin_size=9;
psi_bin_pos=(1:radial_bin_size:254);


Npsi0=histc(psi_pos0(ENERGY_POPULATION_0),psi_bin_pos-0.5*radial_bin_size);
Npsi=post_dist_rescale*histc(alphas_psi(ENERGY_POPULATION),psi_bin_pos-0.5*radial_bin_size);

set(gca,'FontSize',22);
hold on;
grid on;

plot(psi_bin_pos,Npsi0,'b--','LineWidth',2);
plot(psi_bin_pos,Npsi,'r','LineWidth',2);
xlim([1 256])
% plot([psi_rank_q1 psi_rank_q1],[0 max(Npsi0)],'k--','LineWidth',3)
plot([psi_rank_q1 psi_rank_q1],[0 5000],'k--','LineWidth',3)

legend('initial',L2);
title('flux surface position (\psi)');



figure(4)
set(gca,'FontSize',22);
hold on;
grid on;

lambda_bin_size=0.05;
lambda_bins=(0:lambda_bin_size:1.2);
alphas_lambda0=Bavg*alphas_mm0./Ekin0;
alphas_lambda=Bavg*alphas_mm./alphas_Ekin;
TRAPPED=(alphas_lambda0>=0.99);

Nlambda0=histc(alphas_lambda0(ENERGY_POPULATION_0),lambda_bins-0.5*lambda_bin_size);
Nlambda=post_dist_rescale*histc(alphas_lambda(ENERGY_POPULATION),lambda_bins-0.5*lambda_bin_size);

plot(lambda_bins,Nlambda0,'b');
plot(lambda_bins,Nlambda,'r');


legend('initial',L2);
title('pitch angle evolution');


% alphas_vpll0_map=griddata(psi_pos0,X0,abs(vpll0),PSI_xx,XX_psi);
% alphas_vpll0_map(isnan(alphas_vpll0_map))=0;
%alphas_vpll0_map=alphas_vpll0_map.*volume_bin;

% alphas_vpll_map=griddata(alphas_pos_psi,alphas_pos_x,abs(alphas_vpll),PSI_xx,XX_psi);
% alphas_vpll_map(isnan(alphas_vpll_map))=0;
%alphas_vpll_map=alphas_vpll_map.*volume_bin;

% figure(4)
% grid on; 
% hold on;
% plot(scale_X(X_bin_pos),mean(alphas_vpll0_map,1),'b');
% plot(scale_X(X_bin_pos),mean(alphas_vpll_map,1),'r');
% legend('initial','after mixing')
% title('v_{||}');
% xlim([-0.6 0.6]);

% figure(3)
% imagesc(scale_X,scale_Z,alphas_vpll0_map);
% axis xy
% figure(4)
% imagesc(scale_X,scale_Z,alphas_vpll_map);
% axis xy


% alphas_Ekin0_map=griddata(psi_pos0,X0,Ekin0,PSI_xx,XX_psi);
% alphas_Ekin0_map(isnan(alphas_Ekin0_map))=0;
%alphas_Ekin0_map=alphas_Ekin0_map.*volume_bin;

% alphas_Ekin_map=griddata(alphas_pos_psi,alphas_pos_x,alphas_Ekin,PSI_xx,XX_psi);
% alphas_Ekin_map(isnan(alphas_Ekin_map))=0;
%alphas_Ekin_map=alphas_Ekin_map.*volume_bin;
% 
% figure(5)
% grid on; 
% hold on;
% plot(scale_X(X_bin_pos),mean(alphas_Ekin0_map,1),'b');
% plot(scale_X(X_bin_pos),mean(alphas_Ekin_map,1),'r');
% legend('initial','after mixing')
% title('Ekin');
% xlim([-0.6 0.6]);

% figure(6)
% imagesc(scale_X,scale_Z,alphas_Ekin0_map);
% axis xy
% figure(7)
% imagesc(scale_X,scale_Z,alphas_Ekin_map);
% axis xy

figure(5)
set(gca,'FontSize',22);
hold on;
grid on;

p_phi_scale=(-2:0.25:4);
p_phi_bin_size=0.25;

alphas_pphi0=(mHe*(X0+R0).*(vpll0)-ZHe*eV*(psi_value0))/eV;
alphas_pphi=(mHe*(alphas_pos_x+R0).*(alphas_vpll)-ZHe*eV*(alphas_psi_value_corr))/eV;

p_phi0=mHe*(X0(ENERGY_POPULATION)+R0).*(vpll0(ENERGY_POPULATION))-ZHe*eV*(psi_value0(ENERGY_POPULATION));
p_phi=mHe*(alphas_pos_x(ENERGY_POPULATION)+R0).*(alphas_vpll(ENERGY_POPULATION))-ZHe*eV*(alphas_psi_value_corr(ENERGY_POPULATION));

REDIST_TRAPPED=find((alphas_pphi.*TRAPPED.*(Ekin0>16e4)./(alphas_pphi0.*TRAPPED.*(Ekin0>16e4)))<0.9);
disp('number of resonnant particles = ');
disp(length(REDIST_TRAPPED))

radial_pos_ini=(a/257)*interp2(scale_X,scale_Z,radial_XZsmall_map',X0,Z0,'*linear');
radial_pos_end=(a/257)*interp2(scale_X,scale_Z,radial_XZsmall_map',alphas_pos_x,alphas_pos_z,'*linear');
alphas_kappa_ini=sqrt((Ekin0*R0+Bavg*alphas_mm0.*(radial_pos_ini-R0))./(2*alphas_mm0.*radial_pos_ini*Bavg));
alphas_kappa_end=sqrt((alphas_Ekin*R0+Bavg*alphas_mm.*(radial_pos_end-R0))./(2*alphas_mm.*radial_pos_end*Bavg));

Npphi0=histc(p_phi0/eV,p_phi_scale-0.5*p_phi_bin_size);
Npphi=histc(p_phi/eV,p_phi_scale-0.5*p_phi_bin_size);

% hist(p_phi-p_phi0,40)

plot(p_phi_scale,Npphi0,'b','LineWidth',2);
plot(p_phi_scale,Npphi,'r','LineWidth',2);


legend('initial',L2);
title('canonical parallel angular momentum');


Dp_phi_scale=(-0.8:0.04:0.8);


figure(6);
set(gca,'FontSize',22);
hold on;
grid on;

Delta_pphi=(p_phi-p_phi0)/eV;
assignin('base',VARNAME,Delta_pphi);
save(SAVENAME,VARNAME);

NDpphi=histc(Delta_pphi,Dp_phi_scale-0.02)
plot(Dp_phi_scale,NDpphi,'r','LineWidth',2)
xlabel('\Delta p_\phi')

EVALUATED_POPULATION=length(Delta_pphi)

PPHI_IN_LIM=0.15
PPHI_OUT_LIM=-0.15

REDIST_INPULSED=length(find(Delta_pphi>PPHI_IN_LIM))
REDIST_EXPULSED=length(find(Delta_pphi<PPHI_OUT_LIM))
INPULSED_PRECENTAGE=100*REDIST_INPULSED/EVALUATED_POPULATION
EXPULSED_PRECENTAGE=100*REDIST_EXPULSED/EVALUATED_POPULATION

MEAN_DPPHI_IN=mean(Delta_pphi(REDIST_INPULSED))
MEAN_DPPHI_OUT=mean(Delta_pphi(REDIST_EXPULSED))

AVG_DPPHI_IN=0.01*MEAN_DPPHI_IN*INPULSED_PRECENTAGE
AVG_DPPHI_OUT=0.01*MEAN_DPPHI_OUT*EXPULSED_PRECENTAGE

title(strcat('pphi evol (Ekin=',num2str(EKIN_INF),'keV)'));





figure(7)

subplot(2,1,1);

mHist = hist2d ([Z0, X0], -0.7:0.1:0.7, -0.7:0.1:0.7);
plot2Dhist(mHist, -0.7:0.1:0.7, -0.7:0.1:0.7,  'X', 'Z', 'Nb part ini'); 

subplot(2,1,2);

mHist = hist2d ([alphas_pos_z, alphas_pos_x], -0.7:0.1:0.7, -0.7:0.1:0.7);
plot2Dhist(mHist, -0.7:0.1:0.7, -0.7:0.1:0.7,  'X', 'Z', 'Nb part'); 








figure(8)
set(gca,'FontSize',22);
hold on;
grid on;

theta_bin_size=pi/10;
theta_bin_pos=(theta_bin_size:theta_bin_size:21*theta_bin_size)-0.5*theta_bin_size;


Ntheta0=histc(theta_value0(ENERGY_POPULATION_0),theta_bin_pos-0.5*theta_bin_size);
Ntheta=post_dist_rescale*histc(theta_value(ENERGY_POPULATION),theta_bin_pos-0.5*theta_bin_size);


plot(theta_bin_pos(1:end-1),Ntheta0(1:end-1),'b');
plot(theta_bin_pos(1:end-1),Ntheta(1:end-1),'r');
xlim([0 2*pi]);


legend('initial',L2);
title('poloidal angle (\theta)');
ylim([0 max(Ntheta0)+1000])
