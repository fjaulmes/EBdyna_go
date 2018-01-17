close all;
% Bavg=3.46
psi_mix=size_r-4
delta_psi_q1=20
psi_core=round(psi_rank_q1-delta_psi_q1)

% load('map_dimensions.mat')
% load('geomtry_bins.mat')
% load('speed_bins.mat')
% load('XZsmall_fields_tokamak.mat')
% load('physics_constants.mat')
% 
% load('initial_alphas_distributions.mat');
load('initial_NBI_1MEV_D_pre_collapse_all.mat');

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
alphas_pos_psi=interp2(scale_X,scale_Z,psi_norm_XZsmall_map',alphas_pos_x,alphas_pos_z,'*linear');
% alphas_pos_psi=interp1(psi_scale,1:257,psi_value0);
theta_value0=interp2(scale_X,scale_Z,theta_XZsmall_map',alphas_pos_x,alphas_pos_z,'*nearest');
alphas_lambda0=Bavg*alphas_mm0./Ekin0;
radial_pos=(a/257)*interp2(scale_X,scale_Z,radial_XZsmall_map',alphas_pos_x,alphas_pos_z,'*linear');
alphas_kappa=sqrt((alphas_Ekin*R0+Bavg*alphas_mm.*(radial_pos-R0))./(2*alphas_mm.*radial_pos*Bavg));

kappa_scale=[0 0.8 1.6 2.4 3.2 ];

% alphas_pos_psi=interp2(scale_X,scale_Z,psi_norm_XZsmall_map',alphas_pos_x,alphas_pos_z,'*linear');
alphas_psi0=alphas_pos_psi;

clear LOW_ENERGY_POPULATION

%load('alphas_ejected_count.mat');

% alphas_lambda0=alphas_mm0./Ekin0;
% ENERGY_POPULATION_0=find(alphas_Ekin>8*1e3);
% load('alphas_collapse_Glisa_fc1h2_G260713.mat');
load('NBI_1MEV_fc2h2_all.mat');
% alphas_ejected_half=alphas_ejected;
% load('alphas_collapse_Glisa_fc1h2_G280613.mat');
% alphas_ejected=alphas_ejected_half+alphas_ejected;
alphas_ejected(alphas_ejected>1)=1;
EKIN_INF=200
EKIN_SUP=1600
alphas_pos_psi=interp2(scale_X,scale_Z,psi_norm_XZsmall_map',alphas_pos_x,alphas_pos_z,'*linear');

% ENERGY_POPULATION_0=find((~alphas_ejected).*(alphas_lambda0>1).*(Ekin0>240e4).*(Ekin0>160e4));
% ENERGY_POPULATION_0=find((~alphas_ejected).*(alphas_Ekin>1*1e6).*(psi_pos0<psi_mix).*(alphas_lambda0>0.8).*(alphas_lambda0<1.05));
% ENERGY_POPULATION_0=find((~alphas_ejected).*(alphas_Ekin>1*1e6));
% POPULATION_RANGE=(alphas_Ekin>EKIN_INF*1e3).*(alphas_pos_psi>=psi_core).*(alphas_pos_psi<=psi_rank_q1);
% POPULATION_RANGE=(alphas_Ekin>EKIN_INF*1e3).*(psi_pos0>=psi_core).*(psi_pos0<=psi_rank_q1);
% POPULATION_RANGE=(alphas_Ekin>EKIN_INF*1e3).*(psi_pos0>=psi_core).*(psi_pos0<=psi_rank_q1+delta_psi_q1);
% POPULATION_RANGE=(alphas_Ekin>EKIN_INF*1e3).*(alphas_pos_psi>=psi_core).*(alphas_pos_psi<=psi_rank_q1+delta_psi_q1).*(vpll0>0);
POPULATION_RANGE=(alphas_Ekin>EKIN_INF*1e3).*(alphas_Ekin<EKIN_SUP*1e3).*(alphas_psi0>=psi_core).*(alphas_psi0<=psi_outer).*(COUNTER_PASSING_POP);
ENERGY_POPULATION_0=find(POPULATION_RANGE);
ENERGY_POPULATION=find((~alphas_ejected).*POPULATION_RANGE);

Nkappa0=histc(alphas_kappa(ENERGY_POPULATION_0),kappa_scale+0.4);

post_dist_rescale=size(ENERGY_POPULATION_0,1)/size(ENERGY_POPULATION,1);

disp('EVALUATED_POPULATION')
EVALUATED_POPULATION=length(ENERGY_POPULATION);
disp(EVALUATED_POPULATION);
%LOW_ENERGY_POPULATION=find(alphas_Ekin(LOW_ENERGY_POPULATION)<120*1e3);

%load('initial_alphas_pre_collapse');
%load('initial_alphas_pre_collapse');

% time=2e-5;
L2=['after ' num2str(round(time*1e6)) ' \mus']

alphas_value_psi=interp2(scale_X,scale_Z,psi_XZsmall_map',alphas_pos_x,alphas_pos_z,'*linear');
theta_value=interp2(scale_X,scale_Z,theta_XZsmall_map',alphas_pos_x,alphas_pos_z,'*nearest');

figure(1)
set(gca,'FontSize',22);

hold on
grid on;

NXpos0=histc(X0(ENERGY_POPULATION_0),(-15:15)*0.06);xlim([-0.7 0.7]);
NXpos=post_dist_rescale*histc(alphas_pos_x(ENERGY_POPULATION),(-15:15)*0.06);xlim([-0.7 0.7]);

plot((-15:15)*0.06,NXpos0,'b');
plot((-15:15)*0.06,NXpos,'r');

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

AVG_NEG_SPEED_INI=sum(Nvpll0(1:32).*vpll_range(1:32))/sum(Nvpll0(1:32))
AVG_POS_SPEED_INI=sum(Nvpll0(34:end).*vpll_range(34:end))/sum(Nvpll0(34:end))
AVG_NEG_SPEED=sum(Nvpll(1:32).*vpll_range(1:32))/sum(Nvpll(1:32))
AVG_POS_SPEED=sum(Nvpll(34:end).*vpll_range(34:end))/sum(Nvpll(34:end))


% h = findobj(gca,'Type','patch');
% set(h(2),'FaceColor','b','EdgeColor','k');
% set(h(1),'FaceColor','r');

plot(vpll_range,Nvpll0,'b--','LineWidth',2);
plot(vpll_range,Nvpll,'r','LineWidth',2);

legend('initial @q\sim1',L2);
title('v_{||}');


figure(3)
set(gca,'FontSize',22);
hold on;
grid on;

radial_bin_size=8;
psi_bin_pos=(1:radial_bin_size:254);


Npsi0=histc(alphas_psi0(ENERGY_POPULATION_0),psi_bin_pos-0.5*radial_bin_size);
Npsi=histc(alphas_psi(ENERGY_POPULATION),psi_bin_pos-0.5*radial_bin_size);


plot(psi_bin_pos,Npsi0,'b--','LineWidth',2);
plot(psi_bin_pos,Npsi,'r','LineWidth',2);

plot([psi_rank_q1 psi_rank_q1],[0 5000],'k--','LineWidth',3)

legend('initial @q\sim1',L2);
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



figure(5)
set(gca,'FontSize',22);
hold on;
grid on;

Dp_phi_scale=(-0.7:0.04:0.7);
p_phi_scale=(-2:0.25:4);
p_phi_bin_size=0.25;

alphas_pphi0=(mHe*(X0+R0).*(vpll0)-2*eV*(psi_value0))/eV;
alphas_pphi=(mHe*(alphas_pos_x+R0).*(alphas_vpll)-2*eV*(alphas_psi_value_corr))/eV;

p_phi0=mHe*(X0+R0).*(vpll0)-2*eV*(psi_value0);
p_phi=mHe*(alphas_pos_x+R0).*(alphas_vpll)-2*eV*(alphas_psi_value_corr);

REDIST_TRAPPED=find((alphas_pphi.*TRAPPED.*(Ekin0>16e4)-alphas_pphi0.*TRAPPED.*(Ekin0>16e4))<-0.12);

radial_pos_ini=(a/257)*interp2(scale_X,scale_Z,radial_XZsmall_map',X0,Z0,'*linear');
radial_pos_end=(a/257)*interp2(scale_X,scale_Z,radial_XZsmall_map',alphas_pos_x,alphas_pos_z,'*linear');
alphas_kappa_ini=sqrt((Ekin0*R0+Bavg*alphas_mm0.*(radial_pos_ini-R0))./(2*alphas_mm0.*radial_pos_ini*Bavg));
alphas_kappa_end=sqrt((alphas_Ekin*R0+Bavg*alphas_mm.*(radial_pos_end-R0))./(2*alphas_mm.*radial_pos_end*Bavg));

Npphi0=histc(p_phi0(ENERGY_POPULATION)/eV,p_phi_scale-0.5*p_phi_bin_size);
Npphi=histc(p_phi(ENERGY_POPULATION)/eV,p_phi_scale-0.5*p_phi_bin_size);


plot(p_phi_scale,Npphi,'r','LineWidth',2);
plot(p_phi_scale,Npphi0,'b','LineWidth',2);


legend('initial',L2);
title('canonical parallel angular momentum');

figure(6);
set(gca,'FontSize',26);
hold on;
grid on;

hold on;
grid on;

Delta_pphi=(p_phi-p_phi0)/eV;
NDpphi=histc(Delta_pphi(ENERGY_POPULATION),Dp_phi_scale-0.02);
plot(Dp_phi_scale,NDpphi/max(NDpphi),'b--','LineWidth',3)
plot(Dp_phi_scale,NDpphi_co/max(NDpphi_co),'r','LineWidth',3)
xlabel('\Delta p_\phi')

PPHI_IN_LIM=0.16
PPHI_OUT_LIM=-0.16

REDIST_INPULSED=find(Delta_pphi>PPHI_IN_LIM);
REDIST_EXPULSED=find(Delta_pphi<PPHI_OUT_LIM);
INPULSED_PRECENTAGE=100*length(REDIST_INPULSED)/EVALUATED_POPULATION
EXPULSED_PRECENTAGE=100*length(REDIST_EXPULSED)/EVALUATED_POPULATION

MEAN_DPPHI_IN=mean(Delta_pphi(REDIST_INPULSED))
MEAN_DPPHI_OUT=mean(Delta_pphi(REDIST_EXPULSED))

AVG_DPPHI_IN=0.01*MEAN_DPPHI_IN*INPULSED_PRECENTAGE
AVG_DPPHI_OUT=0.01*MEAN_DPPHI_OUT*EXPULSED_PRECENTAGE

title(strcat('co-passing (Ekin>',num2str(EKIN_INF),'keV)'));

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
