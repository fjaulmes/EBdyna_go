reset_data_analysis_environment;

rho_scale=radial_r_value_flux/max(radial_r_value_flux);

psi_mix=size_r-4
% delta_psi_q1=20;
% psi_core=psi_rank_q1-delta_psi_q1

CRITICAL_ENERGY=400*1e3
CRITICAL_ENERGY_PASSING=1500*1e3;
CRITICAL_ENERGY_PASSING=CRITICAL_ENERGY

radial_bin_size_half=8;
radial_bin_size=radial_bin_size_half*2;
psi_bin_pos=(1:radial_bin_size:255-radial_bin_size);
nb_psi_bins=length(psi_bin_pos)
for (psi=1:256)
    volume_radial(psi)=sum(volume_tor_diff(psi,:));
end
volume_radial_bin=zeros(nb_psi_bins,1);
for (psi=1:nb_psi_bins)
    volume_radial_bin(psi)=sum(volume_radial((psi-1)*radial_bin_size+1:psi*radial_bin_size));
end
load('initialG_alphas_vA_all_pre_collapse.mat')
load('initialG_alphas_vA_all_precession_stats.mat')
EKIN_INF=1;
EKIN_SUP=CRITICAL_ENERGY;
EKIN_POP=find(alphas_Ekin<=EKIN_SUP);

Npsi0=histc(alphas_psi(EKIN_POP),psi_bin_pos);
% Npsi=histc(alphas_psi(EKIN_POP),psi_bin_pos);

load('alphas_vA_all_collapse_Glisa_fc0p8h2_G160414.mat')
Npsi=histc(alphas_psi(EKIN_POP),psi_bin_pos);


density_radial_ini=Npsi0./volume_radial_bin;
density_radial_end=Npsi./volume_radial_bin;

figure(3)

subplot(3,1,1)
set(gca,'FontSize',20);
hold on;
grid on;
title('low energies (E_{kin} < E_{crit})')
plot(rho_scale(round(psi_bin_pos+0.5*radial_bin_size)),density_radial_ini/max(density_radial_ini),'b--','LineWidth',2);
plot(rho_scale(round(psi_bin_pos+0.5*radial_bin_size)),density_radial_end/max(density_radial_ini),'r','LineWidth',2);
% xlabel('r (m)')
% ylabel('n/max(n)')

xlim([0 0.6])
h_legend=legend('initial','post collapse');
set(h_legend,'FontSize',16);


 load('initialG_alphas_vA_all_precession_stats.mat')
 load('initialG_alphas_vA_all_pre_collapse.mat');
% load('initialG_alphas_lowEkin_pre_collapse.mat');
% load('initialG_1600keV_flat_pre_collapse.mat');

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


load('alphas_vA_all_collapse_Glisa_fc0p8h2_G160414.mat');

alphas_ejected(alphas_ejected>1)=1;

alphas_pos_psi=interp2(scale_X,scale_Z,psi_norm_XZsmall_map',alphas_pos_x,alphas_pos_z,'*linear');
alphas_value_psi=interp2(scale_X,scale_Z,psi_XZsmall_map',alphas_pos_x,alphas_pos_z,'*linear');
theta_value=interp2(scale_X,scale_Z,theta_XZsmall_map',alphas_pos_x,alphas_pos_z,'*nearest');






figure(3)

subplot(3,1,2)
POPULATION_RANGE=ALL_TRAPPED_POP.*(alphas_Ekin>EKIN_INF).*(alphas_Ekin<EKIN_SUP);
% POPULATION_RANGE=(alphas_Ekin<EKIN_INF);
ENERGY_POPULATION_0=find(POPULATION_RANGE);
ENERGY_POPULATION=find((~alphas_ejected).*POPULATION_RANGE);

set(gca,'FontSize',20);
hold on;
grid on;
title('trapped (E_{kin} > E_{crit})')



Npsi0=histc(psi_pos0(ENERGY_POPULATION_0),psi_bin_pos);
Npsi=histc(alphas_psi(ENERGY_POPULATION),psi_bin_pos);



density_radial_ini=Npsi0./volume_radial_bin;
density_radial_end=Npsi./volume_radial_bin;

% plot(psi_bin_pos+0.5*radial_bin_size,density_radial_ini_co,'k-.','LineWidth',2);
% plot(radial_r_value_flux(round(psi_bin_pos+0.5*radial_bin_size)),density_radial_ini/max(density_radial_ini),'b--','LineWidth',2);
% plot(radial_r_value_flux(round(psi_bin_pos+0.5*radial_bin_size)),density_radial_end/max(density_radial_ini)-density_radial_ini/max(density_radial_ini),'r','LineWidth',2);
plot(rho_scale(round(psi_bin_pos+0.5*radial_bin_size)),density_radial_end./density_radial_ini-1,'r','LineWidth',2);
% plot(psi_bin_pos+0.5*radial_bin_size,density_radial_end_pas/mean(density_radial_ini_pas),'r','LineWidth',3);
% xlabel('r (m)')
% ylabel('n/max(n)')

xlim([0 0.6])
% h_legend=legend('ini','post');
% set(h_legend,'FontSize',16);
% title('flux surface position (\psi)');

ylim([-0.05 0.31])



subplot(3,1,3)
set(gca,'FontSize',20);
hold on;
grid on;
% title('passing E_{crit} <E_{kin} < 1200 keV')
% title('passing (E_{kin} >  1.6MeV)')
% title('passing (E_{kin} > E_{crit})')
title('passing (E_{kin} >  1.5MeV)')

POPULATION_RANGE=CO_PASSING_POP.*(alphas_Ekin>=CRITICAL_ENERGY);
ENERGY_POPULATION_0=find(POPULATION_RANGE);
ENERGY_POPULATION_0_INF=find(POPULATION_RANGE.*(alphas_Ekin<CRITICAL_ENERGY_PASSING));
ENERGY_POPULATION_0_SUP=find(POPULATION_RANGE.*(alphas_Ekin>=CRITICAL_ENERGY_PASSING));
ENERGY_POPULATION_INF=find((~alphas_ejected).*POPULATION_RANGE.*(alphas_Ekin<CRITICAL_ENERGY_PASSING));
ENERGY_POPULATION_SUP=find((~alphas_ejected).*POPULATION_RANGE.*(alphas_Ekin>=CRITICAL_ENERGY_PASSING));


radial_bin_size_half=8;
radial_bin_size=radial_bin_size_half*2;
psi_bin_pos=(1:radial_bin_size:255-radial_bin_size);
nb_psi_bins=length(psi_bin_pos)


Npsi0_inf=histc(psi_pos0(ENERGY_POPULATION_0_INF),psi_bin_pos);
Npsi0_sup=histc(psi_pos0(ENERGY_POPULATION_0_SUP),psi_bin_pos);
Npsi_inf=histc(alphas_psi(ENERGY_POPULATION_INF),psi_bin_pos);
Npsi_sup=histc(alphas_psi(ENERGY_POPULATION_SUP),psi_bin_pos);

for (psi=1:256)
    volume_radial(psi)=sum(volume_tor_diff(psi,:));
end
volume_radial_bin=zeros(nb_psi_bins,1);
for (psi=1:nb_psi_bins)
    volume_radial_bin(psi)=sum(volume_radial((psi-1)*radial_bin_size+1:psi*radial_bin_size));
end

density_radial_ini_inf=Npsi0_inf./volume_radial_bin;
density_radial_ini_sup=Npsi0_sup./volume_radial_bin;
density_radial_end_inf=Npsi_inf./volume_radial_bin;
density_radial_end_sup=Npsi_sup./volume_radial_bin;

% plot(radial_r_value_flux(round(psi_bin_pos+0.5*radial_bin_size)),density_radial_end_sup/max(density_radial_ini_sup)-density_radial_ini_sup/max(density_radial_ini_sup),'r','LineWidth',2);
plot(rho_scale(round(psi_bin_pos+0.5*radial_bin_size)),density_radial_end_sup./density_radial_ini_sup-1,'r','LineWidth',2);


POPULATION_RANGE=COUNTER_PASSING_POP.*(alphas_Ekin>=CRITICAL_ENERGY);
ENERGY_POPULATION_0=find(POPULATION_RANGE);
ENERGY_POPULATION_0_INF=find(POPULATION_RANGE.*(alphas_Ekin<CRITICAL_ENERGY_PASSING));
ENERGY_POPULATION_0_SUP=find(POPULATION_RANGE.*(alphas_Ekin>=CRITICAL_ENERGY_PASSING));
ENERGY_POPULATION_INF=find((~alphas_ejected).*POPULATION_RANGE.*(alphas_Ekin<CRITICAL_ENERGY_PASSING));
ENERGY_POPULATION_SUP=find((~alphas_ejected).*POPULATION_RANGE.*(alphas_Ekin>=CRITICAL_ENERGY_PASSING));


Npsi0_inf=histc(psi_pos0(ENERGY_POPULATION_0_INF),psi_bin_pos);
Npsi0_sup=histc(psi_pos0(ENERGY_POPULATION_0_SUP),psi_bin_pos);
Npsi_inf=histc(alphas_psi(ENERGY_POPULATION_INF),psi_bin_pos);
Npsi_sup=histc(alphas_psi(ENERGY_POPULATION_SUP),psi_bin_pos);

for (psi=1:256)
    volume_radial(psi)=sum(volume_tor_diff(psi,:));
end
volume_radial_bin=zeros(nb_psi_bins,1);
for (psi=1:nb_psi_bins)
    volume_radial_bin(psi)=sum(volume_radial((psi-1)*radial_bin_size+1:psi*radial_bin_size));
end

density_radial_ini_inf=Npsi0_inf./volume_radial_bin;
density_radial_ini_sup=Npsi0_sup./volume_radial_bin;
density_radial_end_inf=Npsi_inf./volume_radial_bin;
density_radial_end_sup=Npsi_sup./volume_radial_bin;

% plot(radial_r_value_flux(round(psi_bin_pos+0.5*radial_bin_size)),density_radial_end_sup/max(density_radial_ini_sup)-density_radial_ini_sup/max(density_radial_ini_sup),'b--','LineWidth',2);
plot(rho_scale(round(psi_bin_pos+0.5*radial_bin_size)),density_radial_end_sup./density_radial_ini_sup-1,'b--','LineWidth',2);




xlabel('[q_0=0.88] r (m)')
% ylabel('n/max(n)')
% ylim([-0.2 0.15])
ylim([-0.2 0.2])

xlim([0 0.6])
h_legend=legend('co-passing','counter-passing');
set(h_legend,'FontSize',14);

% legend('co-passing','counter-passing')
% h_legend=legend('ini','post');
% set(h_legend,'FontSize',16);
% h_legend=legend('ini (E_{crit}<E_{kin}<E_{s})','ini (E_{kin}>E_{s})','post (E_{crit}<E_{kin}<E_{s})','post (E_{kin}>E_{s})');
% set(h_legend,'FontSize',16);
% title('flux surface position (\psi)');

