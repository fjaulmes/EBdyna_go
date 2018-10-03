psi_mix=size_r-4
delta_psi_q1=15;
psi_core=psi_rank_q1-delta_psi_q1

radial_bin_size_half=6;
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

% load('initialG_alphas_lowEkin_pre_collapse.mat')
% Npsi0=histc(alphas_psi,psi_bin_pos);
% Npsi=histc(alphas_psi,psi_bin_pos);
% 
% load('alphas_lowEkin_collapse_Glisa_fc1h2_G290713.mat')
% Npsi=histc(alphas_psi,psi_bin_pos);
% 
% 
% density_radial_ini=Npsi0./volume_radial_bin;
% density_radial_end=Npsi./volume_radial_bin;
% 
% figure(3)
% 
% subplot(3,1,1)
% set(gca,'FontSize',26);
% hold on;
% grid on;
% title('thermal')
% plot(radial_r_value_flux(round(psi_bin_pos+0.5*radial_bin_size)),density_radial_ini/max(density_radial_ini),'b--','LineWidth',3);
% plot(radial_r_value_flux(round(psi_bin_pos+0.5*radial_bin_size)),density_radial_end/max(density_radial_ini),'r','LineWidth',3);
% % xlabel('r (m)')
% ylabel('n/max(n)')
% 
% xlim([0 1])
% legend('initial','post collapse');
% 

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


% alphas_pos_psi=interp2(scale_X,scale_Z,psi_norm_XZsmall_map',alphas_pos_x,alphas_pos_z,'*linear');
psi_pos0=alphas_pos_psi;

clear ENERGY_POPULATION ENERGY_POPULATION_0


load('alphas_vA_all_collapse_Glisa_fc0p8h2_G160414.mat')

alphas_ejected(alphas_ejected>1)=1;
EKIN_INF=2
EKIN_INF=700*1e3;

alphas_pos_psi=interp2(scale_X,scale_Z,psi_norm_XZsmall_map',alphas_pos_x,alphas_pos_z,'*linear');
alphas_value_psi=interp2(scale_X,scale_Z,psi_XZsmall_map',alphas_pos_x,alphas_pos_z,'*linear');
theta_value=interp2(scale_X,scale_Z,theta_XZsmall_map',alphas_pos_x,alphas_pos_z,'*nearest');






figure(3)

subplot(3,1,1)
POPULATION_RANGE=ALL_TRAPPED_POP.*(alphas_Ekin>EKIN_INF);
ENERGY_POPULATION_0=find(POPULATION_RANGE);
ENERGY_POPULATION=find((~alphas_ejected).*POPULATION_RANGE);

set(gca,'FontSize',26);
hold on;
grid on;
title('trapped')



Npsi0=histc(psi_pos0(ENERGY_POPULATION_0),psi_bin_pos);
Npsi=histc(alphas_psi(ENERGY_POPULATION),psi_bin_pos);



density_radial_ini=Npsi0./volume_radial_bin;
density_radial_end=Npsi./volume_radial_bin;

% plot(psi_bin_pos+0.5*radial_bin_size,density_radial_ini_co,'k-.','LineWidth',2);
plot(radial_r_value_flux(round(psi_bin_pos+0.5*radial_bin_size)),density_radial_ini/max(density_radial_ini),'b--','LineWidth',3);
plot(radial_r_value_flux(round(psi_bin_pos+0.5*radial_bin_size)),density_radial_end/max(density_radial_ini),'r','LineWidth',3);

plot([r_value_q1_mean r_value_q1_mean],[-1 1],'g--','linewidth',3)
plot([r_mix r_mix],[-1 1],'g--','linewidth',3)

% plot(psi_bin_pos+0.5*radial_bin_size,density_radial_end_pas/mean(density_radial_ini_pas),'r','LineWidth',3);
% xlabel('r (m)')
ylabel('n/max(n)')

xlim([0.5 1.4])
ylim([0.1 1])
legend('initial','post collapse');
% title('flux surface position (\psi)');




subplot(3,1,2)
POPULATION_RANGE=CO_PASSING_POP.*(alphas_Ekin>EKIN_INF);
ENERGY_POPULATION_0=find(POPULATION_RANGE);
ENERGY_POPULATION=find((~alphas_ejected).*POPULATION_RANGE);

set(gca,'FontSize',26);
hold on;
grid on;
plot([r_value_q1_mean r_value_q1_mean],[-1 1],'g--','linewidth',3)
plot([r_mix r_mix],[-1 1],'g--','linewidth',3)

title('co passing')

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

% plot(psi_bin_pos+0.5*radial_bin_size,density_radial_ini_co,'k-.','LineWidth',2);
plot(radial_r_value_flux(round(psi_bin_pos+0.5*radial_bin_size)),density_radial_ini/max(density_radial_ini),'b--','LineWidth',3);
plot(radial_r_value_flux(round(psi_bin_pos+0.5*radial_bin_size)),density_radial_end/max(density_radial_ini),'r','LineWidth',3);
% plot(psi_bin_pos+0.5*radial_bin_size,density_radial_end_pas/mean(density_radial_ini_pas),'r','LineWidth',3);
% xlabel('r (m)')
ylabel('n/max(n)')

xlim([0.5 1.4])
ylim([0.1 1])



subplot(3,1,3)
POPULATION_RANGE=COUNTER_PASSING_POP.*(alphas_Ekin>EKIN_INF);
ENERGY_POPULATION_0=find(POPULATION_RANGE);
ENERGY_POPULATION=find((~alphas_ejected).*POPULATION_RANGE);

set(gca,'FontSize',26);
hold on;
grid on;

title('counter passing')

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

% plot(psi_bin_pos+0.5*radial_bin_size,density_radial_ini_co,'k-.','LineWidth',2);
plot(radial_r_value_flux(round(psi_bin_pos+0.5*radial_bin_size)),density_radial_ini/max(density_radial_ini),'b--','LineWidth',3);
plot(radial_r_value_flux(round(psi_bin_pos+0.5*radial_bin_size)),density_radial_end/max(density_radial_ini),'r','LineWidth',3);
% plot(psi_bin_pos+0.5*radial_bin_size,density_radial_end_pas/mean(density_radial_ini_pas),'r','LineWidth',3);
xlabel('r (m)')
ylabel('n/max(n)')
plot([r_value_q1_mean r_value_q1_mean],[-1 1],'g--','linewidth',3)
plot([r_mix r_mix],[-1 1],'g--','linewidth',3)

xlim([0.5 1.4])
ylim([0.1 1])
legend('initial','post collapse');
% title('flux surface position (\psi)');

