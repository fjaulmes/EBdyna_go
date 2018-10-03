close all;

psi_mix=size_r-4
delta_psi_q1=20
psi_core=round(psi_rank_q1-delta_psi_q1)
psi_outer=round(psi_rank_q1+delta_psi_q1)


load('initialG_800keV_flat_precession_stats.mat')
load('initialG_800keV_flat_pre_collapse.mat');
pphi_ini=alphas_pphi0;
disp('Total number of particles simulates');
Nalphas_simulated=length(alphas_pos_x);
disp(Nalphas_simulated);
X0=alphas_pos_x;
Z0=alphas_pos_z;
alphas_mm0=alphas_mm;
Bfield0=interp2(scale_X,scale_Z,Btot_XZ_map',alphas_pos_x,alphas_pos_z,'*linear');
alphas_Eperp=alphas_mm0.*Bfield0;
alphas_Epll=0.5*(mHe/eV)*alphas_vpll.^2;
alphas_Ekin=alphas_Eperp+alphas_Epll;
Ekin0=alphas_Ekin;
alphas_lambda0=Bavg*alphas_mm0./Ekin0;
alphas_pos_psi=interp2(scale_X,scale_Z,psi_norm_XZsmall_map',alphas_pos_x,alphas_pos_z,'*linear');
alphas_psi0=alphas_pos_psi;

load('flat800keV_collapse_Glocal_fc1h2_G260713.mat');
pphi_end=alphas_pphi0;
alphas_pos_psi=interp2(scale_X,scale_Z,psi_norm_XZsmall_map',alphas_pos_x,alphas_pos_z,'*linear');
EKIN_INF=200
EKIN_SUP=2200
% ZONE_RANGE=(alphas_psi0>=psi_core).*(alphas_psi0<=psi_outer);
ZONE_RANGE=(alphas_psi0<=psi_outer);
POPULATION_RANGE=(alphas_Ekin>EKIN_INF*1e3).*(alphas_Ekin<EKIN_SUP*1e3).*ZONE_RANGE;
ENERGY_POPULATION_0=find(POPULATION_RANGE);
ENERGY_POPULATION=find((~alphas_ejected).*POPULATION_RANGE);
post_dist_rescale=size(ENERGY_POPULATION_0,1)/size(ENERGY_POPULATION,1);

disp('EVALUATED_POPULATION')
EVALUATED_POPULATION=length(ENERGY_POPULATION);
disp(EVALUATED_POPULATION);


% alphas_pphi0=(mHe*(X0+R0).*(vpll0)-2*eV*(psi_value0))/eV;
% alphas_pphi=(mHe*(alphas_pos_x+R0).*(alphas_vpll)-2*eV*(alphas_psi_value_corr))/eV;

% Npphi0=histc(p_phi0(ENERGY_POPULATION)/eV,p_phi_scale-0.5*p_phi_bin_size);
% Npphi=histc(p_phi(ENERGY_POPULATION)/eV,p_phi_scale-0.5*p_phi_bin_size);

Delta_pphi=pphi_end-pphi_ini;




figure(6);
Dp_phi_bin_size=0.1;
Dp_phi_scale=(-0.9:Dp_phi_bin_size:0.9);

% subplot(2,1,1)
set(gca,'FontSize',26);
hold on;
grid on;

xlim([-0.9 0.9])

NDpphi=histc(Delta_pphi(find(POPULATION_RANGE.*ALL_TRAPPED_POP)),Dp_phi_scale-0.5*Dp_phi_bin_size);
plot(Dp_phi_scale,NDpphi/max(NDpphi),'g--','LineWidth',3);
NDpphi=histc(Delta_pphi(find(POPULATION_RANGE.*CO_PASSING_POP)),Dp_phi_scale-0.5*Dp_phi_bin_size);
plot(Dp_phi_scale,NDpphi/max(NDpphi),'r','LineWidth',3);
NDpphi=histc(Delta_pphi(find(POPULATION_RANGE.*COUNTER_PASSING_POP)),Dp_phi_scale-0.5*Dp_phi_bin_size);
plot(Dp_phi_scale,NDpphi/max(NDpphi),'b-.','LineWidth',3);

xl=xlabel('$$\Delta$$p$$_\varphi$$','interpreter','latex')
set(xl,'Interpreter','latex');

% xlabel('\Delta p_\phi')
ylabel('N/max(N)')


Dp_phi_bin_size=0.04;
Dp_phi_scale=(-0.9:Dp_phi_bin_size:0.9);

legend('trapped','co-passing','counter-passing')

% POPULATION_RANGE=(alphas_Ekin>EKIN_SUP*1e3).*ZONE_RANGE;
% ENERGY_POPULATION_0=find(POPULATION_RANGE);
% ENERGY_POPULATION=find((~alphas_ejected).*POPULATION_RANGE);
% NDpphi=histc(Delta_pphi(ENERGY_POPULATION),Dp_phi_scale-0.02);
% 
% subplot(2,1,2)
% set(gca,'FontSize',26);
% hold on;
% grid on;
% xlim([-0.9 0.9])
% 
% NDpphi=histc(Delta_pphi(find(POPULATION_RANGE.*ALL_TRAPPED_POP)),Dp_phi_scale-0.5*Dp_phi_bin_size);
% plot(Dp_phi_scale,NDpphi/max(NDpphi),'g--','LineWidth',3);
% NDpphi=histc(Delta_pphi(find(POPULATION_RANGE.*CO_PASSING_POP)),Dp_phi_scale-0.5*Dp_phi_bin_size);
% plot(Dp_phi_scale,NDpphi/max(NDpphi),'r','LineWidth',3);
% NDpphi=histc(Delta_pphi(find(POPULATION_RANGE.*COUNTER_PASSING_POP)),Dp_phi_scale-0.5*Dp_phi_bin_size);
% plot(Dp_phi_scale,NDpphi/max(NDpphi),'b-.','LineWidth',3);
% 
% xlabel('\Delta p_\phi')
% ylabel('N/max(N)')
% legend('trapped','co-passing','counter-passing')
