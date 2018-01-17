
close all;

rmix=radial_r_value_flux(size_r-4)
psi_mix=size_r-4
delta_psi_q1=20;
psi_core=round(psi_rank_q1-delta_psi_q1)
psi_outer=max(round(1.1*psi_mix))

PRECESS_LAMBDA_BIN_SIZE=0.018;
precession_lambda_bins=(32*PRECESS_LAMBDA_BIN_SIZE:PRECESS_LAMBDA_BIN_SIZE:62*PRECESS_LAMBDA_BIN_SIZE);
lambda_values=(32.5*PRECESS_LAMBDA_BIN_SIZE:PRECESS_LAMBDA_BIN_SIZE:61.5*PRECESS_LAMBDA_BIN_SIZE);

lambda_precess_value=zeros(length(precession_lambda_bins)-1,5);
lambda_Dpphi_value_pos=zeros(length(precession_lambda_bins)-1,5);
lambda_Dpphi_value_neg=zeros(length(precession_lambda_bins)-1,5);
lambda_omega_r_value=zeros(length(precession_lambda_bins)-1,5);
E_Dpphi_value_pos=zeros(5,1);
E_Dpphi_value_neg=zeros(5,1);
Ekin_values=[5 16 32 48 64 90 ]

ORBIT_POP_name='ALL_TRAPPED_POP'

Ebin=1
load('initialG_flatD_5keV_pre_collapse.mat')
load('initialG_flatD_5keV_precession_stats.mat')
omega_phi_avg=-(omega_precess_avg+omega_psi_avg);
ORBIT_POP=evalin('base',ORBIT_POP_name);
REGION_POP=(alphas_psi<=psi_outer);
REGION_POP1=(alphas_psi<=psi_core);
REGION_POP2=(alphas_psi>=psi_core).*(alphas_psi<=psi_outer);
pphi_ini=alphas_pphi0;
lambda_ini=Bavg*alphas_mm./alphas_Ekin;
for bin=1:(length(precession_lambda_bins)-1)
    PRECESS_LAMBDA_POP=find((alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*ALL_TRAPPED_POP);
    lambda_precess_value(bin,Ebin)=mean((omega_precess_avg(PRECESS_LAMBDA_POP)));
    lambda_omega_r_value(bin,Ebin)=mean(delta_r_avg(PRECESS_LAMBDA_POP));
end
load('flatD_5keV_collapse_Glisa_fc2h2_G290114.mat')
lambda_end=Bavg*alphas_mm./alphas_Ekin;
pphi_end=alphas_pphi0;
delta_pphi=pphi_end-pphi_ini;
delta_lambda=lambda_end-lambda_ini;
ALPHAS_POP=(~alphas_ejected);
for bin=1:(length(precession_lambda_bins)-1)
    PRECESS_LAMBDA_POP=find(REGION_POP1.*(alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*ALL_TRAPPED_POP.*ALPHAS_POP);
    lambda_Dpphi_value_pos(bin,Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP)));
    PRECESS_LAMBDA_POP=find(REGION_POP2.*(alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*ALL_TRAPPED_POP.*ALPHAS_POP);
    lambda_Dpphi_value_neg(bin,Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP)));
end
E_Dpphi_value_pos(Ebin)=mean((delta_pphi(find(ALL_TRAPPED_POP.*REGION_POP1))));
E_Dpphi_value_neg(Ebin)=mean((delta_pphi(find(ALL_TRAPPED_POP.*REGION_POP2))));


% Ebin=1
% load('initialG_flatHe_16keV_pre_collapse.mat')
% load('initialG_flatHe_16keV_precession_stats.mat')
% omega_phi_avg=-(omega_precess_avg+omega_psi_avg);
% ORBIT_POP=evalin('base',ORBIT_POP_name);
% REGION_POP=(alphas_psi<=psi_outer);
% pphi_ini=alphas_pphi0;
% lambda_ini=Bavg*alphas_mm./alphas_Ekin;
% for bin=1:(length(precession_lambda_bins)-1)
%     PRECESS_LAMBDA_POP=find((alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*ALL_TRAPPED_POP);
%     lambda_precess_value(bin,Ebin)=mean((omega_precess_avg(PRECESS_LAMBDA_POP)));
%     lambda_omega_r_value(bin,Ebin)=mean(delta_r_avg(PRECESS_LAMBDA_POP));
% end
% load('flatHe_16keV_collapse_Glisa_fc2h2_G290114.mat')
% lambda_end=Bavg*alphas_mm./alphas_Ekin;
% pphi_end=alphas_pphi0;
% delta_pphi=pphi_end-pphi_ini;
% delta_lambda=lambda_end-lambda_ini;
% ALPHAS_POP=(~alphas_ejected);
% for bin=1:(length(precession_lambda_bins)-1)
%     PRECESS_LAMBDA_POP=find((delta_pphi>=0).*(alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*ALL_TRAPPED_POP.*ALPHAS_POP);
%     lambda_Dpphi_value_pos(bin,Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP)));
%     PRECESS_LAMBDA_POP=find((delta_pphi<=0).*(alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*ALL_TRAPPED_POP.*ALPHAS_POP);
%     lambda_Dpphi_value_neg(bin,Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP)));
% end
% E_Dpphi_value_pos(Ebin)=mean((delta_pphi(find(ALL_TRAPPED_POP.*REGION_POP.*(delta_pphi>0)))));
% E_Dpphi_value_neg(Ebin)=mean((delta_pphi(find(ALL_TRAPPED_POP.*REGION_POP.*(delta_pphi<0)))));


Ebin=Ebin+1
load('initialG_flatD_16keV_pre_collapse.mat')
load('initialG_flatD_16keV_precession_stats.mat')
omega_phi_avg=-(omega_precess_avg+omega_psi_avg);
ORBIT_POP=evalin('base',ORBIT_POP_name);
REGION_POP=(alphas_psi<=psi_outer);
REGION_POP1=(alphas_psi<=psi_core);
REGION_POP2=(alphas_psi>=psi_core).*(alphas_psi<=psi_outer);
pphi_ini=alphas_pphi0;
lambda_ini=Bavg*alphas_mm./alphas_Ekin;
for bin=1:(length(precession_lambda_bins)-1)
    PRECESS_LAMBDA_POP=find((alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*ALL_TRAPPED_POP);
    lambda_precess_value(bin,Ebin)=mean((omega_precess_avg(PRECESS_LAMBDA_POP)));
    lambda_omega_r_value(bin,Ebin)=mean(delta_r_avg(PRECESS_LAMBDA_POP));
end
load('flatD_16keV_collapse_Glisa_fc2h2_G290114.mat')
lambda_end=Bavg*alphas_mm./alphas_Ekin;
pphi_end=alphas_pphi0;
delta_pphi=pphi_end-pphi_ini;
delta_lambda=lambda_end-lambda_ini;
ALPHAS_POP=(~alphas_ejected);
for bin=1:(length(precession_lambda_bins)-1)
    PRECESS_LAMBDA_POP=find(REGION_POP1.*(alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*ALL_TRAPPED_POP.*ALPHAS_POP);
    lambda_Dpphi_value_pos(bin,Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP)));
    PRECESS_LAMBDA_POP=find(REGION_POP2.*(alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*ALL_TRAPPED_POP.*ALPHAS_POP);
    lambda_Dpphi_value_neg(bin,Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP)));
end
E_Dpphi_value_pos(Ebin)=mean((delta_pphi(find(ALL_TRAPPED_POP.*REGION_POP1))));
E_Dpphi_value_neg(Ebin)=mean((delta_pphi(find(ALL_TRAPPED_POP.*REGION_POP2))));



Ebin=Ebin+1
load('initialG_flatD_32keV_pre_collapse.mat')
load('initialG_flatD_32keV_precession_stats.mat')
omega_phi_avg=-(omega_precess_avg+omega_psi_avg);
ORBIT_POP=evalin('base',ORBIT_POP_name);
REGION_POP=(alphas_psi<=psi_outer);
REGION_POP1=(alphas_psi<=psi_core);
REGION_POP2=(alphas_psi>=psi_core).*(alphas_psi<=psi_outer);
pphi_ini=alphas_pphi0;
lambda_ini=Bavg*alphas_mm./alphas_Ekin;
for bin=1:(length(precession_lambda_bins)-1)
    PRECESS_LAMBDA_POP=find((alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*ALL_TRAPPED_POP);
    lambda_precess_value(bin,Ebin)=mean((omega_precess_avg(PRECESS_LAMBDA_POP)));
    lambda_omega_r_value(bin,Ebin)=mean(delta_r_avg(PRECESS_LAMBDA_POP));
end
load('flatD_32keV_collapse_Glisa_fc2h2_G290114.mat')
lambda_end=Bavg*alphas_mm./alphas_Ekin;
pphi_end=alphas_pphi0;
delta_pphi=pphi_end-pphi_ini;
delta_lambda=lambda_end-lambda_ini;
ALPHAS_POP=(~alphas_ejected);
for bin=1:(length(precession_lambda_bins)-1)
    PRECESS_LAMBDA_POP=find(REGION_POP1.*(alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*ALL_TRAPPED_POP.*ALPHAS_POP);
    lambda_Dpphi_value_pos(bin,Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP)));
    PRECESS_LAMBDA_POP=find(REGION_POP2.*(alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*ALL_TRAPPED_POP.*ALPHAS_POP);
    lambda_Dpphi_value_neg(bin,Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP)));
end
E_Dpphi_value_pos(Ebin)=mean((delta_pphi(find(ALL_TRAPPED_POP.*REGION_POP1))));
E_Dpphi_value_neg(Ebin)=mean((delta_pphi(find(ALL_TRAPPED_POP.*REGION_POP2))));


Ebin=Ebin+1
load('initialG_flatD_48keV_pre_collapse.mat')
load('initialG_flatD_48keV_precession_stats.mat')
omega_phi_avg=-(omega_precess_avg+omega_psi_avg);
ORBIT_POP=evalin('base',ORBIT_POP_name);
REGION_POP=(alphas_psi<=psi_outer);
REGION_POP1=(alphas_psi<=psi_core);
REGION_POP2=(alphas_psi>=psi_core).*(alphas_psi<=psi_outer);
pphi_ini=alphas_pphi0;
lambda_ini=Bavg*alphas_mm./alphas_Ekin;
for bin=1:(length(precession_lambda_bins)-1)
    PRECESS_LAMBDA_POP=find((alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*ALL_TRAPPED_POP);
    lambda_precess_value(bin,Ebin)=mean((omega_precess_avg(PRECESS_LAMBDA_POP)));
    lambda_omega_r_value(bin,Ebin)=mean(delta_r_avg(PRECESS_LAMBDA_POP));
end
load('flatD_48keV_collapse_Glisa_fc2h2_G290114.mat')
lambda_end=Bavg*alphas_mm./alphas_Ekin;
pphi_end=alphas_pphi0;
delta_pphi=pphi_end-pphi_ini;
delta_lambda=lambda_end-lambda_ini;
ALPHAS_POP=(~alphas_ejected);
for bin=1:(length(precession_lambda_bins)-1)
    PRECESS_LAMBDA_POP=find(REGION_POP1.*(alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*ALL_TRAPPED_POP.*ALPHAS_POP);
    lambda_Dpphi_value_pos(bin,Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP)));
    PRECESS_LAMBDA_POP=find(REGION_POP2.*(alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*ALL_TRAPPED_POP.*ALPHAS_POP);
    lambda_Dpphi_value_neg(bin,Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP)));
end
E_Dpphi_value_pos(Ebin)=mean((delta_pphi(find(ALL_TRAPPED_POP.*REGION_POP1))));
E_Dpphi_value_neg(Ebin)=mean((delta_pphi(find(ALL_TRAPPED_POP.*REGION_POP2))));




Ebin=Ebin+1
load('initialG_flatD_64keV_pre_collapse.mat')
load('initialG_flatD_64keV_precession_stats.mat')
% load('initialG_flat_He_800keV_pre_collapse_old.mat')
% load('initialG_flat_He_800keV_precession_stats_old.mat')
omega_phi_avg=-(omega_precess_avg+omega_psi_avg);
ORBIT_POP=evalin('base',ORBIT_POP_name);
REGION_POP=(alphas_psi<=psi_outer);
REGION_POP1=(alphas_psi<=psi_core);
REGION_POP2=(alphas_psi>=psi_core).*(alphas_psi<=psi_outer);
pphi_ini=alphas_pphi0;
lambda_ini=Bavg*alphas_mm./alphas_Ekin;
for bin=1:(length(precession_lambda_bins)-1)
    PRECESS_LAMBDA_POP=find((alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*ALL_TRAPPED_POP);
    lambda_precess_value(bin,Ebin)=mean((omega_precess_avg(PRECESS_LAMBDA_POP)));
    lambda_omega_r_value(bin,Ebin)=mean(delta_r_avg(PRECESS_LAMBDA_POP));
end
load('flatD_64keV_collapse_Glisa_fc2h2_G290114.mat')
% load('flat_He_800keV_collapse_Glisa_fc1h2_G090314.mat')
lambda_end=Bavg*alphas_mm./alphas_Ekin;
pphi_end=alphas_pphi0;
delta_pphi=pphi_end-pphi_ini;
delta_lambda=lambda_end-lambda_ini;
ALPHAS_POP=(~alphas_ejected);
for bin=1:(length(precession_lambda_bins)-1)
    PRECESS_LAMBDA_POP=find(REGION_POP1.*(alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*ALL_TRAPPED_POP.*ALPHAS_POP);
    lambda_Dpphi_value_pos(bin,Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP)));
    PRECESS_LAMBDA_POP=find(REGION_POP2.*(alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*ALL_TRAPPED_POP.*ALPHAS_POP);
    lambda_Dpphi_value_neg(bin,Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP)));
end
E_Dpphi_value_pos(Ebin)=mean((delta_pphi(find(ALL_TRAPPED_POP.*REGION_POP1))));
E_Dpphi_value_neg(Ebin)=mean((delta_pphi(find(ALL_TRAPPED_POP.*REGION_POP2))));



Ebin=Ebin+1
load('initialG_flatD_90keV_pre_collapse.mat')
load('initialG_flatD_90keV_precession_stats.mat')
omega_phi_avg=-(omega_precess_avg+omega_psi_avg);
ORBIT_POP=evalin('base',ORBIT_POP_name);
REGION_POP=(alphas_psi<=psi_outer);
REGION_POP1=(alphas_psi<=psi_core);
REGION_POP2=(alphas_psi>=psi_core).*(alphas_psi<=psi_outer);
pphi_ini=alphas_pphi0;
lambda_ini=Bavg*alphas_mm./alphas_Ekin;
for bin=1:(length(precession_lambda_bins)-1)
    PRECESS_LAMBDA_POP=find((alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*ALL_TRAPPED_POP);
    lambda_precess_value(bin,Ebin)=mean((omega_precess_avg(PRECESS_LAMBDA_POP)));
    lambda_omega_r_value(bin,Ebin)=mean(delta_r_avg(PRECESS_LAMBDA_POP));
end
load('flatD_90keV_collapse_Glisa_fc2h2_G290114.mat')
lambda_end=Bavg*alphas_mm./alphas_Ekin;
pphi_end=alphas_pphi0;
delta_pphi=pphi_end-pphi_ini;
delta_lambda=lambda_end-lambda_ini;
ALPHAS_POP=(~alphas_ejected);
for bin=1:(length(precession_lambda_bins)-1)
    PRECESS_LAMBDA_POP=find(REGION_POP1.*(alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*ALL_TRAPPED_POP.*ALPHAS_POP);
    lambda_Dpphi_value_pos(bin,Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP)));
    PRECESS_LAMBDA_POP=find(REGION_POP2.*(alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*ALL_TRAPPED_POP.*ALPHAS_POP);
    lambda_Dpphi_value_neg(bin,Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP)));
end
E_Dpphi_value_pos(Ebin)=mean((delta_pphi(find(ALL_TRAPPED_POP.*REGION_POP1))));
E_Dpphi_value_neg(Ebin)=mean((delta_pphi(find(ALL_TRAPPED_POP.*REGION_POP2))));


E_Dpphi0_pos=abs(E_Dpphi_value_pos(1))
E_Dpphi0_neg=abs(E_Dpphi_value_neg(1))


%%
figure(8)
set(gca,'FontSize',28);
hold on
grid on
plot(Ekin_values,E_Dpphi_value_pos/E_Dpphi0_pos,'r--','Linewidth',3)
plot(Ekin_values,E_Dpphi_value_neg/E_Dpphi0_neg,'b','Linewidth',3)
yl=ylabel('$$  \Delta p_\varphi  / |\Delta p_\varphi _{th}|$$','interpreter','latex');
set(yl,'interpreter','latex');
xl=xlabel('$$  \mathcal{E}_{kin} (keV)  $$','interpreter','latex');
set(xl,'interpreter','latex');
% xlim([65 1600])

legend('region 1 (core)','region 2 (q=1)')

% Ebin=5
% 
% load('initialG_2800keV_flat_pre_collapse.mat')
% load('initialG_2800keV_flat_precession_stats.mat')
% ORBIT_POP=evalin('base',ORBIT_POP_name);
% pphi_ini=alphas_pphi0;
% lambda_ini=Bavg*alphas_mm./alphas_Ekin;
% for bin=1:(length(precession_lambda_bins)-1)
%     PRECESS_LAMBDA_POP=find((alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*ALL_TRAPPED_POP);
%     lambda_precess_value(bin,Ebin)=mean((omega_precess_avg(PRECESS_LAMBDA_POP)));
% end
% load('flat2800keV_collapse_Glisa_fc1h2_G230713.mat')
% lambda_end=Bavg*alphas_mm./alphas_Ekin;
% pphi_end=alphas_pphi0;
% delta_pphi=pphi_end-pphi_ini;
% delta_lambda=lambda_end-lambda_ini;
% ALPHAS_POP=(~alphas_ejected);
% for bin=1:(length(precession_lambda_bins)-1)
%     PRECESS_LAMBDA_POP=find((alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*ALL_TRAPPED_POP.*ALPHAS_POP);
%     lambda_Dpphi_value(bin,Ebin)=mean(abs(delta_pphi(PRECESS_LAMBDA_POP)));
% end

%%
close all;

figure(1)
subplot(4,1,1)
set(gca,'FontSize',20);
% title('Trapped particles');

hold on
grid on


% xlabel('\lambda_{0}')
yl=ylabel('$$  \Delta p_\varphi (1) $$','interpreter','latex');
set(yl,'interpreter','latex');


% plot(lambda_values,lambda_Dpphi_value_pos(:,1),'g-+','LineWidth',2);
plot(lambda_values,lambda_Dpphi_value_pos(:,2),'b--','LineWidth',3);
plot(lambda_values,lambda_Dpphi_value_pos(:,3),'k-.','LineWidth',4);
plot(lambda_values,lambda_Dpphi_value_pos(:,4),'r','LineWidth',5);
plot(lambda_values,lambda_Dpphi_value_pos(:,5),'--','color',[1.0 0.8 0.3],'LineWidth',5);

xlim([0.82 1.08])
% ylim([-0.22 0.02])
% legend('\sim 65 keV','200 keV','400 keV','800 keV')

% ylim([ min(lambda_Dpphi_value_neg(:,1)) max(lambda_Dpphi_value_pos(:,1))])

subplot(4,1,2)
set(gca,'FontSize',20);
% title('Trapped particles');

hold on
grid on


% xlabel('\lambda_{0}')
yl=ylabel('$$  \Delta p_\varphi (2) $$','interpreter','latex');
set(yl,'interpreter','latex');


% plot(lambda_values,lambda_Dpphi_value_neg(:,1),'g-+','LineWidth',2);
plot(lambda_values,lambda_Dpphi_value_neg(:,2),'b--','LineWidth',3);
plot(lambda_values,lambda_Dpphi_value_neg(:,3),'k-.','LineWidth',4);
plot(lambda_values,lambda_Dpphi_value_neg(:,4),'r','LineWidth',5);
plot(lambda_values,lambda_Dpphi_value_neg(:,5),'--','color',[1.0 0.8 0.3],'LineWidth',5);

xlim([0.82 1.08])
% ylim([-0.02 0.11])



figure(1)
omega_crash=0.5*(pi)/(72*1e-6)
% omega_crash_A=0.7*(pi)/(128*1e-6)

subplot(4,1,3)
set(gca,'FontSize',20);
% title('Trapped particles');

hold on
grid on
ylabel('\omega_{v_D} (rad/s)');

% plot(lambda_values,lambda_precess_value(:,1),'g-+','LineWidth',2);
plot(lambda_values,lambda_precess_value(:,2),'b--','LineWidth',3);
plot(lambda_values,lambda_precess_value(:,3),'k-.','LineWidth',4);
plot(lambda_values,lambda_precess_value(:,4),'r','LineWidth',5);
plot(lambda_values,lambda_precess_value(:,5),'--','color',[1.0 0.8 0.3],'LineWidth',5);
% plot(lambda_values,lambda_values.*0+omega_crash,'color',[0.1 0.5 0.1],'LineWidth',3);
% plot(lambda_values,lambda_values.*0-omega_crash,'color',[0.1 0.5 0.1],'LineWidth',3);
plot(lambda_values,lambda_values.*0+omega_crash,'--','color',[0.1 0.9 0.2],'LineWidth',2);
plot(lambda_values,lambda_values.*0-omega_crash,'--','color',[0.1 0.9 0.2],'LineWidth',2);
% plot(lambda_values,lambda_values.*0+omega_crash_A,'color',[0.1 0.8 0.2],'LineWidth',6);
% plot(lambda_values,lambda_values.*0-omega_crash_A,'color',[0.1 0.8 0.2],'LineWidth',6);
% plot(lambda_values,lambda_precess_value(:,1),'g-+','LineWidth',2);
plot(lambda_values,lambda_precess_value(:,2),'b--','LineWidth',3);
plot(lambda_values,lambda_precess_value(:,3),'k-.','LineWidth',4);
plot(lambda_values,lambda_precess_value(:,4),'r','LineWidth',5);
plot(lambda_values,lambda_precess_value(:,5),'--','color',[1.0 0.8 0.3],'LineWidth',5);


annotation('textbox',...
    [0.5 0.448 0.0483 0.02],...
    'EdgeColor','none',...
    'String',{'+\omega_{cr}'},'fontsize',16,...
    'FitBoxToText','off','color',[0 0.4 0.1]);

annotation('textbox',...
    [0.3 0.392 0.0483 0.02],...
    'EdgeColor','none',...
    'String',{'-\omega_{cr}'},'fontsize',16,...
    'FitBoxToText','off','color',[0 0.4 0.1]);



xlim([0.82 1.08])

% ylim([1.1*min(lambda_precess_value(:,5)) 1.05*max(lambda_precess_value(:,5)) ])
% ylim([1.1*min(lambda_precess_value(:,5)) 0.9*max(lambda_precess_value(:,5)) ])


figure(1)
subplot(4,1,4)
set(gca,'FontSize',20);

hold on
grid on
% xlabel('\lambda_{0}')
% yl=ylabel('$$\delta_{r} (m)','interpreter','latex');
% set(yl,'interpreter','latex');
ylabel('\delta_{r} (m)');

% plot(lambda_values,lambda_omega_r_value(:,1),'g-+','LineWidth',3);
plot(lambda_values,lambda_omega_r_value(:,2),'b--','LineWidth',4);
plot(lambda_values,lambda_omega_r_value(:,3),'k-.','LineWidth',4);
plot(lambda_values,lambda_omega_r_value(:,4),'r','LineWidth',5);
plot(lambda_values,lambda_omega_r_value(:,5),'--','color',[1.0 0.8 0.3],'LineWidth',5);
xlim([0.82 1.08])
ylim([0 0.9*max(lambda_omega_r_value(:,5)) ])
hl=legend('16 keV','32 keV','48 keV','64 keV')
set(hl,'fontsize',18)
xlabel('\lambda_{0}')

