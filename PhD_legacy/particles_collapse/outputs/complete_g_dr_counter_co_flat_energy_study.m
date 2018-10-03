% reset_data_analysis_environment

close all;

rmix=radial_r_value_flux(size_r-4)

PRECESS_LAMBDA_BIN_SIZE=0.1;
precession_lambda_bins=(1*PRECESS_LAMBDA_BIN_SIZE:PRECESS_LAMBDA_BIN_SIZE:12*PRECESS_LAMBDA_BIN_SIZE);
lambda_values=(1.5*PRECESS_LAMBDA_BIN_SIZE:PRECESS_LAMBDA_BIN_SIZE:11.5*PRECESS_LAMBDA_BIN_SIZE);

lambda_dr_value=zeros(length(precession_lambda_bins)-1,5);
lambda_g_value=zeros(length(precession_lambda_bins)-1,5);
lambda_Dpphi_value=zeros(length(precession_lambda_bins)-1,5);

g_Ekin_counter=zeros(5,1);
g_Ekin_co=zeros(5,1);
dr_Ekin_counter=zeros(5,1);
dr_Ekin_co=zeros(5,1);
omega_vD_Ekin_counter=zeros(5,1);
omega_vD_Ekin_cor=zeros(5,1);

Dpphi_Ekin_counter=zeros(5,1);
Dpphi_Ekin_co=zeros(5,1);

ORBIT_POP_name='ALL_PASSING_POP'

Ebin=1
load('initialG_alphas_lowEkin_precession_stats.mat')
load('initialG_alphas_lowEkin_pre_collapse.mat')
ORBIT_POP=evalin('base',ORBIT_POP_name);
pphi_ini=alphas_pphi0;
lambda_ini=Bavg*alphas_mm./alphas_Ekin;
for bin=1:(length(precession_lambda_bins)-1)
    PRECESS_LAMBDA_POP_CO=find((alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*CO_PASSING_POP);
    PRECESS_LAMBDA_POP_COUNTER=find((alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*COUNTER_PASSING_POP);
    lambda_g_value(bin,Ebin)=(mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER))/mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)))));
    lambda_g_value(bin,Ebin)=lambda_g_value(bin,Ebin)*(mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_CO))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)))));
    lambda_dr_value(bin,Ebin)=mean(delta_r_avg(PRECESS_LAMBDA_POP_COUNTER))/mean(delta_r_avg(PRECESS_LAMBDA_POP_CO));
end
PRECESS_LAMBDA_POP_COUNTER=find(COUNTER_PASSING_POP);
PRECESS_LAMBDA_POP_CO=find(CO_PASSING_POP);
g_Ekin_counter(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)));
g_Ekin_co(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_CO)));
dr_Ekin_counter(Ebin)=mean(delta_r_avg(PRECESS_LAMBDA_POP_COUNTER));
dr_Ekin_co(Ebin)=mean(delta_r_avg(PRECESS_LAMBDA_POP_CO));
omega_vD_Ekin_counter(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)));
omega_vD_Ekin_co(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)));
load('alphas_lowEkin_collapse_Glisa_fc1h2_G250713.mat')
lambda_end=Bavg*alphas_mm./alphas_Ekin;
pphi_end=alphas_pphi0;
delta_pphi=pphi_end-pphi_ini;
delta_lambda=lambda_end-lambda_ini;
ALPHAS_POP=(~alphas_ejected);
for bin=1:(length(precession_lambda_bins)-1)
    PRECESS_LAMBDA_POP_CO=find((alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*CO_PASSING_POP.*ALPHAS_POP);
    PRECESS_LAMBDA_POP_COUNTER=find((alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*COUNTER_PASSING_POP.*ALPHAS_POP);
    lambda_Dpphi_value(bin,Ebin)=mean(abs(delta_pphi(PRECESS_LAMBDA_POP_COUNTER)))/mean(abs(delta_pphi(PRECESS_LAMBDA_POP_CO)));
end
PRECESS_LAMBDA_POP_COUNTER=find(COUNTER_PASSING_POP);
PRECESS_LAMBDA_POP_CO=find(CO_PASSING_POP);
Dpphi_Ekin_counter(Ebin)=mean(abs(delta_pphi(PRECESS_LAMBDA_POP_COUNTER)));
Dpphi_Ekin_co(Ebin)=mean(abs(delta_pphi(PRECESS_LAMBDA_POP_CO)));



Ebin=2
load('initialG_400keV_flat_pre_collapse.mat')
load('initialG_400keV_flat_precession_stats.mat')
ORBIT_POP=evalin('base',ORBIT_POP_name);
pphi_ini=alphas_pphi0;
lambda_ini=Bavg*alphas_mm./alphas_Ekin;
for bin=1:(length(precession_lambda_bins)-1)
    PRECESS_LAMBDA_POP_CO=find((alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*CO_PASSING_POP);
    PRECESS_LAMBDA_POP_COUNTER=find((alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*COUNTER_PASSING_POP);
    lambda_g_value(bin,Ebin)=(mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER))/mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)))));
    lambda_g_value(bin,Ebin)=lambda_g_value(bin,Ebin)*(mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_CO))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)))));
    lambda_dr_value(bin,Ebin)=mean(delta_r_avg(PRECESS_LAMBDA_POP_COUNTER))/mean(delta_r_avg(PRECESS_LAMBDA_POP_CO));
end
PRECESS_LAMBDA_POP_COUNTER=find(COUNTER_PASSING_POP);
PRECESS_LAMBDA_POP_CO=find(CO_PASSING_POP);
g_Ekin_counter(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)));
g_Ekin_co(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_CO)));
dr_Ekin_counter(Ebin)=mean(delta_r_avg(PRECESS_LAMBDA_POP_COUNTER));
dr_Ekin_co(Ebin)=mean(delta_r_avg(PRECESS_LAMBDA_POP_CO));
omega_vD_Ekin_counter(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)));
omega_vD_Ekin_co(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)));
load('flat400keV_collapse_Glisa_fc1h2_G250713.mat')
lambda_end=Bavg*alphas_mm./alphas_Ekin;
pphi_end=alphas_pphi0;
delta_pphi=pphi_end-pphi_ini;
delta_lambda=lambda_end-lambda_ini;
ALPHAS_POP=(~alphas_ejected);
for bin=1:(length(precession_lambda_bins)-1)
    PRECESS_LAMBDA_POP_CO=find((alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*CO_PASSING_POP.*ALPHAS_POP);
    PRECESS_LAMBDA_POP_COUNTER=find((alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*COUNTER_PASSING_POP.*ALPHAS_POP);
    lambda_Dpphi_value(bin,Ebin)=mean(abs(delta_pphi(PRECESS_LAMBDA_POP_COUNTER)))/mean(abs(delta_pphi(PRECESS_LAMBDA_POP_CO)));
end
PRECESS_LAMBDA_POP_COUNTER=find(COUNTER_PASSING_POP);
PRECESS_LAMBDA_POP_CO=find(CO_PASSING_POP);
Dpphi_Ekin_counter(Ebin)=mean(abs(delta_pphi(PRECESS_LAMBDA_POP_COUNTER)));
Dpphi_Ekin_co(Ebin)=mean(abs(delta_pphi(PRECESS_LAMBDA_POP_CO)));



Ebin=3

load('initialG_800keV_flat_pre_collapse.mat')
load('initialG_800keV_flat_precession_stats.mat')
ORBIT_POP=evalin('base',ORBIT_POP_name);
pphi_ini=alphas_pphi0;
lambda_ini=Bavg*alphas_mm./alphas_Ekin;
for bin=1:(length(precession_lambda_bins)-1)
    PRECESS_LAMBDA_POP_CO=find((alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*CO_PASSING_POP);
    PRECESS_LAMBDA_POP_COUNTER=find((alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*COUNTER_PASSING_POP);
    lambda_g_value(bin,Ebin)=(mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER))/mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)))));
    lambda_g_value(bin,Ebin)=lambda_g_value(bin,Ebin)*(mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_CO))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)))));
    lambda_dr_value(bin,Ebin)=mean(delta_r_avg(PRECESS_LAMBDA_POP_COUNTER))/mean(delta_r_avg(PRECESS_LAMBDA_POP_CO));
end
PRECESS_LAMBDA_POP_COUNTER=find(COUNTER_PASSING_POP);
PRECESS_LAMBDA_POP_CO=find(CO_PASSING_POP);
g_Ekin_counter(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)));
g_Ekin_co(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_CO)));
dr_Ekin_counter(Ebin)=mean(delta_r_avg(PRECESS_LAMBDA_POP_COUNTER));
dr_Ekin_co(Ebin)=mean(delta_r_avg(PRECESS_LAMBDA_POP_CO));
omega_vD_Ekin_counter(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)));
omega_vD_Ekin_co(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)));
load('flat800keV_collapse_Glisa_fc1h2_G260713.mat')
lambda_end=Bavg*alphas_mm./alphas_Ekin;
pphi_end=alphas_pphi0;
delta_pphi=pphi_end-pphi_ini;
delta_lambda=lambda_end-lambda_ini;
ALPHAS_POP=(~alphas_ejected);
for bin=1:(length(precession_lambda_bins)-1)
    PRECESS_LAMBDA_POP_CO=find((alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*CO_PASSING_POP.*ALPHAS_POP);
    PRECESS_LAMBDA_POP_COUNTER=find((alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*COUNTER_PASSING_POP.*ALPHAS_POP);
    lambda_Dpphi_value(bin,Ebin)=mean(abs(delta_pphi(PRECESS_LAMBDA_POP_COUNTER)))/mean(abs(delta_pphi(PRECESS_LAMBDA_POP_CO)));
end
PRECESS_LAMBDA_POP_COUNTER=find(COUNTER_PASSING_POP);
PRECESS_LAMBDA_POP_CO=find(CO_PASSING_POP);
Dpphi_Ekin_counter(Ebin)=mean(abs(delta_pphi(PRECESS_LAMBDA_POP_COUNTER)));
Dpphi_Ekin_co(Ebin)=mean(abs(delta_pphi(PRECESS_LAMBDA_POP_CO)));


Ebin=4

load('initialG_1600keV_flat_pre_collapse.mat')
load('initialG_1600keV_flat_precession_stats.mat')
ORBIT_POP=evalin('base',ORBIT_POP_name);
pphi_ini=alphas_pphi0;
lambda_ini=Bavg*alphas_mm./alphas_Ekin;
for bin=1:(length(precession_lambda_bins)-1)
    PRECESS_LAMBDA_POP_CO=find((alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*CO_PASSING_POP);
    PRECESS_LAMBDA_POP_COUNTER=find((alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*COUNTER_PASSING_POP);
    lambda_g_value(bin,Ebin)=(mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER))/mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)))));
    lambda_g_value(bin,Ebin)=lambda_g_value(bin,Ebin)*(mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_CO))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)))));
    lambda_dr_value(bin,Ebin)=mean(delta_r_avg(PRECESS_LAMBDA_POP_COUNTER))/mean(delta_r_avg(PRECESS_LAMBDA_POP_CO));
end
PRECESS_LAMBDA_POP_COUNTER=find(COUNTER_PASSING_POP);
PRECESS_LAMBDA_POP_CO=find(CO_PASSING_POP);
g_Ekin_counter(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)));
g_Ekin_co(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_CO)));
dr_Ekin_counter(Ebin)=mean(delta_r_avg(PRECESS_LAMBDA_POP_COUNTER));
dr_Ekin_co(Ebin)=mean(delta_r_avg(PRECESS_LAMBDA_POP_CO));
omega_vD_Ekin_counter(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)));
omega_vD_Ekin_co(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)));
load('flat1600keV_collapse_Glisa_fc1h2_G260713.mat')
lambda_end=Bavg*alphas_mm./alphas_Ekin;
pphi_end=alphas_pphi0;
delta_pphi=pphi_end-pphi_ini;
delta_lambda=lambda_end-lambda_ini;
ALPHAS_POP=(~alphas_ejected);
for bin=1:(length(precession_lambda_bins)-1)
    PRECESS_LAMBDA_POP_CO=find((alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*CO_PASSING_POP.*ALPHAS_POP);
    PRECESS_LAMBDA_POP_COUNTER=find((alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*COUNTER_PASSING_POP.*ALPHAS_POP);
    lambda_Dpphi_value(bin,Ebin)=mean(abs(delta_pphi(PRECESS_LAMBDA_POP_COUNTER)))/mean(abs(delta_pphi(PRECESS_LAMBDA_POP_CO)));
end
PRECESS_LAMBDA_POP_COUNTER=find(COUNTER_PASSING_POP);
PRECESS_LAMBDA_POP_CO=find(CO_PASSING_POP);
Dpphi_Ekin_counter(Ebin)=mean(abs(delta_pphi(PRECESS_LAMBDA_POP_COUNTER)));
Dpphi_Ekin_co(Ebin)=mean(abs(delta_pphi(PRECESS_LAMBDA_POP_CO)));


Ebin=5

load('initialG_2800keV_flat_pre_collapse.mat')
load('initialG_2800keV_flat_precession_stats.mat')
ORBIT_POP=evalin('base',ORBIT_POP_name);
pphi_ini=alphas_pphi0;
lambda_ini=Bavg*alphas_mm./alphas_Ekin;
for bin=1:(length(precession_lambda_bins)-1)
    PRECESS_LAMBDA_POP_CO=find((alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*CO_PASSING_POP);
    PRECESS_LAMBDA_POP_COUNTER=find((alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*COUNTER_PASSING_POP);
    lambda_g_value(bin,Ebin)=(mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER))/mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)))));
    lambda_g_value(bin,Ebin)=lambda_g_value(bin,Ebin)*(mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_CO))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)))));
    lambda_dr_value(bin,Ebin)=mean(delta_r_avg(PRECESS_LAMBDA_POP_COUNTER))/mean(delta_r_avg(PRECESS_LAMBDA_POP_CO));
end
PRECESS_LAMBDA_POP_COUNTER=find(COUNTER_PASSING_POP);
PRECESS_LAMBDA_POP_CO=find(CO_PASSING_POP);
g_Ekin_counter(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)))
g_Ekin_co(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_CO)))
dr_Ekin_counter(Ebin)=mean(delta_r_avg(PRECESS_LAMBDA_POP_COUNTER))
dr_Ekin_co(Ebin)=mean(delta_r_avg(PRECESS_LAMBDA_POP_CO))
omega_vD_Ekin_counter(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)));
omega_vD_Ekin_co(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)));
load('flat2800keV_collapse_Glisa_fc1h2_G230713.mat')
lambda_end=Bavg*alphas_mm./alphas_Ekin;
pphi_end=alphas_pphi0;
delta_pphi=pphi_end-pphi_ini;
delta_lambda=lambda_end-lambda_ini;
ALPHAS_POP=(~alphas_ejected);
for bin=1:(length(precession_lambda_bins)-1)
    PRECESS_LAMBDA_POP_CO=find((alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*CO_PASSING_POP.*ALPHAS_POP);
    PRECESS_LAMBDA_POP_COUNTER=find((alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*COUNTER_PASSING_POP.*ALPHAS_POP);
    lambda_Dpphi_value(bin,Ebin)=mean(abs(delta_pphi(PRECESS_LAMBDA_POP_COUNTER)))/mean(abs(delta_pphi(PRECESS_LAMBDA_POP_CO)));
end
PRECESS_LAMBDA_POP_COUNTER=find(COUNTER_PASSING_POP);
PRECESS_LAMBDA_POP_CO=find(CO_PASSING_POP);
Dpphi_Ekin_counter(Ebin)=mean(abs(delta_pphi(PRECESS_LAMBDA_POP_COUNTER)))
Dpphi_Ekin_co(Ebin)=mean(abs(delta_pphi(PRECESS_LAMBDA_POP_CO)))


figure(1)
subplot(3,1,1)
set(gca,'FontSize',26);
% title('Trapped particles');

hold on
grid on


% xlabel('\lambda_{0}')
ylabel('|\Deltap_\phi|_{coutner}/|\Deltap_\phi|_{co}');


plot(lambda_values,lambda_Dpphi_value(:,1),'g-+','LineWidth',2);
% plot(lambda_values,lambda_Dpphi_value(:,2),'k--','LineWidth',3);
plot(lambda_values,lambda_Dpphi_value(:,3),'b-.','LineWidth',4);
plot(lambda_values,lambda_Dpphi_value(:,4),'r','LineWidth',5);
plot(lambda_values,lambda_Dpphi_value(:,5),'y--','LineWidth',6);
plot(lambda_values,lambda_values.*0+1,'k','LineWidth',2);

% xlim([0.83 1.17])
% legend('\sim 65 keV','400 keV','800 keV','1600 keV','2800 keV')

xlim([0.2 0.8])


figure(1)
subplot(3,1,2)
set(gca,'FontSize',26);
% title('Trapped particles');

hold on
grid on


xlabel('\lambda_{0}')
ylabel('g_{coutner}/g_{co}');


plot(lambda_values,lambda_g_value(:,1),'g-+','LineWidth',2);
% plot(lambda_values,lambda_g_value(:,2),'k--','LineWidth',3);
plot(lambda_values,lambda_g_value(:,3),'b-.','LineWidth',4);
plot(lambda_values,lambda_g_value(:,4),'r','LineWidth',5);
plot(lambda_values,lambda_g_value(:,5),'y--','LineWidth',6);
plot(lambda_values,lambda_values.*0+1,'k','LineWidth',2);
% plot(lambda_values,lambda_values.*0+0.4*rmix,'color',[0.1 0.5 0.1],'LineWidth',2);

xlim([0.2 0.8])
ylim([0.5 2.5])

subplot(3,1,3)
set(gca,'FontSize',26);
% title('Trapped particles');

hold on
grid on


% xlabel('\lambda_{0}')
ylabel('\deltar_{counter}/\deltar_{co}');


plot(lambda_values,lambda_dr_value(:,1),'g-+','LineWidth',2);
% plot(lambda_values,lambda_dr_value(:,2),'k--','LineWidth',3);
plot(lambda_values,lambda_dr_value(:,3),'b-.','LineWidth',4);
plot(lambda_values,lambda_dr_value(:,4),'r','LineWidth',5);
plot(lambda_values,lambda_dr_value(:,5),'y--','LineWidth',6);
plot(lambda_values,lambda_values.*0+1,'k','LineWidth',2);

% xlim([0.83 1.17])
% legend('\sim 65 keV','400 keV','800 keV','1600 keV','2800 keV')

xlim([0.2 0.8])
ylim([0.95 max(lambda_dr_value(:,1))])

% legend('\sim 65 keV','400 keV','800 keV','1600 keV','2800 keV')
legend('\sim 65 keV','800 keV','1600 keV','2800 keV')



figure(2)
subplot(3,1,1)
set(gca,'FontSize',26);
hold on
grid on
plot([65 400 800 1600 2800],Dpphi_Ekin_counter,'b-.','LineWidth',4)
plot([65 400 800 1600 2800],Dpphi_Ekin_co,'r','LineWidth',4)
% xlabel('Ekin');
yl=ylabel('$$|\Delta p_\varphi|$$','interpreter','latex');
set(yl,'interpreter','latex')

figure(2)
subplot(3,1,2)
set(gca,'FontSize',26);
hold on
grid on
plot([65 400 800 1600 2800],dr_Ekin_counter,'b-.','LineWidth',4)
plot([65 400 800 1600 2800],dr_Ekin_co,'r','LineWidth',4)
% xlabel('Ekin');
ylabel('\deltar (m)');

figure(2)
subplot(3,1,3)
set(gca,'FontSize',26);
hold on
grid on
plot([65 400 800 1600 2800],g_Ekin_counter,'b-.','LineWidth',4)
plot([65 400 800 1600 2800],g_Ekin_co,'r','LineWidth',4)
xlabel('Ekin (keV)');
ylabel('\omega_{v_D}/\omega_{\psi}');




figure(3)
set(gca,'FontSize',26);
hold on
grid on
plot([65 400 800 1600 2800],omega_vD_Ekin_counter,'b-.','LineWidth',4)
plot([65 400 800 1600 2800],omega_vD_Ekin_co,'r','LineWidth',4)
xlabel('Ekin (keV)');
ylabel('|\omega_{v_D}|');
