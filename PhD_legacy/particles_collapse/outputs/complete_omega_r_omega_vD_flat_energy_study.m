% reset_data_analysis_environment

close all;

rmix=radial_r_value_flux(size_r-4)

PRECESS_LAMBDA_BIN_SIZE=0.03;
precession_lambda_bins=(21*PRECESS_LAMBDA_BIN_SIZE:PRECESS_LAMBDA_BIN_SIZE:50*PRECESS_LAMBDA_BIN_SIZE);
lambda_values=(21.5*PRECESS_LAMBDA_BIN_SIZE:PRECESS_LAMBDA_BIN_SIZE:49.5*PRECESS_LAMBDA_BIN_SIZE);

lambda_omega_vD_value=zeros(length(precession_lambda_bins)-1,5);
lambda_omega_r_value=zeros(length(precession_lambda_bins)-1,5);
lambda_Dpphi_value=zeros(length(precession_lambda_bins)-1,5);

ORBIT_POP_name='ALL_TRAPPED_POP'

Ebin=1
load('initialG_alphas_lowEkin_precession_stats.mat')
load('initialG_alphas_lowEkin_pre_collapse.mat')
ORBIT_POP=evalin('base',ORBIT_POP_name);
pphi_ini=alphas_pphi0;
lambda_ini=Bavg*alphas_mm./alphas_Ekin;
for bin=1:(length(precession_lambda_bins)-1)
    PRECESS_LAMBDA_POP=find((alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*ALL_TRAPPED_POP);
    lambda_omega_vD_value(bin,Ebin)=mean((omega_precess_avg(PRECESS_LAMBDA_POP)));
    lambda_omega_r_value(bin,Ebin)=mean((omega_r_avg(PRECESS_LAMBDA_POP)))/mean((omega_phi_avg(PRECESS_LAMBDA_POP)));
end
load('alphas_lowEkin_collapse_Glisa_fc1h2_G250713.mat')
lambda_end=Bavg*alphas_mm./alphas_Ekin;
pphi_end=alphas_pphi0;
delta_pphi=pphi_end-pphi_ini;
delta_lambda=lambda_end-lambda_ini;
ALPHAS_POP=(~alphas_ejected);
for bin=1:(length(precession_lambda_bins)-1)
    PRECESS_LAMBDA_POP=find((alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*ALL_TRAPPED_POP.*ALPHAS_POP);
    lambda_Dpphi_value(bin,Ebin)=mean(abs(delta_pphi(PRECESS_LAMBDA_POP)));
end



Ebin=2
load('initialG_400keV_flat_pre_collapse.mat')
load('initialG_400keV_flat_precession_stats.mat')
ORBIT_POP=evalin('base',ORBIT_POP_name);
pphi_ini=alphas_pphi0;
lambda_ini=Bavg*alphas_mm./alphas_Ekin;
for bin=1:(length(precession_lambda_bins)-1)
    PRECESS_LAMBDA_POP=find((alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*ALL_TRAPPED_POP);
    lambda_omega_vD_value(bin,Ebin)=mean((omega_precess_avg(PRECESS_LAMBDA_POP)));
    lambda_omega_r_value(bin,Ebin)=mean((omega_r_avg(PRECESS_LAMBDA_POP)))/mean((omega_phi_avg(PRECESS_LAMBDA_POP)));
end
load('flat400keV_collapse_Glisa_fc1h2_G250713.mat')
lambda_end=Bavg*alphas_mm./alphas_Ekin;
pphi_end=alphas_pphi0;
delta_pphi=pphi_end-pphi_ini;
delta_lambda=lambda_end-lambda_ini;
ALPHAS_POP=(~alphas_ejected);
for bin=1:(length(precession_lambda_bins)-1)
    PRECESS_LAMBDA_POP=find((alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*ALL_TRAPPED_POP.*ALPHAS_POP);
    lambda_Dpphi_value(bin,Ebin)=mean(abs(delta_pphi(PRECESS_LAMBDA_POP)));
end



Ebin=3

load('initialG_800keV_flat_pre_collapse.mat')
load('initialG_800keV_flat_precession_stats.mat')
ORBIT_POP=evalin('base',ORBIT_POP_name);
pphi_ini=alphas_pphi0;
lambda_ini=Bavg*alphas_mm./alphas_Ekin;
for bin=1:(length(precession_lambda_bins)-1)
    PRECESS_LAMBDA_POP=find((alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*ALL_TRAPPED_POP);
    lambda_omega_vD_value(bin,Ebin)=mean((omega_precess_avg(PRECESS_LAMBDA_POP)));
    lambda_omega_r_value(bin,Ebin)=mean((omega_r_avg(PRECESS_LAMBDA_POP)))/mean((omega_phi_avg(PRECESS_LAMBDA_POP)));
end
load('flat800keV_collapse_Glisa_fc1h2_G260713.mat')
lambda_end=Bavg*alphas_mm./alphas_Ekin;
pphi_end=alphas_pphi0;
delta_pphi=pphi_end-pphi_ini;
delta_lambda=lambda_end-lambda_ini;
ALPHAS_POP=(~alphas_ejected);
for bin=1:(length(precession_lambda_bins)-1)
    PRECESS_LAMBDA_POP=find((alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*ALL_TRAPPED_POP.*ALPHAS_POP);
    lambda_Dpphi_value(bin,Ebin)=mean(abs(delta_pphi(PRECESS_LAMBDA_POP)));
end


Ebin=4

load('initialG_1600keV_flat_pre_collapse.mat')
load('initialG_1600keV_flat_precession_stats.mat')
ORBIT_POP=evalin('base',ORBIT_POP_name);
pphi_ini=alphas_pphi0;
lambda_ini=Bavg*alphas_mm./alphas_Ekin;
for bin=1:(length(precession_lambda_bins)-1)
    PRECESS_LAMBDA_POP=find((alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*ALL_TRAPPED_POP);
    lambda_omega_vD_value(bin,Ebin)=mean((omega_precess_avg(PRECESS_LAMBDA_POP)));
    lambda_omega_r_value(bin,Ebin)=mean((omega_r_avg(PRECESS_LAMBDA_POP)))/mean((omega_phi_avg(PRECESS_LAMBDA_POP)));
end
load('flat1600keV_collapse_Glisa_fc1h2_G260713.mat')
lambda_end=Bavg*alphas_mm./alphas_Ekin;
pphi_end=alphas_pphi0;
delta_pphi=pphi_end-pphi_ini;
delta_lambda=lambda_end-lambda_ini;
ALPHAS_POP=(~alphas_ejected);
for bin=1:(length(precession_lambda_bins)-1)
    PRECESS_LAMBDA_POP=find((alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*ALL_TRAPPED_POP.*ALPHAS_POP);
    lambda_Dpphi_value(bin,Ebin)=mean(abs(delta_pphi(PRECESS_LAMBDA_POP)));
end


Ebin=5

load('initialG_2800keV_flat_pre_collapse.mat')
load('initialG_2800keV_flat_precession_stats.mat')
ORBIT_POP=evalin('base',ORBIT_POP_name);
pphi_ini=alphas_pphi0;
lambda_ini=Bavg*alphas_mm./alphas_Ekin;
for bin=1:(length(precession_lambda_bins)-1)
    PRECESS_LAMBDA_POP=find((alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*ALL_TRAPPED_POP);
    lambda_omega_vD_value(bin,Ebin)=mean((omega_precess_avg(PRECESS_LAMBDA_POP)));
    lambda_omega_r_value(bin,Ebin)=mean((omega_r_avg(PRECESS_LAMBDA_POP)))/mean((omega_phi_avg(PRECESS_LAMBDA_POP)));
end
load('flat2800keV_collapse_Glisa_fc1h2_G230713.mat')
lambda_end=Bavg*alphas_mm./alphas_Ekin;
pphi_end=alphas_pphi0;
delta_pphi=pphi_end-pphi_ini;
delta_lambda=lambda_end-lambda_ini;
ALPHAS_POP=(~alphas_ejected);
for bin=1:(length(precession_lambda_bins)-1)
    PRECESS_LAMBDA_POP=find((alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*ALL_TRAPPED_POP.*ALPHAS_POP);
    lambda_Dpphi_value(bin,Ebin)=mean(abs(delta_pphi(PRECESS_LAMBDA_POP)));
end

lambda_omega_r_value=abs(lambda_omega_r_value);

figure(1)
subplot(3,1,1)
set(gca,'FontSize',26);
% title('Trapped particles');

hold on
grid on


% xlabel('\lambda_{0}')
ylabel('|\Deltap_\phi|');


plot(lambda_values,lambda_Dpphi_value(:,1),'g-+','LineWidth',2);
plot(lambda_values,lambda_Dpphi_value(:,2),'k--','LineWidth',2);
plot(lambda_values,lambda_Dpphi_value(:,3),'b-.','LineWidth',4);
plot(lambda_values,lambda_Dpphi_value(:,4),'r','LineWidth',5);
% plot(lambda_values,lambda_Dpphi_value(:,5),'y--','LineWidth',6);

xlim([0.83 1.15])
% legend('\sim 65 keV','400 keV','800 keV','1600 keV','2800 keV')

ylim([0.0 0.12])


figure(1)
subplot(3,1,2)
set(gca,'FontSize',26);

hold on
grid on
% xlabel('\lambda_{0}')
ylabel('\omega_{r} / \omega_\phi');

plot(lambda_values,lambda_omega_r_value(:,1),'g-+','LineWidth',2);
plot(lambda_values,lambda_omega_r_value(:,2),'k--','LineWidth',2);
plot(lambda_values,lambda_omega_r_value(:,3),'b-.','LineWidth',4);
plot(lambda_values,lambda_omega_r_value(:,4),'r','LineWidth',5);
% plot(lambda_values,lambda_omega_r_value(:,5),'y--','LineWidth',6);
xlim([0.83 1.15])
ylim([0.0 20])



figure(1)
subplot(3,1,3)
set(gca,'FontSize',26);
% title('Trapped particles');

hold on
grid on


xlabel('\lambda_{0}')
ylabel('\omega_{v_D} (m)');


plot(lambda_values,lambda_omega_vD_value(:,1),'g-+','LineWidth',2);
plot(lambda_values,lambda_omega_vD_value(:,2),'k--','LineWidth',2);
plot(lambda_values,lambda_omega_vD_value(:,3),'b-.','LineWidth',4);
plot(lambda_values,lambda_omega_vD_value(:,4),'r','LineWidth',5);
% plot(lambda_values,lambda_omega_vD_value(:,5),'y--','LineWidth',6);
% plot(lambda_values,lambda_values.*0+rmix,'color',[0.1 0.5 0.1],'LineWidth',2);
plot(lambda_values,lambda_values.*0-omega_crash,'color',[0.1 0.5 0.1],'LineWidth',3);
plot(lambda_values,lambda_values.*0+omega_crash,'color',[0.1 0.5 0.1],'LineWidth',3);

xlim([0.83 1.15])
legend('65 keV','400 keV','800 keV','1600 keV')

ylim([min(lambda_omega_vD_value(:,4)) max(lambda_omega_vD_value(:,4)) ])
