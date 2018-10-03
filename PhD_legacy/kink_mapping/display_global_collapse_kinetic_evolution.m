reset_data_analysis_environment
load('../data_tokamak/psi_profiles_kadomtsev.mat')
load('../data_tokamak/Pkinetic_final_profiles.mat')
load('../data_tokamak/alphas_kinetic_density_profiles.mat')

psi_mix_pos=size_r-4
delta_psi_q1=15;
rho_scale=radial_r_value_flux/max(radial_r_value_flux);
rho_q1=interp1(psi_scale,rho_scale,psi_q1)
rho_mix=interp1(1:Nradial,rho_scale,psi_mix_pos)

figure(2);
subplot(3,1,1);
set(gca,'FontSize',16);
grid on;
hold on;
plot(rho_scale,q_initial,'b--','LineWidth',3);
% xlabel('r (m)');
ylabel('q');
%plot(pos_Bstar_final,q_final_interp,'r');
plot(rho_scale,q_final,'r','LineWidth',3);
plot(rho_scale,q_final_profile_diff,'--','Color',[1 0.6 0.3],'LineWidth',3);


plot([rho_q1 rho_q1],[-1 2],'g--','linewidth',3)
plot([rho_mix rho_mix],[-1 2],'g--','linewidth',3)

ylim([0.85 1.33]);
xlim([0 0.6])

hl=legend('before collapse','after collapse','after diffusion');
set(hl,'fontsize',14)

% subplot(4,1,2);
% set(gca,'FontSize',16);
% plot(rho_scale,psi_star_initial,'b--','LineWidth',3);
% grid on;
% hold on;
% plot(rho_scale,psi_star_final,'r','LineWidth',3);
% 
% plot([rho_q1 rho_q1],[-1 1],'g--','linewidth',3)
% plot([rho_mix rho_mix],[-1 1],'g--','linewidth',3)
% 
% % xlabel('r (m)');
% ylabel('\Psi_* (T.m^{-2})');
% % legend('\psi_{*-}','\psi_{*+}');
% hl=legend('before collapse','after collapse');
% set(hl,'fontsize',14)
% ylim([-0.015 0.2])
% xlim([0 0.6])


subplot(3,1,2);
set(gca,'FontSize',16);
plot(rho_scale,P_initial_profile/P0,'b--','LineWidth',3);
% plot(rho_scale,P_alphas_ini/P0,'g--','LineWidth',3);
grid on;
hold on;
P_final_kprofile(isnan(P_final_kprofile))=P_final_profile(isnan(P_final_kprofile));
% P_final_kprofile(size_r+8:end)=P_final_profile(size_r+8:end);
plot(rho_scale,P_final_kprofile/P0,'r','LineWidth',3);

plot([rho_q1 rho_q1],[-1 1],'g--','linewidth',3)
plot([rho_mix rho_mix],[-1 1],'g--','linewidth',3)

% plot(rho_scale,P_alphas_end/P0,'k','LineWidth',3);
xlim([0 r_value_q1_mean*2.0]);
xlabel('\rho');
ylabel('P/P_0');
hl=legend('before collapse','after collapse');
set(hl,'fontsize',14)
ylim([0 1.02])
xlim([0 0.6])

P_alphas_ini=n_alphas_initial_profile.*T_alphas_initial_profile;
P_alphas_end=n_alphas_final_profile.*T_alphas_final_profile;

subplot(3,1,3);
set(gca,'FontSize',16);
% plot(rho_scale,P_initial_profile/P0,'b--','LineWidth',3);
plot(rho_scale,P_alphas_ini/P0,'b--','LineWidth',3);
grid on;
hold on;
% plot(rho_scale,P_final_profile/P0,'r','LineWidth',3);
plot(rho_scale,P_alphas_end/P0,'r','LineWidth',3);
plot([rho_q1 rho_q1],[-1 1],'g--','linewidth',3)
plot([rho_mix rho_mix],[-1 1],'g--','linewidth',3)

xlim([0 r_value_q1_mean*2.0]);
xlabel('\rho');
ylabel('P_\alpha/P_0');
hl=legend('before collapse','after collapse');
set(hl,'fontsize',14)
ylim([0 0.15])
xlim([0 0.6])
% %%
% figure(8)
% set(gca,'FontSize',24);
% plot(rho_scale,psi_star_initial,'b--','LineWidth',3);
% grid on;
% hold on;
% plot(rho_scale,psi_star_final,'r','LineWidth',3);
% xlim([0 1]);
% % xlabel('r (m)');
% ylabel('\Psi_* (T.m^{-2})');
% legend('\psi_{*-}','\psi_{*+}');
% ylim([-0.05 0.05])
