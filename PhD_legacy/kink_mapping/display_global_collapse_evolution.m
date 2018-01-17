rho_scale=radial_r_value_flux/max(radial_r_value_flux);


figure(1)
plot(rho_scale,P_initial_profile/P0,'b--','LineWidth',3);
% plot(rho_scale,P_alphas_ini/P0,'g--','LineWidth',3);
grid on;
hold on;
plot(rho_scale,P_final_profile/P0,'r','LineWidth',3);




figure(2);
subplot(4,1,1);
set(gca,'FontSize',16);
grid on;
hold on;
plot(rho_scale,q_initial,'b--','LineWidth',3);
% xlabel('r (m)');
ylabel('q');
%plot(pos_Bstar_final,q_final_interp,'r');
plot(rho_scale,q_final,'r','LineWidth',3);
plot(rho_scale,q_final_profile_diff,'--','Color',[1 0.6 0.3],'LineWidth',3);
ylim([0.85 1.2]);
xlim([0 0.6])

hl=legend('before collapse','after collapse','after diffusion');
set(hl,'fontsize',14)

subplot(4,1,2);
set(gca,'FontSize',16);
plot(rho_scale,psi_star_initial,'b--','LineWidth',3);
grid on;
hold on;
plot(rho_scale,psi_star_final,'r','LineWidth',3);
% xlabel('r (m)');
ylabel('\Psi_* (T.m^{-2})');
% legend('\psi_{*-}','\psi_{*+}');
hl=legend('before collapse','after collapse');
set(hl,'fontsize',14)
ylim([-0.005 0.08])
xlim([0 0.6])

% P_alphas_ini=n_alphas_initial_profile.*T_alphas_initial_profile;
% P_alphas_end=n_alphas_final_profile.*T_alphas_final_profile;

subplot(4,1,3);
set(gca,'FontSize',16);
plot(rho_scale,P_initial_profile/P0,'b--','LineWidth',3);
% plot(rho_scale,P_alphas_ini/P0,'g--','LineWidth',3);
grid on;
hold on;
plot(rho_scale,P_final_profile/P0,'r','LineWidth',3);
% plot(rho_scale,P_alphas_end/P0,'k','LineWidth',3);
xlim([0 r_value_q1_mean*2.0]);
xlabel('\rho');
ylabel('P/P_0');
hl=legend('before collapse','after collapse');
set(hl,'fontsize',14)
ylim([0 1.02])
xlim([0 0.6])


subplot(4,1,4);
set(gca,'FontSize',16);
% plot(rho_scale,P_initial_profile/P0,'b--','LineWidth',3);
plot(rho_scale,P_alphas_ini/P0,'b--','LineWidth',3);
grid on;
hold on;
% plot(rho_scale,P_final_profile/P0,'r','LineWidth',3);
plot(rho_scale,P_alphas_end/P0,'r','LineWidth',3);
xlim([0 r_value_q1_mean*2.0]);
xlabel('\rho');
ylabel('P_\alpha/P_0');
hl=legend('before collapse','after collapse');
set(hl,'fontsize',14)
ylim([0 0.12])
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
