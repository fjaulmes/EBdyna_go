close all
Delta_SP_recalc=2*delta_evol.*vA_out./ksi_dot_recalc(1:TRANSITION_FRAME);

figure(1)
subplot(3,1,1)
set(gca,'fontsize',20)

grid on
hold on;
% plot(ksi0_evol_lin(1:TRANSITION_FRAME),f1(1:TRANSITION_FRAME).*E1,'r','linewidth',2)
% plot(ksi0_evol_lin(1:TRANSITION_FRAME),f2(1:TRANSITION_FRAME).*E2,'b--','linewidth',2)
plot(ksi0_evol_lin(1:TRANSITION_FRAME),0.5*(E1+E2)-Efinal,'b--','linewidth',2)
% plot(ksi0_evol_lin(1:TRANSITION_FRAME),E1,'r','linewidth',2)
% plot(ksi0_evol_lin(1:TRANSITION_FRAME),E2,'b--','linewidth',2)
% plot(ksi0_evol_lin(1:TRANSITION_FRAME),Efinal,'k-.','linewidth',2)
plot(ksi0_evol_lin(1:TRANSITION_FRAME),joule_heating./dvolume3_recalc,'g','linewidth',2)

xlim([0.1 2.2]*1e-1)
ylim([0 2800])

ylabel('J / m^3')
xlabel('\xi (m)')


% hl=legend('$E_1$','$E_2$','$E_3$','joule heating');
hl=legend('$\Delta E_{mag}$','joule heating');
set(hl,'interpreter','latex')

subplot(3,1,2)
set(gca,'fontsize',20)

grid on
hold on;
plot(ksi0_evol_lin(1:TRANSITION_FRAME),delta_evol,'r','linewidth',2)
plot(ksi0_evol_lin(1:TRANSITION_FRAME),delta_evol*0+rholi,'k-.','linewidth',2)
plot(ksi0_evol_lin(1:TRANSITION_FRAME),Delta_SP_recalc/20,'b--','linewidth',2)

xlim([0.1 2.2]*1e-1)
ylim([0 0.02])

% ylabel('\delta (cm) \Delta (m)')
ylabel('m')
xlabel('\xi (m)')

hl=legend('$\delta$','$\rho_{Li}$','$\Delta/20$');
set(hl,'interpreter','latex')



subplot(3,1,3)
set(gca,'fontsize',20)
grid on
hold on;

plot(ksi0_evol_lin(1:TRANSITION_FRAME),0.05*vA_out,'r','linewidth',2)
plot(ksi0_evol_lin(1:TRANSITION_FRAME),ksi_dot_recalc(1:TRANSITION_FRAME),'b--','linewidth',2)
plot(ksi0_evol_lin(1:TRANSITION_FRAME),rx_dot_recalc_evol(1:TRANSITION_FRAME),'k-.','linewidth',2)

xlabel('\xi (m)')
ylabel('m/s')

xlim([0.1 2.2]*1e-1)
hl=legend('$v_{A*}/20$','$\dot{\xi}$','$ \dot{r_x}$');
set(hl,'interpreter','latex')

ylim([0 11000])
