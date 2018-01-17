reset_data_analysis_environment
r_max=max(radial_r_value_flux)
rho_scale=radial_r_value_flux/r_max;

figure(2);grid on;hold on;
set(gca,'FontSize',22);
plot(rho_scale,q_initial_profile,'b','LineWidth',3)
xlabel('\rho')
ylabel('q')
xlim([0 0.95]);
ylim([0.6 4.5]);
pause(0.1)

load('../data_tokamak/psi_star_evol.mat')


%%
figure(3);grid on;hold on;

set(gca,'FontSize',26);

tau_cr=atan(2)/pi;
rmix=rx_evol(end);

plot(time_scale_lin, rx_evol_lin,'b--','LineWidth',2);
plot(time_scale_lin, ksi0_evol_lin,'r','LineWidth',2);
plot(time_scale_lin, time_scale_lin.*0+rmix,'k-.','LineWidth',3);

%set( findobj(gcf,'Type','line'),'LineWidth',1.5);
ylabel('r (m)')
xlabel('t (\tau_{coll})')
xlim([0 1]);
ylim([0 0.4]);

legend('\rho_x','\xi','\rho_{mix}');