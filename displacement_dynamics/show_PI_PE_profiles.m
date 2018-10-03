reset_data_analysis_environment

TE_profile_radial_ini=interp1((0:Nradial-1)/(Nradial-1),TE_profile_interp_ini,rho_tor_scale);
TI_profile_radial_ini=interp1((0:Nradial-1)/(Nradial-1),TI_profile_interp_ini,rho_tor_scale);
NE_profile_radial_ini=interp1((0:Nradial-1)/(Nradial-1),NE_profile_interp_ini,rho_tor_scale);
NI_profile_radial_ini=interp1((0:Nradial-1)/(Nradial-1),NI_profile_interp_ini,rho_tor_scale);
PTOT_profile_radial_ini=interp1((0:Nradial-1)/(Nradial-1),PTOT_profile_interp_ini,rho_tor_scale);

figure(1)
subplot(2,1,1)
set(gca,'fontsize',20)
hold on
grid on
plot(radial_r_value_flux, NI_profile_radial_ini,'r','linewidth',3)
plot(radial_r_value_flux, NE_profile_radial_ini,'b--','linewidth',3)
xlabel('r (m)')
ylabel('m^{-3}')
legend('n_i','n_e')


subplot(2,1,2)
set(gca,'fontsize',20)
hold on
grid on
plot(radial_r_value_flux, NI_profile_radial_ini.*TI_profile_radial_ini*eV,'r','linewidth',3)
plot(radial_r_value_flux, NE_profile_radial_ini.*TE_profile_radial_ini*eV,'b--','linewidth',3)
xlabel('r (m)')
ylabel('Pa')
legend('P_i','P_e')