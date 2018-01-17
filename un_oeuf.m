load('./data_tokamak/physics_constants.mat')
load('./data_tokamak/tokamak_PR_map.mat')
load('./data_tokamak/tokamak_map_dimensions.mat')
load('./data_tokamak/flux_geometry.mat')
load('./data_tokamak/pressure_profile.mat')
load('./data_tokamak/q_profile.mat')
rho_scale=radial_r_value_flux/max(radial_r_value_flux);

close all
figure(1)
hold on; grid on
contour(X_scale+R0,-Z_scale,psi_XZ_map',psi_scale(2:22:end),'k')
contour(X_scale+R0,-Z_scale,psi_XZ_map',[0 0],'k','linewidth',4)
set(gca,'FontSize',26)
xlabel('R (m)')
ylabel('Z (m)')
% colorbar
axis square
% title('Contour of poloidal flux surfaces of AUG')

figure(2)
% % pol_Te=[-8*1e3 -8.0*1e3 -4.0*1e3 -0.8*1e3 2.8*1e3 -8.0*1e2 5.0*1e3]
% pol_Ne=Ne0*[-6*1e4 -3*1e3 -3*1e3 -2*1e3 -1.0*1e3 -1.0*1e3 -0.8*1e3 1.8*1e3 -8.0*1e2 5.0*1e3]/(5.0*1e3)
% % pol_Te=fliplr(pol_Te);
% % Te_prof_r=polyval(pol_Te,radial_r_value_flux);
% Ne_prof_r=polyval(pol_Ne,radial_r_value_flux);
% Te_prof_r=0.5*P_initial_profile./Ne_prof_r;

subplot(3,1,1)
hold on; grid on

plot(rho_scale,Te_profile/(1e3*eV),'linewidth',4)
set(gca,'FontSize',25)
xlim([0 0.8])
ylabel('Te (keV)')
% xlabel('r (m)')
% plot(radial_r_value_flux,P_initial_profile,'linewidth',4)


subplot(3,1,2)
hold on; grid on

plot(rho_scale,Ne_profile,'linewidth',4)
set(gca,'FontSize',25)
xlim([0 0.8])
ylabel('Ne (m^{-3})')
% xlabel('r (m)')


subplot(3,1,3)
hold on; grid on

plot(rho_scale,P_initial_profile,'linewidth',4)
set(gca,'FontSize',25)
ylabel('P (Pa)')
xlim([0 0.8])
xlabel('\rho')



figure(3)

subplot(2,1,1)
hold on; grid on

plot(rho_scale,P_initial_profile,'linewidth',4)
set(gca,'FontSize',25)
ylabel('P (Pa)')
xlim([0 0.8])



subplot(2,1,2)
hold on; grid on

plot(rho_scale,q_initial_profile,'linewidth',4)
set(gca,'FontSize',25)
ylabel('q')
xlim([0 0.8])
xlabel('\rho')
