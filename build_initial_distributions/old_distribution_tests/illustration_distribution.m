close all

figure(1);
set(gca,'FontSize',22);

hist(alphas_Ekin,Evalues);
xlabel('Ekin (eV)');
ylabel('Npart')


figure(2);
set(gca,'FontSize',22);

hist(alphas_vpll,40);
xlabel('v_{||} (m/s)');
ylabel('Npart')

figure(3);
set(gca,'FontSize',22);
plot(alphas_pos_x,alphas_pos_z,'b.');
xlabel('X (m)');
ylabel('Z (m)');
title('spread of particles ')


figure(4);
set(gca,'FontSize',22);
axis xy;
imagesc(radial_bins,Evalues,f_alpha_energy_percentage_norm');
xlabel('\psi (a.u.)');
ylabel('Ekin (eV)');
title('distribution of particles (normalized)')
colorbar;

