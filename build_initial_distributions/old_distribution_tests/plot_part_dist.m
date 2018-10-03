figure(1);
hold on

psi_bin_sup-257;

plot(alphas_pos_x,alphas_pos_z,'b.')
plot(scale_X,Z_psi_fit_up_small(psi_bin_sup,:),'r');
plot(scale_X,Z_psi_fit_down_small(psi_bin_sup,:),'r');