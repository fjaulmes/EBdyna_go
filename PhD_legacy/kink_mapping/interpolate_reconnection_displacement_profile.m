for t=1:length(time_scale_lin)
    rx_evol(t)=interp1(radial_r_value_flux,psi_bar_rho_pol,rx_evol(t));
    xi_evol(t)=interp1(radial_r_value_flux,psi_bar_rho_pol,ksi0_evol(t));
end