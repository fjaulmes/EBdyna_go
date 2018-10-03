    vperp0=sqrt(max(2*w0/mD-vpll_range.^2,0));
    vperp0_inf=sqrt(max(2*w0_inf/mD-vpll_range.^2,0));
    vperp0_sup=sqrt(max(2*w0_sup/mD-vpll_range.^2,0));

    volume_shell=(4/3)*pi*(theta0_sup)^3-(4/3)*pi*(theta0_inf)^3;
    %DV=vpll_bin_size*theta0./vperp0;
    volume_rings=vpll_bin_size.*(pi*(vperp0_sup).^2-pi*(vperp0_inf).^2);
    f_alpha_vpll=volume_rings/volume_shell;
