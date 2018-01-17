function values=interp3_phi_theta_psi(scale_theta,scale_phi,scale_psi,PTP_map,thetas,phis,psis,PART_LIST)


values=interp3(scale_theta,scale_phi,scale_psi,PTP_map,thetas(PART_LIST),phis(PART_LIST),psis(PART_LIST),'*linear');

% values=lininterp3(scale_phi,scale_theta,scale_psi,PTP_map,phis(PART_LIST),thetas(PART_LIST),psis(PART_LIST));
