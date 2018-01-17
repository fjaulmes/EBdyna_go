
vD_X_map_theta=squeeze(mean(vD_X_PR_map_phi(1:end-1,:,:),1));
vD_Z_map_theta=squeeze(mean(vD_Z_PR_map_phi(1:end-1,:,:),1));
vD_phi_map_theta=squeeze(mean(vD_phi_PR_map_phi(1:end-1,:,:),1));

Efield_X_map_theta=squeeze(mean(Efield_X_map_phi(1:end-1,:,:),1));
Efield_Z_map_theta=squeeze(mean(Efield_Z_map_phi(1:end-1,:,:),1));
Efield_phi_map_theta=squeeze(mean(Efield_phi_map_phi(1:end-1,:,:),1));

vD_dot_E_map_phi=vD_X_PR_map_phi.*Efield_X_map_phi+vD_Z_PR_map_phi.*Efield_Z_map_phi+vD_phi_PR_map_phi.*Efield_phi_map_phi;
vD_dot_E_map_theta=squeeze(mean(vD_dot_E_map_phi(1:end-1,:,:),1));
% vD_dot_E_map_theta=squeeze(mean(vD_dot_E_map_phi(:,1:end-1,:),2));

dBtot_map_theta=squeeze(mean(dBtot_map_phi(1:end-1,:,:),1));
% dBtot_map_theta=squeeze(mean(dBtot_map_phi(:,1:end-1,:),2));

theta_scale=((1:NP)-1)*DTHETA;
% theta_scale=((1:129)-1)*DPHI;