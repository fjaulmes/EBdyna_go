Ni_XZsmall_map=interp1(1:257,Ne_profile,psi_norm_XZsmall_map);
vA_map=(Btot_XZ_map./sqrt(mu0*Ni_XZsmall_map*mDT));
omegaA_map=vA_map./(2*(1.3)*R0);
grid on; hold on
contour(scale_X,scale_Z,q_initial_XZsmall_map',((1:30)+0.5)/5,'k--');
contour(scale_X,scale_Z,omegaA_map',(6:1:18)*1e5,'LineWidth',3);
axis equal
colorbar
xlabel('X (m)')
ylabel('Z (m)')