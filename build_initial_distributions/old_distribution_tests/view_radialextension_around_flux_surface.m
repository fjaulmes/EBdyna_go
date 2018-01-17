
Eth=eV*2e6;
vth_perp=sqrt(2*(2/3)*Eth/mHe);
Bpol_tot_XZ_map=sqrt(BpolX_XZ_map.^2+BpolZ_XZ_map.^2);

wci_XZ_map=2*eV*Btot_XZ_map/mHe;
wci_pol_XZ_map=2*eV*Bpol_tot_XZ_map/mHe;
rhoL_XZ_map=vth_perp./wci_XZ_map;
rhoL_pol_XZ_map=vth_perp./wci_pol_XZ_map;
r_value_XZ_map=interp1(1:257,radial_r_value_flux,radial_XZsmall_map);
r_value_XZ_map=max(r_value_XZ_map,DX);



%extension_XZ_map=min(2*rhoL_XZ_map.*sqrt(R0./r_value_XZ_map),r_value_XZ_map);
%extension_XZ_map=2*rhoL_XZ_map.*sqrt(2*R0./r_value_XZ_map);
extension_XZ_map=rhoL_pol_XZ_map;
extension_XZ_map_lim=extension_XZ_map;

extension_XZ_map_lim=extension_XZ_map+(abs(extension_XZ_map-r_value_XZ_map)<0.005);

figure(1);
set(gca,'FontSize',22);
imagesc(scale_X,scale_Z,extension_XZ_map_lim',[0 1.2]);
xlabel('X (m)');
ylabel('Z (m)');
% xlim([-0.4 0.5]);
% ylim([-0.6 0.6]);
title('approximative orbit width')
h=colorbar;
xlabel(h,'m');
axis xy square
