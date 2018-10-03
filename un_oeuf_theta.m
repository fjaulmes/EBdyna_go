close all
figure(1)
hold on; grid on
imagesc(scale_X,scale_Z,theta_XZsmall_map')
contour(scale_X,scale_Z,psi_XZsmall_map',[0 0],'k','linewidth',4)
contour(scale_X,scale_Z,psi_XZsmall_map',[psi_q1 psi_q1],'g','linewidth',1.5)
contour(scale_X,scale_Z,theta_XZsmall_map',[0 0],'k','linewidth',4)
contour(scale_X,scale_Z,theta_XZsmall_map',[0.5*pi 0.5*pi],'k','linewidth',4)
contour(scale_X,scale_Z,theta_XZsmall_map',[pi pi],'k','linewidth',4)
contour(scale_X,scale_Z,theta_XZsmall_map',[1.5*pi  1.5*pi],'k','linewidth',4)
set(gca,'FontSize',26)
xlabel('X (m)')
ylabel('Z (m)')
xlim([-0.56 0.56])
ylim([-0.86 0.86])
% colorbar
axis square
title('Poloidal angle \theta')

