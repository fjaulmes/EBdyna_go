load('quiver_data.mat')
load('deltaW_NBI_data.mat')

for n=1:128
   imagesc(scale_X,scale_Z,squeeze(deltaW_XZ_map_phi(n,:,:))',[-4.4 2.2]*1e5);
   title(strcat('\phi = ',num2str((n-1)*2*pi/128)))
   colorbar
   pause(0.05)
end

close all
deltaW_NBI_XZ_map=squeeze(mean(deltaW_XZ_map_phi(:,:,:),1));

imagesc(scale_X,scale_Z,deltaW_NBI_XZ_map',[-2.8 1.0]*1e5);


%%
figure(1)
set(gca,'fontsize',20)

VECSP=20

hold on
imagesc(scale_X,scale_Z,-(mu0*(R0^3)/(6*pi^2*Bphi0^2*r_value_q1_mean^4))*deltaW_NBI_XZ_map',[-5 5]*0.1);
contour(scale_X,scale_Z,-deltaW_NBI_XZ_map',[0 0],'g','linewidth',5)
colormap summer
% colormap('summer')
%
quiver(scXX(1:VECSP:end,1:VECSP:end-VECSP),scZZ(1:VECSP:end,1:VECSP:end-VECSP),gradP_X(1:VECSP:end-VECSP,1:VECSP:end)',gradP_Z(1:VECSP:end-VECSP,1:VECSP:end)','k')
quiver(scXX(1:VECSP:end,1:VECSP:end-VECSP),scZZ(1:VECSP:end,1:VECSP:end-VECSP),kappa_X_XZ_map(1:VECSP:end-VECSP,1:VECSP:end)',kappa_Z_XZ_map(1:VECSP:end-VECSP,1:VECSP:end)',0.5,'color',[0.3 0.3 1.0])

colorbar

xlim([-0.15 0.26])
ylim([-0.25 0.25])
xlabel('X (m)')
ylabel('Z (m)')


%%
figure(2)
plot((0:127)*2*pi/128,(mu0*(R0^3)/(6*pi^2*Bphi0^2*r_value_q1_mean^4))*deltaW_est_phi);