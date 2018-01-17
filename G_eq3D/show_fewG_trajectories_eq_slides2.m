% time_stamp_end = 1231
% find(Ekin_outputG(:,1)<60*1e3)


% % 
P1=1;
P2=4;


% vperp_sq_outputG=vDsq_outputG;
% Bfield_part_outputG=interp2(scale_X,scale_Z,Btot_XZ_map',Xpos_outputG,Zpos_outputG,'*linear');
mm_outputG=Eperp_outputG./Bfield_outputG;
lambda_outputG=Bavg*mm_outputG./Ekin_outputG;
% psi_value_outputG=interp2(scale_X,scale_Z,psi_XZsmall_map',Xpos_outputG,Zpos_outputG,'*linear');
% lap_psi_outputG=interp2(scale_X,scale_Z,lap_psi_XZsmall_map',Xpos_outputG,Zpos_outputG,'*linear');

% for (n=1:24)
% psi_value_corr_outputG(n,:)=psipos_outputG(n,:)-alphas_psi_star_ini(n)+psi_star_outputG(n,:);
% end
lambda_P1=lambda_outputG(P1,1);
lambda_P2=lambda_outputG(P2,1);

Ekin_P1=round(Ekin_outputG(P1,1));
Ekin_P2=round(Ekin_outputG(P2,1));

lP1=strcat('\lambda _{ini} = ',num2str(lambda_P1));
lP2=strcat('\lambda _{ini} = ',num2str(lambda_P2));

EP1=strcat('Ekin _{ini} = ',num2str(Ekin_P1));
EP2=strcat('Ekin _{ini} = ',num2str(Ekin_P2));



phipos_outputG_wrap=wrap2pi(phipos_outputG);
omega_outputG=theta_outputG-phipos_outputG_wrap;
omega_outputG=wrap2pi(omega_outputG);

size_phi=257;
scale_tor=40*pi*((0:size_phi-1)/(size_phi-1));
Btot_X_vector=mean(Btot_XZ_map,2);
[Btot_Xphi_map tmp]=meshgrid(Btot_X_vector,ones(size_phi,1));
Btot_Xphi_map=Btot_Xphi_map';

end_time=size(time_scale_G,2);
ejection_time=size(time_scale_G,2);
% end_time=250;
%end_time=150000;
time_end=time_scale_G(end_time);
time_stamp_end=end_time

%time0=round(0.75*end_time);
time0=1;


close all
time_orbit=10000

figure(1);
set(gca,'FontSize',16);
hold on;

plot(Xpos_outputG(P1,time0:time_orbit),Zpos_outputG(P1,time0:time_orbit),'b','LineWidth',2);
plot(Xpos_outputG(P2,time0:time_orbit),Zpos_outputG(P2,time0:time_orbit),'r','LineWidth',2);
legend(lP1,lP2);

contour(scale_X,scale_Z,(psi_norm_XZsmall_map'-257)*psi_global/257,14,'y');
xlim([-0.4 0.7]);
ylim([-0.7 0.7]);

contour(scale_X,scale_Z,Btot_XZ_map',(2:0.2:4.4));axis xy square
colormap('summer')
colorbar;

plot(Xpos_outputG(P1,time0:time_orbit),Zpos_outputG(P1,time0:time_orbit),'b','LineWidth',2);
plot(Xpos_outputG(P2,time0:time_orbit),Zpos_outputG(P2,time0:time_orbit),'r','LineWidth',2);

xlabel('X (m)');
ylabel('Z (m)');






figure(3);

subplot(2,1,1);
set(gca,'FontSize',16);
grid on;
hold on;
xlim([time_scale_G(1) time_scale_G(end-1)])
plot(time_scale_G(time0:end_time),phipos_outputG(P1,time0:end_time),'b','LineWidth',2);
plot(time_scale_G(time0:end_time),phipos_outputG(P2,time0:end_time),'r','LineWidth',2);
xlabel('time (s)');
ylabel('\phi (rad)');
xlim([0 time_end]);
legend(lP1,lP2);

subplot(2,1,2);
set(gca,'FontSize',16);
grid on;
hold on
xlim([time_scale_G(1) time_scale_G(end-1)])
plot(time_scale_G(time0:end_time),vparallel_outputG(P1,time0:end_time),'b','LineWidth',2);
plot(time_scale_G(time0:end_time),vparallel_outputG(P2,time0:end_time),'r','LineWidth',2);
xlabel('time (s)');
ylabel('v_{||} (m/s)');
xlim([0 time_end]);









figure(5);
set(gca,'FontSize',16);

grid on;
hold on;
xlim([time_scale_G(1) time_scale_G(end-1)])
plot(time_scale_G(1:end_time),(mHe/eV)*(Xpos_outputG(P1,1:end_time)+R0).*(vphi_outputG(P1,1:end_time))-ZHe*(psi_value_outputG(P1,1:end_time)),'b','LineWidth',2);
plot(time_scale_G(1:end_time),pphi_outputG(P1,1:end_time),'k--','LineWidth',2);
plot(time_scale_G(1:end_time),(mHe/eV)*(Xpos_outputG(P2,1:end_time)+R0).*(vphi_outputG(P2,1:end_time))-ZHe*(psi_value_outputG(P2,1:end_time)),'r','LineWidth',2);
plot(time_scale_G(1:end_time),pphi_outputG(P2,1:end_time),'k--','LineWidth',2);
xlabel('time (s)');
ylabel('p_\phi');
xlim([0 time_end]);





figure(6)

subplot(2,1,1);
set(gca,'FontSize',16);

grid on;
hold on
xlim([time_scale_G(1) time_scale_G(end-1)])
% plot(time_scale_G(1:end_time),Eperp_outputG(P1,1:end_time),'r');
plot(time_scale_G(1:end_time),Eperp_outputG(P1,1:end_time),'b','LineWidth',2);
plot(time_scale_G(1:end_time),Eperp_outputG(P2,1:end_time),'r','LineWidth',2);
xlabel('time (s)');
ylabel('Eperp (eV)');
legend(lP1,lP2);

subplot(2,1,2);
set(gca,'FontSize',16);

grid on;
hold on
xlim([time_scale_G(1) time_scale_G(end-1)])
% plot(time_scale_G(1:end_time),vparallel_outputG(P1,1:end_time),'r');
plot(time_scale_G(1:end_time),vparallel_outputG(P1,1:end_time),'b','LineWidth',2);
plot(time_scale_G(1:end_time),vparallel_outputG(P2,1:end_time),'r','LineWidth',2);
xlabel('time (s)');
ylabel('v_{||} (m/s)');





figure(7)

subplot(2,1,1);
set(gca,'FontSize',16);
grid on;
hold on

plot(time_scale_G(time0:end_time),mm_outputG(P1,time0:end_time),'b','LineWidth',1.5);
plot(time_scale_G(time0:end_time),mm_outputG(P2,time0:end_time),'r','LineWidth',1.5);
xlabel('time (s)');
ylabel('\mu (eV T^{-1})');
% xlim([time_scale_G(1) time_scale_G(end-1)])
xlim([0 time_end]);
legend(lP1,lP2);

subplot(2,1,2);
set(gca,'FontSize',16);
grid on;
hold on
plot(time_scale_G(time0:end_time),rhoL_outputG(P1,time0:end_time),'b','LineWidth',1.5);
plot(time_scale_G(time0:end_time),rhoL_outputG(P2,time0:end_time),'r','LineWidth',1.5);
xlabel('time (s)')
ylabel('\rho_L (m)');
% xlim([time_scale_G(1) time_scale_G(end-1)])
xlim([0 time_end]);



figure(8)

subplot(2,1,1);
set(gca,'FontSize',16);
grid on;
hold on

plot(time_scale_G(time0:end_time),omega_outputG(P1,time0:end_time),'b','LineWidth',1.5);
plot(time_scale_G(time0:end_time),omega_outputG(P2,time0:end_time),'r','LineWidth',1.5);
xlabel('time (s)');
ylabel('\omega (rad)');
% xlim([time_scale_G(1) time_scale_G(end-1)])
xlim([0 time_end]);

subplot(2,1,2);
set(gca,'FontSize',16);
grid on;
hold on
plot(time_scale_G(time0:end_time),rhoL_outputG(P1,time0:end_time),'b','LineWidth',1.5);
plot(time_scale_G(time0:end_time),rhoL_outputG(P2,time0:end_time),'r','LineWidth',1.5);
xlabel('time (s)')
ylabel('\rho_L (m)');
% xlim([time_scale_G(1) time_scale_G(end-1)])
xlim([0 time_end]);

legend(lP1,lP2);



figure(9)

subplot(2,1,1);
set(gca,'FontSize',16);
grid on;
hold on

plot(time_scale_G(time0:end_time),psi_value_outputG(P1,time0:end_time),'b','LineWidth',1.5);
plot(time_scale_G(time0:end_time),psi_value_outputG(P2,time0:end_time),'r','LineWidth',1.5);
xlabel('time (s)');
ylabel('\psi (T.m^{-2})');
% xlim([time_scale_G(1) time_scale_G(end-1)])
xlim([0 time_end]);

subplot(2,1,2);
set(gca,'FontSize',16);
grid on;
hold on
plot(time_scale_G(time0:end_time),rhoL_outputG(P1,time0:end_time),'b','LineWidth',1.5);
plot(time_scale_G(time0:end_time),rhoL_outputG(P2,time0:end_time),'r','LineWidth',1.5);
xlabel('time (s)')
ylabel('\rho_L (m)');
% xlim([time_scale_G(1) time_scale_G(end-1)])
xlim([0 time_end]);




figure(10)

subplot(2,1,1);
set(gca,'FontSize',16);
grid on;
hold on

% plot(time_scale_G(time0:end_time),0.5*(mHe/eV)*vparallel_outputG(P1,time0:end_time).^2+Eperp_outputG(P1,time0:end_time),'b','LineWidth',1.5);
% plot(time_scale_G(time0:end_time),0.5*(mHe/eV)*vparallel_outputG(P2,time0:end_time).^2+Eperp_outputG(P2,time0:end_time),'r','LineWidth',1.5);
plot(time_scale_G(time0:end_time),Ekin_outputG(P1,time0:end_time),'b','LineWidth',1.5);
plot(time_scale_G(time0:end_time),Ekin_outputG(P2,time0:end_time),'r','LineWidth',1.5);
xlabel('time (s)');
ylabel('Ekin (eV)');
% xlim([time_scale_G(1) time_scale_G(end-1)])
xlim([0 time_end]);

subplot(2,1,2);
set(gca,'FontSize',16);
grid on;
hold on
plot(time_scale_G(time0:end_time),rhoL_outputG(P1,time0:end_time),'b','LineWidth',1.5);
plot(time_scale_G(time0:end_time),rhoL_outputG(P2,time0:end_time),'r','LineWidth',1.5);
xlabel('time (s)')
ylabel('\rho_L (m)');
% xlim([time_scale_G(1) time_scale_G(end-1)])
xlim([0 time_end]);

legend(lP1,lP2);




