% find(Ekin0_output(:,1)<60*1e3)



% % 
P1=1;
P2=2;
P5=9;
P6=10;
% 

% vperp_sq_output=vDsq_output;
% Bfield_part_output=interp2(scale_X,scale_Z,Btot_XZ_map',Xpos_output,Zpos_output,'*linear');
mm_output=Eperp_output./Bfield_output;
lambda_output=Bavg*mm_output./Ekin_output;
% psi_value_output=interp2(scale_X,scale_Z,psi_XZsmall_map',Xpos_output,Zpos_output,'*linear');
% lap_psi_output=interp2(scale_X,scale_Z,lap_psi_XZsmall_map',Xpos_output,Zpos_output,'*linear');

% for (n=1:24)
% psi_value_corr_output(n,:)=psipos_output(n,:)-alphas_psi_star_ini(n)+psi_star_output(n,:);
% end
lambda_P1=lambda_output(P1,1);
lambda_P2=lambda_output(P2,1);
lambda_P5=lambda_output(P5,1);
lambda_P6=lambda_output(P6,1);

Ekin_P1=round(Ekin_output(P1,1));
Ekin_P2=round(Ekin_output(P2,1));
Ekin_P5=round(Ekin_output(P5,1));
Ekin_P6=round(Ekin_output(P6,1));

lP1=strcat('\lambda _{ini} = ',num2str(lambda_P1));
lP2=strcat('\lambda _{ini} = ',num2str(lambda_P2));
lP5=strcat('\lambda _{ini} = ',num2str(lambda_P5));
lP6=strcat('\lambda _{ini} = ',num2str(lambda_P6));

EP1=strcat('Ekin _{ini} = ',num2str(Ekin_P1));
EP2=strcat('Ekin _{ini} = ',num2str(Ekin_P2));
EP5=strcat('Ekin _{ini} = ',num2str(Ekin_P5));
EP6=strcat('Ekin _{ini} = ',num2str(Ekin_P6));


phipos_output_wrap=wrap2pi(phipos_output);
omega_output=theta_output-phipos_output_wrap;
omega_output=wrap2pi(omega_output);

size_phi=257;
scale_tor=40*pi*((0:size_phi-1)/(size_phi-1));
Btot_X_vector=mean(Btot_XZ_map,2);
[Btot_Xphi_map tmp]=meshgrid(Btot_X_vector,ones(size_phi,1));
Btot_Xphi_map=Btot_Xphi_map';

end_time=size(time_scale,2);
ejection_time=size(time_scale,2);
% end_time=250;
%end_time=150000;
time_end=time_scale(end_time);

%time0=round(0.75*end_time);
time0=1;


close all

figure(1);
set(gca,'FontSize',16);
hold on;

plot(Xpos_output(P1,time0:end_time),Zpos_output(P1,time0:end_time),'b','LineWidth',2);
plot(Xpos_output(P2,time0:end_time),Zpos_output(P2,time0:end_time),'Color',[0.1 0.1 0.6],'LineWidth',2);
plot(Xpos_output(P5,time0:end_time),Zpos_output(P5,time0:end_time),'Color',[1.0 0.2 0.2],'LineWidth',1);
plot(Xpos_output(P6,time0:end_time),Zpos_output(P6,time0:end_time),'Color',[0.6 0.2 0.2],'LineWidth',1);
legend(lP1,lP2,lP5,lP6);

contour(scale_X,scale_Z,(radial_XZsmall_map'-257)*psi_global/257,12,'y');
xlim([-0.5 0.8]);
ylim([-0.6 0.7]);

contour(scale_X,scale_Z,Btot_XZ_map',(2:0.2:4.4));axis xy square
colormap('summer')
colorbar;

plot(Xpos_output(P1,time0:end_time),Zpos_output(P1,time0:end_time),'b','LineWidth',2);
plot(Xpos_output(P2,time0:end_time),Zpos_output(P2,time0:end_time),'Color',[0.1 0.1 0.6],'LineWidth',2);
plot(Xpos_output(P5,time0:end_time),Zpos_output(P5,time0:end_time),'Color',[1.0 0.2 0.2],'LineWidth',1);
plot(Xpos_output(P6,time0:end_time),Zpos_output(P6,time0:end_time),'Color',[0.6 0.2 0.2],'LineWidth',1);

xlabel('X (m)');
ylabel('Z (m)');



% figure(2);
% set(gca,'FontSize',16);
% 
% imagesc(scale_X,scale_tor,Btot_Xphi_map',[1 4]);
% axis xy square
% xlim([-0.5 0.5]);
% ylim([0 40*pi]);
% colormap('summer')
% %brighten(0.6);
% hold on;
% plot(Xpos_output(P4,time0:end_time),phipos_output(P4,time0:end_time),'b','LineWidth',2);
% plot(Xpos_output(P5,time0:end_time),phipos_output(P5,time0:end_time),'Color',[0.2 0.8 0.2],'LineWidth',1);
% plot(Xpos_output(P6,time0:end_time),phipos_output(P6,time0:end_time),'Color',[0.8 0.2 0.2],'LineWidth',1);
% 
% xlabel('X (m)');
% ylabel('\phi (rad)');
% colorbar;
% 




figure(3);

subplot(2,1,1);
set(gca,'FontSize',16);
grid on;
hold on;
xlim([time_scale(1) time_scale(end-1)])
plot(time_scale(time0:end_time),phipos_output(P1,time0:end_time),'b','LineWidth',2);
plot(time_scale(time0:end_time),phipos_output(P2,time0:end_time),'Color',[0.1 0.1 0.6],'LineWidth',2);
plot(time_scale(time0:end_time),phipos_output(P5,time0:end_time),'Color',[0.2 0.8 0.2],'LineWidth',2);
plot(time_scale(time0:end_time),phipos_output(P6,time0:end_time),'Color',[0.8 0.2 0.2],'LineWidth',2);
xlabel('time (s)');
ylabel('\phi (rad)');
xlim([0 time_end]);
legend(lP1,lP2,lP5,lP6);

subplot(2,1,2);
set(gca,'FontSize',16);
grid on;
hold on
xlim([time_scale(1) time_scale(end-1)])
plot(time_scale(time0:end_time),vparallel_output(P1,time0:end_time),'b','LineWidth',2);
plot(time_scale(time0:end_time),vparallel_output(P2,time0:end_time),'Color',[0.1 0.1 0.6],'LineWidth',2);
plot(time_scale(time0:end_time),vparallel_output(P5,time0:end_time),'Color',[0.2 0.8 0.2],'LineWidth',2);
plot(time_scale(time0:end_time),vparallel_output(P6,time0:end_time),'Color',[0.8 0.2 0.2],'LineWidth',2);
xlabel('time (s)');
ylabel('v_{||} (m/s)');
xlim([0 time_end]);



figure(4)
set(gca,'FontSize',16);

grid on;
hold on;
xlim([time_scale(1) time_scale(end-1)])
% plot(time_scale(1:end_time),0.5*(mHe/eV)*(vparallel_output(P1,1:end_time).^2+vdriftsq_output(P1,1:end_time))+Eperp_output(P1,1:end_time),'b','LineWidth',2);
% plot(time_scale(1:end_time),0.5*(mHe/eV)*(vparallel_output(P2,1:end_time).^2+vdriftsq_output(P2,1:end_time))+Eperp_output(P2,1:end_time),'Color',[0.1 0.1 0.6],'LineWidth',2);
% plot(time_scale(1:end_time),0.5*(mHe/eV)*(vparallel_output(P5,1:end_time).^2+vdriftsq_output(P5,1:end_time))+Eperp_output(P5,1:end_time),'Color',[0.2 1.0 0.2],'LineWidth',2);
% plot(time_scale(1:end_time),0.5*(mHe/eV)*(vparallel_output(P6,1:end_time).^2+vdriftsq_output(P6,1:end_time))+Eperp_output(P6,1:end_time),'Color',[0.8 0.2 0.2],'LineWidth',2);
plot(time_scale(1:end_time),0.5*(mHe/eV)*(vparallel_output(P1,1:end_time).^2)+Eperp_output(P1,1:end_time),'b','LineWidth',2);
plot(time_scale(1:end_time),0.5*(mHe/eV)*(vparallel_output(P2,1:end_time).^2)+Eperp_output(P2,1:end_time),'Color',[0.1 0.1 0.6],'LineWidth',2);
plot(time_scale(1:end_time),0.5*(mHe/eV)*(vparallel_output(P5,1:end_time).^2)+Eperp_output(P5,1:end_time),'Color',[0.2 1.0 0.2],'LineWidth',2);
plot(time_scale(1:end_time),0.5*(mHe/eV)*(vparallel_output(P6,1:end_time).^2)+Eperp_output(P6,1:end_time),'Color',[0.8 0.2 0.2],'LineWidth',2);
xlabel('time (s)');
ylabel('Ekin (eV)');
xlim([0 time_end]);

xlim([time_scale(1) time_scale(end-1)])

legend(lP1,lP2,lP5,lP6);








figure(5);
set(gca,'FontSize',16);

grid on;
hold on;
xlim([time_scale(1) time_scale(end-1)])
plot(time_scale(1:end_time),(mHe/eV)*(Xpos_output(P1,1:end_time)+R0).*(vphi_output(P1,1:end_time))-ZHe*(psi_value_output(P1,1:end_time)),'b','LineWidth',2);
plot(time_scale(1:end_time),(mHe/eV)*(Xpos_output(P2,1:end_time)+R0).*(vphi_output(P2,1:end_time))-ZHe*(psi_value_output(P2,1:end_time)),'Color',[0.1 0.1 0.6],'LineWidth',2);
plot(time_scale(1:end_time),(mHe/eV)*(Xpos_output(P5,1:end_time)+R0).*(vphi_output(P5,1:end_time))-ZHe*(psi_value_output(P5,1:end_time)),'Color',[0.2 0.8 0.2],'LineWidth',2);
plot(time_scale(1:end_time),(mHe/eV)*(Xpos_output(P6,1:end_time)+R0).*(vphi_output(P6,1:end_time))-ZHe*(psi_value_output(P6,1:end_time)),'Color',[0.8 0.2 0.2],'LineWidth',2);
xlabel('time (s)');
ylabel('p_\phi');
xlim([0 time_end]);





figure(6)

subplot(2,1,1);
set(gca,'FontSize',16);

grid on;
hold on
xlim([time_scale(1) time_scale(end-1)])
% plot(time_scale(1:end_time),Eperp_output(P1,1:end_time),'r');
plot(time_scale(1:end_time),Eperp_output(P1,1:end_time),'b','LineWidth',2);
plot(time_scale(1:end_time),Eperp_output(P2,1:end_time),'Color',[0.1 0.1 0.6],'LineWidth',2);
plot(time_scale(1:end_time),Eperp_output(P5,1:end_time),'Color',[0.2 0.8 0.2],'LineWidth',2);
plot(time_scale(1:end_time),Eperp_output(P6,1:end_time),'Color',[0.8 0.2 0.2],'LineWidth',2);
xlabel('time (s)');
ylabel('Eperp (eV)');
legend(lP1,lP2,lP5,lP6);

subplot(2,1,2);
set(gca,'FontSize',16);

grid on;
hold on
xlim([time_scale(1) time_scale(end-1)])
% plot(time_scale(1:end_time),vparallel_output(P1,1:end_time),'r');
plot(time_scale(1:end_time),vparallel_output(P1,1:end_time),'b','LineWidth',2);
plot(time_scale(1:end_time),vparallel_output(P2,1:end_time),'Color',[0.1 0.1 0.6],'LineWidth',2);
plot(time_scale(1:end_time),vparallel_output(P5,1:end_time),'Color',[0.2 0.8 0.2],'LineWidth',2);
plot(time_scale(1:end_time),vparallel_output(P6,1:end_time),'Color',[0.8 0.2 0.2],'LineWidth',2);
xlabel('time (s)');
ylabel('v_{||} (m/s)');





figure(7)

subplot(2,1,1);
set(gca,'FontSize',16);
grid on;
hold on

plot(time_scale(time0:end_time),mm_output(P1,time0:end_time),'b','LineWidth',1.5);
plot(time_scale(time0:end_time),mm_output(P2,time0:end_time),'Color',[0.1 0.1 0.6],'LineWidth',1.5);
plot(time_scale(time0:end_time),mm_output(P5,time0:end_time),'Color',[0.2 0.8 0.2],'LineWidth',1.5);
plot(time_scale(time0:end_time),mm_output(P6,time0:end_time),'Color',[0.8 0.2 0.2],'LineWidth',1.5);
xlabel('time (s)');
ylabel('\mu (eV T^{-1})');
% xlim([time_scale(1) time_scale(end-1)])
xlim([0 time_end]);
legend(lP1,lP2,lP5,lP6);

subplot(2,1,2);
set(gca,'FontSize',16);
grid on;
hold on
plot(time_scale(time0:end_time),lambda_output(P1,time0:end_time),'b','LineWidth',1.5);
plot(time_scale(time0:end_time),lambda_output(P2,time0:end_time),'Color',[0.1 0.1 0.6],'LineWidth',1.5);
plot(time_scale(time0:end_time),lambda_output(P5,time0:end_time),'Color',[0.2 0.8 0.2],'LineWidth',1.5);
plot(time_scale(time0:end_time),lambda_output(P6,time0:end_time),'Color',[0.8 0.2 0.2],'LineWidth',1.5);
xlabel('time (s)')
ylabel('\lambda (Bavg)');
% xlim([time_scale(1) time_scale(end-1)])
xlim([0 time_end]);



figure(8)

subplot(2,1,1);
set(gca,'FontSize',16);
grid on;
hold on

plot(time_scale(time0:end_time),omega_output(P1,time0:end_time),'b','LineWidth',1.5);
plot(time_scale(time0:end_time),omega_output(P2,time0:end_time),'Color',[0.1 0.1 0.6],'LineWidth',1.5);
plot(time_scale(time0:end_time),omega_output(P5,time0:end_time),'Color',[0.2 0.8 0.2],'LineWidth',1.5);
plot(time_scale(time0:end_time),omega_output(P6,time0:end_time),'Color',[0.8 0.2 0.2],'LineWidth',1.5);
xlabel('time (s)');
ylabel('\omega (rad)');
% xlim([time_scale(1) time_scale(end-1)])
xlim([0 time_end]);

subplot(2,1,2);
set(gca,'FontSize',16);
grid on;
hold on
plot(time_scale(time0:end_time),(vEradial_output(P1,time0:end_time)),'b','LineWidth',1.5);
plot(time_scale(time0:end_time),(vEradial_output(P2,time0:end_time)),'Color',[0.1 0.1 0.6],'LineWidth',1.5);
plot(time_scale(time0:end_time),(vEradial_output(P5,time0:end_time)),'Color',[0.2 0.8 0.2],'LineWidth',1.5);
plot(time_scale(time0:end_time),(vEradial_output(P6,time0:end_time)),'Color',[0.8 0.2 0.2],'LineWidth',1.5);
xlabel('time (s)')
ylabel('v_E (radial) (m/s)');
% xlim([time_scale(1) time_scale(end-1)])
xlim([0 time_end]);

legend(EP1,EP2,EP5,EP6);


figure(9)

subplot(2,1,1);
set(gca,'FontSize',16);
grid on;
hold on

plot(time_scale(time0:end_time),psi_value_output(P1,time0:end_time),'b','LineWidth',1.5);
plot(time_scale(time0:end_time),psi_value_output(P2,time0:end_time),'Color',[0.1 0.1 0.6],'LineWidth',1.5);
plot(time_scale(time0:end_time),psi_value_output(P5,time0:end_time),'Color',[0.2 0.8 0.2],'LineWidth',1.5);
plot(time_scale(time0:end_time),psi_value_output(P6,time0:end_time),'Color',[0.8 0.2 0.2],'LineWidth',1.5);
xlabel('time (s)');
ylabel('\psi (T.m^{-2})');
% xlim([time_scale(1) time_scale(end-1)])
xlim([0 time_end]);



subplot(2,1,2);
set(gca,'FontSize',16);
grid on;
hold on
plot(time_scale(time0:end_time),sqrt(vEsq_output(P1,time0:end_time)),'b','LineWidth',1.5);
plot(time_scale(time0:end_time),sqrt(vEsq_output(P2,time0:end_time)),'Color',[0.1 0.1 0.6],'LineWidth',1.5);
plot(time_scale(time0:end_time),sqrt(vEsq_output(P5,time0:end_time)),'Color',[0.2 0.8 0.2],'LineWidth',1.5);
plot(time_scale(time0:end_time),sqrt(vEsq_output(P6,time0:end_time)),'Color',[0.8 0.2 0.2],'LineWidth',1.5);
xlabel('time (s)')
ylabel('v_E (m/s)');
% xlim([time_scale(1) time_scale(end-1)])
xlim([0 time_end]);

legend(lP1,lP2,lP5,lP6);



figure(10)

subplot(2,1,1);
set(gca,'FontSize',16);
grid on;
hold on
xlim([time_scale(1) time_scale(end-1)])
% plot(time_scale(time0:end_time),psi_star_output(P1,time0:end_time),'b');
% plot(time_scale(time0:end_time),psi_star_output(P2,time0:end_time),'Color',[0.1 0.1 0.6],'LineWidth',1.5);
% plot(time_scale(time0:end_time),psi_star_output(P5,time0:end_time),'Color',[0.2 0.8 0.2],'LineWidth',1.5);
% plot(time_scale(time0:end_time),psi_star_output(P6,time0:end_time),'Color',[0.8 0.2 0.2],'LineWidth',1.5);
% xlabel('time (s)');
% ylabel('\psi*');
plot(time_scale(time0:end_time),Epot_output(P1,time0:end_time),'b');
plot(time_scale(time0:end_time),Epot_output(P2,time0:end_time),'Color',[0.1 0.1 0.6],'LineWidth',1.5);
plot(time_scale(time0:end_time),Epot_output(P5,time0:end_time),'Color',[0.2 0.8 0.2],'LineWidth',1.5);
plot(time_scale(time0:end_time),Epot_output(P6,time0:end_time),'Color',[0.8 0.2 0.2],'LineWidth',1.5);
xlabel('time (s)');
ylabel('Epot (eV)');

subplot(2,1,2);
set(gca,'FontSize',16);
grid on;
hold on
plot(time_scale(1:end_time),Epot_output(P1,1:end_time)+0.5*(mHe/eV)*(vparallel_output(P1,1:end_time).^2)+Eperp_output(P1,1:end_time),'b','LineWidth',2);
plot(time_scale(1:end_time),Epot_output(P2,1:end_time)+0.5*(mHe/eV)*(vparallel_output(P2,1:end_time).^2)+Eperp_output(P2,1:end_time),'Color',[0.1 0.1 0.6],'LineWidth',2);
plot(time_scale(1:end_time),Epot_output(P5,1:end_time)+0.5*(mHe/eV)*(vparallel_output(P5,1:end_time).^2)+Eperp_output(P5,1:end_time),'Color',[0.2 1.0 0.2],'LineWidth',2);
plot(time_scale(1:end_time),Epot_output(P6,1:end_time)+0.5*(mHe/eV)*(vparallel_output(P6,1:end_time).^2)+Eperp_output(P6,1:end_time),'Color',[0.8 0.2 0.2],'LineWidth',2);
xlabel('time (s)');
ylabel('Etot (eV)');



figure(11)

subplot(2,1,1);
set(gca,'FontSize',16);
grid on;
hold on
xlim([time_scale(1) time_scale(end-1)])
plot(time_scale(time0:end_time),psi_star_output(P1,time0:end_time),'b');
plot(time_scale(time0:end_time),psi_star_output(P2,time0:end_time),'Color',[0.1 0.1 0.6],'LineWidth',1.5);
plot(time_scale(time0:end_time),psi_star_output(P5,time0:end_time),'Color',[0.2 0.8 0.2],'LineWidth',1.5);
plot(time_scale(time0:end_time),psi_star_output(P6,time0:end_time),'Color',[0.8 0.2 0.2],'LineWidth',1.5);
xlabel('time (s)');
ylabel('\psi*');

subplot(2,1,2);
set(gca,'FontSize',16);
grid on;
hold on
plot(time_scale(time0:end_time),(vphi_output(P1,time0:end_time)),'b','LineWidth',1.5);
plot(time_scale(time0:end_time),(vphi_output(P2,time0:end_time)),'Color',[0.1 0.1 0.6],'LineWidth',1.5);
plot(time_scale(time0:end_time),(vphi_output(P5,time0:end_time)),'Color',[0.2 0.8 0.2],'LineWidth',1.5);
plot(time_scale(time0:end_time),(vphi_output(P6,time0:end_time)),'Color',[0.8 0.2 0.2],'LineWidth',1.5);
xlabel('time (s)')
ylabel('v_\phi (m/s)');
% xlim([time_scale(1) time_scale(end-1)])
xlim([0 time_end]);

legend(EP1,EP2,EP5,EP6);


% 
% figure(12)
% 
% subplot(2,1,1);
% set(gca,'FontSize',16);
% grid on;
% hold on
% xlim([time_scale(1) time_scale(end-1)])
% % plot(time_scale(time0:end_time),psi_star_output(P1,time0:end_time),'b');
% % plot(time_scale(time0:end_time),psi_star_output(P2,time0:end_time),'Color',[0.1 0.1 0.6],'LineWidth',1.5);
% % plot(time_scale(time0:end_time),psi_star_output(P5,time0:end_time),'Color',[0.2 0.8 0.2],'LineWidth',1.5);
% % plot(time_scale(time0:end_time),psi_star_output(P6,time0:end_time),'Color',[0.8 0.2 0.2],'LineWidth',1.5);
% % xlabel('time (s)');
% % ylabel('\psi*');
% plot(time_scale(time0:end_time),theta_output(P1,time0:end_time),'b');
% plot(time_scale(time0:end_time),theta_output(P2,time0:end_time),'Color',[0.1 0.1 0.6],'LineWidth',1.5);
% plot(time_scale(time0:end_time),theta_output(P5,time0:end_time),'Color',[0.2 0.8 0.2],'LineWidth',1.5);
% plot(time_scale(time0:end_time),theta_output(P6,time0:end_time),'Color',[0.8 0.2 0.2],'LineWidth',1.5);
% xlabel('time (s)');
% ylabel('Epot (eV)');
% 
% subplot(2,1,2);
% set(gca,'FontSize',16);
% grid on;
% hold on
% plot(time_scale(1:end_time),vparallel_output(P1,1:end_time),'b','LineWidth',2);
% plot(time_scale(1:end_time),vparallel_output(P2,1:end_time),'Color',[0.1 0.1 0.6],'LineWidth',2);
% plot(time_scale(1:end_time),vparallel_output(P5,1:end_time),'Color',[0.2 0.8 0.2],'LineWidth',2);
% plot(time_scale(1:end_time),vparallel_output(P6,1:end_time),'Color',[0.8 0.2 0.2],'LineWidth',2);
% 
% xlabel('time (s)');
% ylabel('Etot (eV)');
% 
% legend(EP1,EP2,EP5,EP6);
