
% % 
P4=1;
P5=4;
P6=7;
% 


% vperp_sq_output=vDsq_output;
% Bfield_part_output=interp2(scale_X,scale_Z,Btot_XZ_map',Xpos_output,Zpos_output,'*linear');
mm_output=Eperp_output./Bfield_output;
% psi_value_output=interp2(scale_X,scale_Z,psi_XZsmall_map',Xpos_output,Zpos_output,'*linear');
% lap_psi_output=interp2(scale_X,scale_Z,lap_psi_XZsmall_map',Xpos_output,Zpos_output,'*linear');

% for (n=1:24)
% psi_value_corr_output(n,:)=psipos_output(n,:)-alphas_psi_star_ini(n)+psi_star_output(n,:);
% end


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

contour(scale_X,scale_Z,(radial_XZsmall_map'-257)*psi_global/257,12,'y');
hold on;

xlim([-0.5 1.0]);
ylim([-0.7 0.8]);

contour(scale_X,scale_Z,Btot_XZ_map',(2:0.2:4.4));axis xy square
colormap('summer')
colorbar;

plot(Xpos_output(P4,time0:end_time),Zpos_output(P4,time0:end_time),'b','LineWidth',2);
plot(Xpos_output(P5,time0:end_time),Zpos_output(P5,time0:end_time),'g--','LineWidth',1);
plot(Xpos_output(P6,time0:end_time),Zpos_output(P6,time0:end_time),'Color',[0.4 0.4 0.4],'LineWidth',1);

xlabel('X (m)');
ylabel('Z (m)');



figure(2);
set(gca,'FontSize',16);

imagesc(scale_X,scale_tor,Btot_Xphi_map',[1 4]);
axis xy square
xlim([-0.5 0.5]);
ylim([0 40*pi]);
colormap('summer')
%brighten(0.6);
hold on;
plot(Xpos_output(P4,time0:end_time),phipos_output(P4,time0:end_time),'b','LineWidth',2);
plot(Xpos_output(P5,time0:end_time),phipos_output(P5,time0:end_time),'g--','LineWidth',1);
plot(Xpos_output(P6,time0:end_time),phipos_output(P6,time0:end_time),'Color',[0.4 0.4 0.4],'LineWidth',1);

xlabel('X (m)');
ylabel('\phi (rad)');
colorbar;





figure(3);

subplot(2,1,1);
set(gca,'FontSize',16);
grid on;
hold on;
xlim([time_scale(1) time_scale(end-1)])
plot(time_scale(time0:end_time),phipos_output(P4,time0:end_time),'b','LineWidth',2);
plot(time_scale(time0:end_time),phipos_output(P5,time0:end_time),'Color',[0.2 0.8 0.2],'LineWidth',1);
plot(time_scale(time0:end_time),phipos_output(P6,time0:end_time),'Color',[0.4 0.4 0.4],'LineWidth',1);
xlabel('time (s)');
ylabel('\phi (rad)');
xlim([0 time_end]);

subplot(2,1,2);
set(gca,'FontSize',16);
grid on;
hold on
xlim([time_scale(1) time_scale(end-1)])
plot(time_scale(time0:end_time),vparallel_output(P4,time0:end_time),'b','LineWidth',2);
plot(time_scale(time0:end_time),vparallel_output(P5,time0:end_time),'Color',[0.2 0.8 0.2],'LineWidth',1);
plot(time_scale(time0:end_time),vparallel_output(P6,time0:end_time),'Color',[0.4 0.4 0.4],'LineWidth',1);
xlabel('time (s)');
ylabel('v_{||} (m/s)');
xlim([0 time_end]);



figure(4)
set(gca,'FontSize',16);

grid on;
hold on;
xlim([time_scale(1) time_scale(end-1)])
plot(time_scale(1:end_time),0.5*(mHe/eV)*(vparallel_output(P4,1:end_time).^2+vdriftsq_output(P4,1:end_time))+Eperp_output(P4,1:end_time),'b','LineWidth',2);
plot(time_scale(1:end_time),0.5*(mHe/eV)*(vparallel_output(P5,1:end_time).^2+vdriftsq_output(P5,1:end_time))+Eperp_output(P5,1:end_time),'Color',[0.2 0.8 0.2],'LineWidth',1);
plot(time_scale(1:end_time),0.5*(mHe/eV)*(vparallel_output(P6,1:end_time).^2+vdriftsq_output(P6,1:end_time))+Eperp_output(P6,1:end_time),'Color',[0.4 0.4 0.4],'LineWidth',1);
xlabel('time (s)');
ylabel('Ekin (eV)');
xlim([0 time_end]);

xlim([time_scale(1) time_scale(end-1)])



% figure(5)
% set(gca,'FontSize',22);
% 
% grid on;
% hold on
% xlim([time_scale(1) time_scale(end-1)])
% plot(time_scale(time0:end_time),psipos_output(P1,time0:end_time),'r','LineWidth',2);
% % plot(time_scale(time0:end_time),psipos_output(P2,time0:end_time),'--','Color',[0.2 0.8 0.2]);
% % plot(time_scale(time0:end_time),psipos_output(P3,time0:end_time),'b--');
% plot(time_scale(time0:end_time),psipos_output(P4,time0:end_time),'b','LineWidth',2);
% plot(time_scale(time0:end_time),psipos_output(P5,time0:end_time),'Color',[0.2 0.8 0.2],'LineWidth',2);
% plot(time_scale(time0:end_time),psipos_output(P6,time0:end_time),'Color',[0.4 0.4 0.4],'LineWidth',2);
% xlabel('time (s)');
% ylabel('\psi');
% 
% xlim([0 time_end]);


% figure(6);
% 
% subplot(2,1,1);
% set(gca,'FontSize',16);
% grid on;
% hold on;
% xlim([time_scale(1) time_scale(end-1)])
% % plot(time_scale(1:end_time),sqrt(Omega_output(P1,1:end_time)),'r','LineWidth',2);
% % plot(time_scale(1:end_time),sqrt(vDsq_output(P2,1:end_time)),'--','Color',[0.2 0.8 0.2]);
% % plot(time_scale(1:end_time),sqrt(vDsq_output(P3,1:end_time)),'b--');
% plot(time_scale(1:end_time),sqrt(Omega_output(P4,1:end_time)),'b','LineWidth',2);
% plot(time_scale(1:end_time),sqrt(Omega_output(P5,1:end_time)),'Color',[0.2 0.8 0.2],'LineWidth',2);
% plot(time_scale(1:end_time),sqrt(Omega_output(P6,1:end_time)),'Color',[0.4 0.4 0.4],'LineWidth',2);
% xlabel('time (s)');
% ylabel('\Omega (rad/s)');
% xlim([0 time_end]);
% 
% subplot(2,1,2);
% set(gca,'FontSize',16);
% grid on;
% hold on
% xlim([time_scale(1) time_scale(end-1)])
% % plot(time_scale(1:end_time),omega_output(P1,1:end_time),'r','LineWidth',2);
% % plot(time_scale(1:end_time),omega_output(P2,1:end_time),'--','Color',[0.2 0.8 0.2]);
% % plot(time_scale(1:end_time),omega_output(P3,1:end_time),'b--');
% plot(time_scale(1:end_time),omega_output(P4,1:end_time),'b','LineWidth',2);
% plot(time_scale(1:end_time),omega_output(P5,1:end_time),'Color',[0.2 0.8 0.2],'LineWidth',2);
% plot(time_scale(1:end_time),omega_output(P6,1:end_time),'Color',[0.4 0.4 0.4],'LineWidth',2);
% xlabel('time (s)');
% ylabel('\omega (rad)');
% xlim([0 time_end]);



figure(7);
set(gca,'FontSize',16);

grid on;
hold on;
xlim([time_scale(1) time_scale(end-1)])
plot(time_scale(1:end_time),(mHe/eV)*(Xpos_output(P4,1:end_time)+R0).*(vphi_output(P4,1:end_time))-ZHe*(psi_value_output(P4,1:end_time)),'b','LineWidth',2);
plot(time_scale(1:end_time),(mHe/eV)*(Xpos_output(P5,1:end_time)+R0).*(vphi_output(P5,1:end_time))-ZHe*(psi_value_output(P5,1:end_time)),'Color',[0.2 0.8 0.2],'LineWidth',1);
plot(time_scale(1:end_time),(mHe/eV)*(Xpos_output(P6,1:end_time)+R0).*(vphi_output(P6,1:end_time))-ZHe*(psi_value_output(P6,1:end_time)),'Color',[0.4 0.4 0.4],'LineWidth',1);
xlabel('time (s)');
ylabel('p_\phi');
xlim([0 time_end]);



% figure(8)
% grid on;
% hold on
% xlim([time_scale(1) time_scale(end-1)])
% plot(time_scale(time0:end_time),psi_star_output(P1,time0:end_time),'k--');
% plot(time_scale(time0:end_time),psi_star_output(P2,time0:end_time),'--','Color',[0.2 0.8 0.2]);
% plot(time_scale(time0:end_time),psi_star_output(P3,time0:end_time),'b--');
% plot(time_scale(time0:end_time),psi_star_output(P4,time0:end_time),'b','LineWidth',2);
% plot(time_scale(time0:end_time),psi_star_output(P5,time0:end_time),'Color',[0.2 0.8 0.2],'LineWidth',2);
% plot(time_scale(time0:end_time),psi_star_output(P6,time0:end_time),'Color',[0.4 0.4 0.4],'LineWidth',2);
% xlabel('time (s)');
% ylabel('\psi');


figure(9)

subplot(2,1,1);
set(gca,'FontSize',16);

grid on;
hold on
xlim([time_scale(1) time_scale(end-1)])
% plot(time_scale(1:end_time),Eperp_output(P1,1:end_time),'r');
plot(time_scale(1:end_time),Eperp_output(P4,1:end_time),'b','LineWidth',2);
plot(time_scale(1:end_time),Eperp_output(P5,1:end_time),'Color',[0.2 0.8 0.2],'LineWidth',1);
plot(time_scale(1:end_time),Eperp_output(P6,1:end_time),'Color',[0.4 0.4 0.4],'LineWidth',1);
xlabel('time (s)');
ylabel('Eperp (eV)');

subplot(2,1,2);
set(gca,'FontSize',16);

grid on;
hold on
xlim([time_scale(1) time_scale(end-1)])
% plot(time_scale(1:end_time),vparallel_output(P1,1:end_time),'r');
plot(time_scale(1:end_time),vparallel_output(P4,1:end_time),'b','LineWidth',2);
plot(time_scale(1:end_time),vparallel_output(P5,1:end_time),'Color',[0.2 0.8 0.2],'LineWidth',1);
plot(time_scale(1:end_time),vparallel_output(P6,1:end_time),'Color',[0.4 0.4 0.4],'LineWidth',1);
xlabel('time (s)');
ylabel('v_{||} (m/s)');


% figure(10)
% set(gca,'FontSize',16);
% grid on;
% hold on
% xlim([time_scale(1) time_scale(end-1)])
% plot(time_scale(time0:end_time),pphi_output(P1,time0:end_time),'r');
% plot(time_scale(time0:end_time),pphi_output(P4,time0:end_time),'b','LineWidth',2);
% plot(time_scale(time0:end_time),pphi_output(P5,time0:end_time),'Color',[0.2 0.8 0.2],'LineWidth',2);
% plot(time_scale(time0:end_time),pphi_output(P6,time0:end_time),'Color',[0.4 0.4 0.4],'LineWidth',2);
% xlabel('time (s)');
% ylabel('pphi0');


figure(11)

subplot(2,1,1);
set(gca,'FontSize',16);
grid on;
hold on
xlim([time_scale(1) time_scale(end-1)])
plot(time_scale(time0:end_time),psi_star_output(P4,time0:end_time),'b','LineWidth',2);
plot(time_scale(time0:end_time),psi_star_output(P5,time0:end_time),'Color',[0.2 0.8 0.2],'LineWidth',2);
plot(time_scale(time0:end_time),psi_star_output(P6,time0:end_time),'Color',[0.4 0.4 0.4],'LineWidth',2);
xlabel('time (s)');
ylabel('\psi _*');

subplot(2,1,2);
set(gca,'FontSize',16);
grid on;
hold on
plot(time_scale(time0:end_time),sqrt(vDsq_output(P4,time0:end_time)),'b','LineWidth',1.5);
plot(time_scale(time0:end_time),sqrt(vDsq_output(P5,time0:end_time)),'Color',[0.2 0.8 0.2],'LineWidth',1.5);
plot(time_scale(time0:end_time),sqrt(vDsq_output(P6,time0:end_time)),'Color',[0.4 0.4 0.4],'LineWidth',1.5);
xlabel('time (s)')
ylabel('v_D (m/s)');


% curl_corr_output=interp2(scale_X,scale_Z,rotational_b_pll',Xpos_output,Zpos_output,'*linear');
% 
% mm_recalc_output=(1.55e4)*(1+(mHe/(eV*ZHe))*vparallel_output.*curl_corr_output);

% figure(12)
% set(gca,'FontSize',16);
% grid on;
% hold on
% xlim([time_scale(1) time_scale(end-1)])
% plot(time_scale(time0:end_time),mm_output(P4,time0:end_time),'b','LineWidth',2);
% plot(time_scale(time0:end_time),mm_output(P5,time0:end_time),'Color',[0.2 0.8 0.2],'LineWidth',1);
% plot(time_scale(time0:end_time),mm_output(P6,time0:end_time),'Color',[0.4 0.4 0.4],'LineWidth',1);
% xlabel('time (s)');
% ylabel('\mu (eV T^{-1})');
% 



figure(7)

subplot(2,1,1);
set(gca,'FontSize',16);
grid on;
hold on

plot(time_scale(time0:end_time),mm_output(P4,time0:end_time),'b','LineWidth',1.5);
plot(time_scale(time0:end_time),mm_output(P5,time0:end_time),'Color',[0.2 0.8 0.2],'LineWidth',1.5);
plot(time_scale(time0:end_time),mm_output(P6,time0:end_time),'Color',[0.4 0.4 0.4],'LineWidth',1.5);
xlabel('time (s)');
ylabel('\mu (eV T^{-1})');
% xlim([time_scale(1) time_scale(end-1)])
xlim([0 time_end]);

subplot(2,1,2);
set(gca,'FontSize',16);
grid on;
hold on
plot(time_scale(time0:end_time),sqrt(vEsq_output(P4,time0:end_time)),'b','LineWidth',1.5);
plot(time_scale(time0:end_time),sqrt(vEsq_output(P5,time0:end_time)),'Color',[0.2 0.8 0.2],'LineWidth',1.5);
plot(time_scale(time0:end_time),sqrt(vEsq_output(P6,time0:end_time)),'Color',[0.4 0.4 0.4],'LineWidth',1.5);
xlabel('time (s)')
ylabel('v_E (m/s)');
% xlim([time_scale(1) time_scale(end-1)])
xlim([0 time_end]);
