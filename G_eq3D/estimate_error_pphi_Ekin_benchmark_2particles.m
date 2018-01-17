% time_stamp_end = 1231
% find(Ekin_outputG(:,1)<60*1e3)
TIME_INI=1000

% % 
P1=60;
P2=2+32;
P5=1;
P6=2;


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
lambda_P5=lambda_outputG(P5,1);
lambda_P6=lambda_outputG(P6,1);

Ekin_P1=round(Ekin_outputG(P1,1));
Ekin_P2=round(Ekin_outputG(P2,1));
Ekin_P5=round(Ekin_outputG(P5,1));
Ekin_P6=round(Ekin_outputG(P6,1));

lP1=strcat('\lambda _{0} = ',num2str(lambda_P1));
lP2=strcat('\lambda _{0} = ',num2str(lambda_P2));
lP5=strcat('\lambda _{0} = ',num2str(lambda_P5));
lP6=strcat('\lambda _{0} = ',num2str(lambda_P6));

EP1=strcat('Ekin _{ini} = ',num2str(Ekin_P1));
EP2=strcat('Ekin _{ini} = ',num2str(Ekin_P2));
EP5=strcat('Ekin _{ini} = ',num2str(Ekin_P5));
EP6=strcat('Ekin _{ini} = ',num2str(Ekin_P6));



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


figure(1);

subplot(2,1,1);
set(gca,'FontSize',16);
grid on;
hold on;
xlim([time_scale_G(1) time_scale_G(end-1)])
plot(time_scale_G(time0:end_time),phipos_outputG(P1,time0:end_time),'b','LineWidth',2);
plot(time_scale_G(time0:end_time),phipos_outputG(P5,time0:end_time),'Color',[0.2 0.8 0.2],'LineWidth',2);
xlabel('time (s)');
ylabel('\phi (rad)');
xlim([0 time_end]);
legend(lP1,lP5);

subplot(2,1,2);
set(gca,'FontSize',16);
grid on;
hold on
xlim([time_scale_G(1) time_scale_G(end-1)])
plot(time_scale_G(time0:end_time),Xpos_outputG(P1,time0:end_time),'b','LineWidth',2);
plot(time_scale_G(time0:end_time),Xpos_outputG(P5,time0:end_time),'Color',[0.2 0.8 0.2],'LineWidth',2);
xlabel('time (s)');
ylabel('X (m)');
xlim([0 time_end]);







figure(3);

subplot(2,1,1);
set(gca,'FontSize',16);
grid on;
hold on;
xlim([time_scale_G(1) time_scale_G(end-1)])
plot(time_scale_G(time0:end_time),mm_outputG(P1,time0:end_time),'b','LineWidth',2);
plot(time_scale_G(time0:end_time),mm_outputG(P5,time0:end_time),'Color',[0.2 0.8 0.2],'LineWidth',2);
xlabel('time (s)');
ylabel('\mu)');
xlim([0 time_end]);
legend(lP1,lP5);

subplot(2,1,2);
set(gca,'FontSize',16);
grid on;
hold on
xlim([time_scale_G(1) time_scale_G(end-1)])
plot(time_scale_G(time0:end_time),vparallel_outputG(P1,time0:end_time),'b','LineWidth',2);
plot(time_scale_G(time0:end_time),vparallel_outputG(P5,time0:end_time),'Color',[0.2 0.8 0.2],'LineWidth',2);
xlabel('time (s)');
ylabel('v_{||} (m/s)');
xlim([0 time_end]);




%%
Ekin_outputG1=Ekin_outputG;
Ekin_outputG=Ekin_output2G;
figure(4);

subplot(2,1,1);
set(gca,'FontSize',16);
grid on;
hold on;


Ekin_error=Ekin_outputG*0;
Ekin_value=squeeze(mean(Ekin_outputG(:,TIME_INI,end),2));


for n=1:Nalphas_simulated
    Ekin_error(n,:)=(squeeze(Ekin_value(n))-Ekin_outputG(n,:))./squeeze(Ekin_value(n));
end

plot(time_scale_G(1:end_time),Ekin_outputG(P1,1:end_time),'b','LineWidth',2);
plot(time_scale_G(1:end_time),Ekin_outputG(P5,1:end_time),'Color',[0.2 0.8 0.2],'LineWidth',2);

legend(lP1,lP5);

xlim([time_scale_G(1) time_scale_G(end-1)])
plot(time_scale_G(1:end_time),Ekin_outputG(P1,1:end_time),'b','LineWidth',2);
plot([time_scale_G(1) time_scale_G(end_time)],squeeze(Ekin_value(P1))*[1 1],'k--','LineWidth',1);
plot(time_scale_G(1:end_time),Ekin_outputG(P1,1:end_time),'Color',[0.2 0.9 0.2],'LineWidth',2);
plot([time_scale_G(1) time_scale_G(end_time)],squeeze(Ekin_value(P5))*[1 1],'k--','LineWidth',1);
xlabel('time (s)');
ylabel('Ekin');
xlim([0 time_end]);


subplot(2,1,2);
set(gca,'FontSize',16);
grid on;
hold on
plot(time_scale_G(TIME_INI:end_time),Ekin_error(P1,TIME_INI:end_time),'b','LineWidth',1.5);
plot(time_scale_G(TIME_INI:end_time),Ekin_error(P5,TIME_INI:end_time),'Color',[0.2 0.9 0.2],'LineWidth',1.5);
xlabel('time (s)')
ylabel('\delta Ekin / Ekin');
% xlim([time_scale_G(1) time_scale_G(end-1)])
xlim([0 time_end]);

%

tr_error=0.5*(mean(mean(abs(Ekin_error(1:7,TIME_INI:end_time))))+mean(mean(abs(Ekin_error(25:32,TIME_INI:end_time)))));
passing_error=mean(mean(abs(Ekin_error(33:64,TIME_INI:end_time))));

disp('------------------------------------')
disp('TRAPPED particles Ekin error:')
tr_error
disp('PASSING particles Ekin error:')
passing_error
disp('------------------------------------')





%%
figure(5);

subplot(2,1,1);
set(gca,'FontSize',16);
grid on;
hold on;

pphi_recalc_outputG=(mHe/eV)*(Xpos_outputG(:,1:end_time)+R0).*(vphi_outputG(:,1:end_time))-ZHe*(psi_value_outputG(:,1:end_time));

pphi_outputG=squeeze(mean(pphi_recalc_outputG(:,TIME_INI,end),2));
pphi_error=pphi_recalc_outputG*0;

for n=1:size(pphi_recalc_outputG,1)
    pphi_error(n,:)=(squeeze(pphi_outputG(n))-pphi_recalc_outputG(n,:))./squeeze(pphi_outputG(n));
end

plot(time_scale_G(1:end_time),pphi_recalc_outputG(P1,1:end_time),'b','LineWidth',2);
plot(time_scale_G(1:end_time),pphi_recalc_outputG(P5,1:end_time),'Color',[0.2 0.8 0.2],'LineWidth',2);

legend(lP1,lP5);

xlim([time_scale_G(1) time_scale_G(end-1)])
plot(time_scale_G(1:end_time),(mHe/eV)*(Xpos_outputG(P1,1:end_time)+R0).*(vphi_outputG(P1,1:end_time))-ZHe*(psi_value_outputG(P1,1:end_time)),'b','LineWidth',2);
plot([time_scale_G(1) time_scale_G(end_time)],squeeze(pphi_outputG(P1))*[1 1],'k--','LineWidth',1);
plot(time_scale_G(1:end_time),(mHe/eV)*(Xpos_outputG(P5,1:end_time)+R0).*(vphi_outputG(P5,1:end_time))-ZHe*(psi_value_outputG(P5,1:end_time)),'Color',[0.2 0.9 0.2],'LineWidth',2);
plot([time_scale_G(1) time_scale_G(end_time)],squeeze(pphi_outputG(P5))*[1 1],'k--','LineWidth',1);
xlabel('time (s)');
ylabel('p_\phi');
xlim([0 time_end]);


subplot(2,1,2);
set(gca,'FontSize',16);
grid on;
hold on
plot(time_scale_G(time0:end_time),pphi_error(P1,time0:end_time),'b','LineWidth',1.5);
plot(time_scale_G(time0:end_time),pphi_error(P5,time0:end_time),'Color',[0.2 0.9 0.2],'LineWidth',1.5);
xlabel('time (s)')
ylabel('\delta p_\phi / p_\phi');
% xlim([time_scale_G(1) time_scale_G(end-1)])
xlim([0 time_end]);

%
tr_error2=0.5*(mean(mean(abs(pphi_error(1:7,TIME_INI:end))))+mean(mean(abs(pphi_error(25:32,TIME_INI:end)))));
passing_error2=mean(mean(abs(pphi_error(33:64,TIME_INI:end))));

disp('------------------------------------')
disp('TRAPPED particles pphi error:')
tr_error2
disp('PASSING particles pphi error:')
passing_error2
disp('------------------------------------')

%%
figure(6);

subplot(2,1,1);
set(gca,'FontSize',16);
grid on;
hold on;


mm_avg_outputG=squeeze(mean(mm_outputG(:,:),2));
mm_error=mm_outputG*0;

for n=1:size(pphi_recalc_outputG,1)
    mm_error(n,:)=(squeeze(mm_avg_outputG(n))-mm_outputG(n,:))./squeeze(mm_avg_outputG(n));
end

plot(time_scale_G(1:end_time),mm_outputG(P1,1:end_time),'b','LineWidth',2);
plot(time_scale_G(1:end_time),mm_outputG(P5,1:end_time),'Color',[0.2 0.8 0.2],'LineWidth',2);

legend(lP1,lP5);

xlim([time_scale_G(1) time_scale_G(end-1)])
plot(time_scale_G(time0:end_time),mm_outputG(P1,time0:end_time),'b','LineWidth',2);
plot([time_scale_G(1) time_scale_G(end_time)],squeeze(mm_avg_outputG(P1))*[1 1],'k--','LineWidth',1);
plot(time_scale_G(time0:end_time),mm_outputG(P5,time0:end_time),'Color',[0.2 0.9 0.2],'LineWidth',2);
plot([time_scale_G(1) time_scale_G(end_time)],squeeze(mm_avg_outputG(P5))*[1 1],'k--','LineWidth',1);
xlabel('time (s)');
ylabel('\mu (eV/T)');
xlim([0 time_end]);


subplot(2,1,2);
set(gca,'FontSize',16);
grid on;
hold on
plot(time_scale_G(time0:end_time),mm_error(P1,time0:end_time),'b','LineWidth',1.5);
plot(time_scale_G(time0:end_time),mm_error(P5,time0:end_time),'Color',[0.2 0.9 0.2],'LineWidth',1.5);
xlabel('time (s)')
ylabel('\delta \mu / \mu');
% xlim([time_scale_G(1) time_scale_G(end-1)])
xlim([0 time_end]);

%
tr_error3=0.5*(mean(mean(abs(mm_error(1:7,:))))+mean(mean(abs(mm_error(25:32,:)))));
passing_error3=mean(mean(abs(mm_error(33:64,:))));

% disp('------------------------------------')
% disp('TRAPPED particles mm error:')
% tr_error3
% disp('PASSING particles mm error:')
% passing_error3
% disp('------------------------------------')






figure(7);

subplot(3,1,1);
set(gca,'FontSize',16);
grid on;
hold on;



plot(time_scale_G(1:end_time),mm_outputG(P1,1:end_time),'b','LineWidth',2);
plot(time_scale_G(1:end_time),mm_outputG(P5,1:end_time),'Color',[0.2 0.8 0.2],'LineWidth',2);

legend(lP1,lP5);

xlim([time_scale_G(1) time_scale_G(end-1)])
plot(time_scale_G(time0:end_time),mm_outputG(P1,time0:end_time),'b','LineWidth',2);
plot([time_scale_G(1) time_scale_G(end_time)],squeeze(mm_avg_outputG(P1))*[1 1],'k--','LineWidth',1);
plot(time_scale_G(time0:end_time),mm_outputG(P5,time0:end_time),'Color',[0.2 0.9 0.2],'LineWidth',2);
plot([time_scale_G(1) time_scale_G(end_time)],squeeze(mm_avg_outputG(P5))*[1 1],'k--','LineWidth',1);
xlabel('time (s)');
ylabel('\mu (eV/T)');
xlim([0 time_end]);


subplot(3,1,2);
set(gca,'FontSize',16);
grid on;
hold on

set(gca,'FontSize',16);
grid on;
hold on
plot(time_scale_G(time0:end_time),pphi_error(P1,time0:end_time),'b','LineWidth',1.5);
plot(time_scale_G(time0:end_time),pphi_error(P5,time0:end_time),'Color',[0.2 0.9 0.2],'LineWidth',1.5);
xlabel('time (s)')
yl=ylabel('$\delta p_\varphi / p_\varphi$');
set(yl,'Interpreter','latex')
% xlim([time_scale_G(1) time_scale_G(end-1)])
xlim([0 time_end]);

subplot(3,1,3);
set(gca,'FontSize',16);
grid on;
hold on

plot(time_scale_G(time0:end_time),Ekin_error(P1,time0:end_time),'b','LineWidth',1.5);
plot(time_scale_G(time0:end_time),Ekin_error(P5,time0:end_time),'Color',[0.2 0.9 0.2],'LineWidth',1.5);
xlabel('time (s)')
yl=ylabel('$\delta \mathcal{E}_{kin} / \mathcal{E}_{kin}$');
set(yl,'Interpreter','latex')
% xlim([time_scale_G(1) time_scale_G(end-1)])
xlim([0 time_end]);


%
tr_error3=0.5*(mean(mean(abs(mm_error(1:7,:))))+mean(mean(abs(mm_error(25:32,:)))));
passing_error3=mean(mean(abs(mm_error(33:64,:))));


% disp('------------------------------------')
% disp('TRAPPED particles mm error:')
% tr_error3
% disp('PASSING particles mm error:')
% passing_error3
% disp('------------------------------------')

