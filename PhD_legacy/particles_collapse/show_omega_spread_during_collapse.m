% load('particles_omega_spread_corr_Glisa140213_fc1h2.mat')
phipos_outputG_wrap=wrap2pi(phipos_outputG);
omega_outputG=theta_outputG-phipos_outputG_wrap;
omega_outputG=wrap2pi(omega_outputG);
B0=3.5
mm_outputG=Eperp_outputG./Bfield_outputG;
lambda_outputG=B0*mm_outputG./Ekin_outputG;

Ekin_avg=800
NApsi=28
NL=16*NApsi;
close all
end_ts=(size(omega_outputG,2))-1



alphas_omega0=omega_outputG(:,1);
alphas_omega=omega_outputG(:,end);
% alphas_omega0=omega_output(:,1);
% alphas_omega=omega_output(:,7000);


figure(1);
grid on;
hold on;
set(gca,'FontSize',16);

plot(alphas_omega0(1:NL),alphas_omega(1:NL),'b.');
plot(alphas_omega0(5*NL+1:6*NL),alphas_omega(5*NL+1:6*NL),'r.');
plot(alphas_omega0(2*NL+1:3*NL),alphas_omega(2*NL+1:3*NL),'g+')
legend('counter passing','co passing','trapped')

plot(alphas_omega0(NL+1:2*NL),alphas_omega(NL+1:2*NL),'b+')

plot(alphas_omega0(3*NL+1:4*NL),alphas_omega(3*NL+1:4*NL),'g+')

plot(alphas_omega0(4*NL+1:5*NL),alphas_omega(4*NL+1:5*NL),'r+')

ylim([0 2*pi])
xlim([0 2*pi])
xlabel('\omega_{ini}')
ylabel('\omega_{final}')

plot([0 2*pi],0.5*[pi pi],'k--','LineWidth',2)
plot([0 2*pi],1.5*[pi pi],'k--','LineWidth',2)

titre=strcat(['Redistribution in \omega of ' num2str(Ekin_avg)],' keV helium ions');
title(titre)



% alphas_psi0=interp1(psi_scale,1:257,psipos_output(:,1));
% alphas_psi=interp1(psi_scale,1:257,psipos_output(:,7));
alphas_psi0=psipos_outputG(:,1);
alphas_psi=psipos_outputG(:,end);

figure(2);
grid on;
hold on;
set(gca,'FontSize',16);

plot(alphas_psi0(1:NL),alphas_psi(1:NL),'b.');
plot(alphas_psi0(5*NL+1:6*NL),alphas_psi(5*NL+1:6*NL),'r.');
plot(alphas_psi0(2*NL+1:3*NL),alphas_psi(2*NL+1:3*NL),'g+')
legend('counter passing','co passing','trapped')

plot(alphas_psi0(NL+1:2*NL),alphas_psi(NL+1:2*NL),'b+')
plot(alphas_psi0(4*NL+1:5*NL),alphas_psi(4*NL+1:5*NL),'r+')
% 
plot(alphas_psi0(3*NL+1:4*NL),alphas_psi(3*NL+1:4*NL),'g+')
% 

ylim([1 160])
xlim([1 160])
xlabel('\psi_{ini}')
ylabel('\psi_{final}')

plot([0 160],[size_r-2 size_r-2],'k--','LineWidth',2)
plot([0 160],[85 85],'k--','LineWidth',2)




titre=strcat(['Redistribution in \psi of ' num2str(Ekin_avg)],' keV helium ions');
title(titre)


alphas_mm0=mm_outputG(:,1);
alphas_mm=mm_outputG(:,end);
figure(3);
grid on;
hold on;
set(gca,'FontSize',16);


plot(alphas_mm0(1:NL),alphas_mm(1:NL),'b.');
plot(alphas_mm0(5*NL+1:6*NL),alphas_mm(5*NL+1:6*NL),'r.');
plot(alphas_mm0(2*NL+1:3*NL),alphas_mm(2*NL+1:3*NL),'g+')
legend('counter passing','co passing','trapped')

plot(alphas_mm0(NL+1:2*NL),alphas_mm(NL+1:2*NL),'b+')
plot(alphas_mm0(3*NL+1:4*NL),alphas_mm(3*NL+1:4*NL),'g+')

plot(alphas_mm0(4*NL+1:5*NL),alphas_mm(4*NL+1:5*NL),'r+')

% ylim([0 2*pi])
% xlim([0 2*pi])
xlabel('\mu_{ini}')
ylabel('\mu_{final}')

plot([0 2*pi],0.5*[pi pi],'k--','LineWidth',2)
plot([0 2*pi],1.5*[pi pi],'k--','LineWidth',2)

titre=strcat(['Redistribution in \mu of ' num2str(Ekin_avg)],' keV helium ions');
title(titre)





figure(5);
grid on;
hold on;
set(gca,'FontSize',16);

 pphi_recalc_outputG=(mHe/eV)*(Xpos_outputG+R0).*(vphi_outputG)-ZHe*(psi_value_outputG);
alphas_pphi_ini=pphi_recalc_outputG(:,1);
alphas_pphi=pphi_recalc_outputG(:,end_ts);

plot(alphas_pphi_ini(1:NL),alphas_pphi(1:NL),'b.');
plot(alphas_pphi_ini(5*NL+1:6*NL),alphas_pphi(5*NL+1:6*NL),'r.');
plot(alphas_pphi_ini(2*NL+1:3*NL),alphas_pphi(2*NL+1:3*NL),'g+')
legend('counter passing','co passing','trapped')

plot(alphas_pphi_ini(NL+1:2*NL),alphas_pphi(NL+1:2*NL),'b+')
plot(alphas_pphi_ini(4*NL+1:5*NL),alphas_pphi(4*NL+1:5*NL),'r+');

plot(alphas_pphi_ini(3*NL+1:4*NL),alphas_pphi(3*NL+1:4*NL),'g+');


% ylim([0 2*pi])
% xlim([0 2*pi])
xlabel('p\phi_{ini}')
ylabel('p\phi_{final}')
plot([1 3.2],[1 3.2],'k--','LineWidth',2)
xlim([1 3.2])
ylim([1 3.2])

titre=strcat(['Redistribution in pphi of ' num2str(Ekin_avg)],' keV helium ions');
title(titre)


% 
% 
% figure(1);
% grid on;
% hold on;
% set(gca,'FontSize',16);
% 
% plot(alphas_omega0(1:NL),alphas_omega(1:NL),'b.');
% plot(alphas_omega0(5*NL+1:6*NL),alphas_omega(5*NL+1:6*NL),'r.');
% legend('counter passing','co passing')
% 
% plot(alphas_omega0(NL+1:2*NL),alphas_omega(NL+1:2*NL),'b.')
% 
% plot(alphas_omega0(2*NL+1:3*NL),alphas_omega(2*NL+1:3*NL),'b.')
% plot(alphas_omega0(3*NL+1:4*NL),alphas_omega(3*NL+1:4*NL),'r.')
% 
% plot(alphas_omega0(4*NL+1:5*NL),alphas_omega(4*NL+1:5*NL),'r.')
% 
% ylim([0 2*pi])
% xlim([0 2*pi])
% xlabel('\omega_{ini}')
% ylabel('\omega_{final}')
% 
% plot([0 2*pi],0.5*[pi pi],'k--','LineWidth',2)
% plot([0 2*pi],1.5*[pi pi],'k--','LineWidth',2)
% 
% titre=strcat(['Redistribution in \omega of ' num2str(Ekin_avg)],' keV helium ions');
% title(titre)
% 
% 
% 
% % alphas_psi0=interp1(psi_scale,1:257,psipos_output(:,1));
% % alphas_psi=interp1(psi_scale,1:257,psipos_output(:,7));
% % alphas_psi0=psipos_outputG(:,1);
% % alphas_psi=psipos_outputG(:,7000);
% alphas_psi0=psipos_outputG(:,1);
% alphas_psi=psipos_outputG(:,end);
% 
% figure(2);
% grid on;
% hold on;
% set(gca,'FontSize',16);
% 
% plot(alphas_psi0(1:NL),alphas_psi(1:NL),'b.');
% plot(alphas_psi0(5*NL+1:6*NL),alphas_psi(5*NL+1:6*NL),'r.');
% legend('counter passing','co passing')
% 
% % plot(alphas_psi0(NL+1:2*NL),alphas_psi(NL+1:2*NL),'b.')
% % 
% plot(alphas_psi0(2*NL+1:3*NL),alphas_psi(2*NL+1:3*NL),'b.')
% plot(alphas_psi0(3*NL+1:4*NL),alphas_psi(3*NL+1:4*NL),'r.')
% % 
% % plot(alphas_psi0(4*NL+1:5*NL),alphas_psi(4*NL+1:5*NL),'r.')
% 
% ylim([1 150])
% xlim([1 150])
% xlabel('\psi_{ini}')
% ylabel('\psi_{final}')
% 
% plot([0 150],[133 133],'k--','LineWidth',2)
% 
% 
% titre=strcat(['Redistribution in \psi of ' num2str(Ekin_avg)],' keV helium ions');
% title(titre)
% 
% 
% alphas_mm0=mm_outputG(:,1);
% alphas_mm=mm_outputG(:,end);
% figure(3);
% grid on;
% hold on;
% set(gca,'FontSize',16);
% 
% 
% plot(alphas_mm0(1:NL),alphas_mm(1:NL),'b.');
% plot(alphas_mm0(5*NL+1:6*NL),alphas_mm(5*NL+1:6*NL),'r.');
% legend('counter passing','co passing')
% 
% plot(alphas_mm0(NL+1:2*NL),alphas_mm(NL+1:2*NL),'b.')
% plot(alphas_mm0(2*NL+1:3*NL),alphas_mm(2*NL+1:3*NL),'b.')
% 
% plot(alphas_mm0(4*NL+1:5*NL),alphas_mm(4*NL+1:5*NL),'r.')
% plot(alphas_mm0(3*NL+1:4*NL),alphas_mm(3*NL+1:4*NL),'r.')
% 
% % ylim([0 2*pi])
% % xlim([0 2*pi])
% xlabel('\mu_{ini}')
% ylabel('\mu_{final}')
% 
% plot([0 2*pi],0.5*[pi pi],'k--','LineWidth',2)
% plot([0 2*pi],1.5*[pi pi],'k--','LineWidth',2)
% 
% titre=strcat(['Redistribution in \mu of ' num2str(Ekin_avg)],' keV helium ions');
% title(titre)
