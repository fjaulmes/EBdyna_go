% load('particles_omega_spread_corr_Glisa140213_fc1h2.mat')
% load('../../data_tokamak/q_profile.mat', 'psi_rank_q1')


alphas_mm=alphas_mm_part;

phipos_outputG_wrap=wrap2pi(phipos_outputG);
omega_outputG=theta_outputG-phipos_outputG_wrap;
omega_outputG=wrap2pi(omega_outputG);

Ekin_avg=round(mean(Ekin_outputG(:,end))*0.001)   % in keV
r_output=sqrt((Xpos_outputG-X_axis).^2+Zpos_outputG.^2);
q_output=interp2(scale_X,scale_Z,q_initial_XZsmall_map',Xpos_outputG,Zpos_outputG,'*linear');

alphas_v=sqrt(2*(eV/mHe)*Ekin_avg);
r_ini=r_output(find(~alphas_ejected_G),end);
q_ini=q_output(find(~alphas_ejected_G),end);
% alphas_kappa=sqrt((Ekin_avg*R0+Bavg*alphas_mm.*(r_ini-R0))./(2*alphas_mm.*r_ini*Bavg));
% alphas_kappa=sqrt(((alphas_Ekin+Bavg*alphas_mm)*(R0-r_avg))./(2*alphas_mm.*r_avg*Bavg));
delta_orbit=0.5*(mHe/eV)*(alphas_v.*q_ini)./(ZHe*Bavg*sqrt(r_ini/R0));
delta_potato=R0*(2*mHe*q_ini.*alphas_v./(R0*ZHe*eV*Bavg)).^(2/3);


% Ekin_avg=200
NApsi=28
NL=16*NApsi;
close all
end_ts=(size(omega_outputG,2))-1
% end_ts=round(end_ts/3)
% end_ts=1800
alphas_omega0=omega_outputG(:,1);
alphas_omega=omega_outputG(:,end);
% alphas_omega0=omega_output(:,1);
% alphas_omega=omega_output(:,7000);


% figure(1);
% grid on;
% hold on;
% set(gca,'FontSize',16);
% 
% plot(alphas_omega0(1:NL),alphas_omega(1:NL),'b.');
% plot(alphas_omega0(7*NL+1:8*NL),alphas_omega(7*NL+1:8*NL),'r.');
% plot(alphas_omega0(4*NL+1:5*NL),alphas_omega(4*NL+1:5*NL),'g+')
% legend('counter passing','co passing','trapped')
% 
% plot(alphas_omega0(NL+1:2*NL),alphas_omega(NL+1:2*NL),'b.')
% plot(alphas_omega0(2*NL+1:3*NL),alphas_omega(2*NL+1:3*NL),'b+');
% plot(alphas_omega0(5*NL+1:6*NL),alphas_omega(5*NL+1:6*NL),'r+');
% plot(alphas_omega0(6*NL+1:7*NL),alphas_omega(6*NL+1:7*NL),'r.');
% 
% plot(alphas_omega0(3*NL+1:4*NL),alphas_omega(3*NL+1:4*NL),'g+');
% 
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



% alphas_psi0=interp1(psi_scale,1:257,psipos_outputG(:,1));
% alphas_psi=interp1(psi_scale,1:257,psipos_outputG(:,7));
alphas_psi0=psipos_outputG(:,1);
alphas_psi=psipos_outputG(:,end_ts);

figure(2);
grid on;
hold on;
set(gca,'FontSize',16);

plot(alphas_psi0(1:NL),alphas_psi(1:NL),'b.');
plot(alphas_psi0(7*NL+1:8*NL),alphas_psi(7*NL+1:8*NL),'r.');
plot(alphas_psi0(4*NL+1:5*NL),alphas_psi(4*NL+1:5*NL),'g+')
legend('counter passing','co passing','trapped')

plot(alphas_psi0(NL+1:2*NL),alphas_psi(NL+1:2*NL),'b.')
plot(alphas_psi0(2*NL+1:3*NL),alphas_psi(2*NL+1:3*NL),'b+');
plot(alphas_psi0(5*NL+1:6*NL),alphas_psi(5*NL+1:6*NL),'r+');
plot(alphas_psi0(6*NL+1:7*NL),alphas_psi(6*NL+1:7*NL),'r.');

plot(alphas_psi0(3*NL+1:4*NL),alphas_psi(3*NL+1:4*NL),'g+');

ylim([1 160])
xlim([1 160])
xlabel('\psi_{ini}')
ylabel('\psi_{final}')

plot([0 160],[size_r-2 size_r-2],'k--','LineWidth',2)
plot([0 160],[psi_rank_q1 psi_rank_q1],'k--','LineWidth',2)
plot([psi_rank_q1 psi_rank_q1],[0 160],'k--','LineWidth',2)


titre=strcat(['Redistribution in \psi of ' num2str(Ekin_avg)],' keV helium ions');
title(titre)

B0=Bavg
mm_outputG=Eperp_outputG./Bfield_outputG;
lambda_outputG=B0*mm_outputG./Ekin_outputG;

% end_ts=1287;
alphas_mm0=mm_outputG(:,1);
alphas_mm=mm_outputG(:,end_ts);

figure(6)
hold on
grid on
radial_bin_size=30;
psi_bin_pos=(11:radial_bin_size:160);
Npsi0=histc(alphas_psi0,psi_bin_pos-0.5*radial_bin_size);
Npsi=histc(alphas_psi,psi_bin_pos-0.5*radial_bin_size);
plot(psi_bin_pos,Npsi0,'b');
plot(psi_bin_pos,Npsi,'r');
legend('before collapse','after collapse')



figure(3);
grid on;
hold on;
set(gca,'FontSize',16);


plot(alphas_mm0(1:NL),alphas_mm(1:NL),'b.');
plot(alphas_mm0(7*NL+1:8*NL),alphas_mm(7*NL+1:8*NL),'r.');
plot(alphas_mm0(4*NL+1:5*NL),alphas_mm(4*NL+1:5*NL),'g+')
legend('counter passing','co passing','trapped')

plot(alphas_mm0(NL+1:2*NL),alphas_mm(NL+1:2*NL),'b.')
plot(alphas_mm0(2*NL+1:3*NL),alphas_mm(2*NL+1:3*NL),'b+');
plot(alphas_mm0(5*NL+1:6*NL),alphas_mm(5*NL+1:6*NL),'r+');
plot(alphas_mm0(6*NL+1:7*NL),alphas_mm(6*NL+1:7*NL),'r.');

plot(alphas_mm0(3*NL+1:4*NL),alphas_mm(3*NL+1:4*NL),'g+');

% ylim([0 2*pi])
% xlim([0 2*pi])
xlabel('\mu_{ini}')
ylabel('\mu_{final}')

plot([0 2*pi],0.5*[pi pi],'k--','LineWidth',2)
plot([0 2*pi],1.5*[pi pi],'k--','LineWidth',2)

titre=strcat(['Redistribution in \mu of ' num2str(Ekin_avg)],' keV helium ions');
title(titre)



figure(4);
grid on;
hold on;
set(gca,'FontSize',16);

% end_ts=1287;
alphas_lambda0=lambda_outputG(:,1);
alphas_lambda=lambda_outputG(:,end_ts);

plot(alphas_lambda0(1:NL),alphas_lambda(1:NL),'b.');
plot(alphas_lambda0(7*NL+1:8*NL),alphas_lambda(7*NL+1:8*NL),'r.');
plot(alphas_lambda0(4*NL+1:5*NL),alphas_lambda(4*NL+1:5*NL),'g+')
legend('counter passing','co passing','trapped')

plot(alphas_lambda0(NL+1:2*NL),alphas_lambda(NL+1:2*NL),'b.')
plot(alphas_lambda0(2*NL+1:3*NL),alphas_lambda(2*NL+1:3*NL),'b+');
plot(alphas_lambda0(5*NL+1:6*NL),alphas_lambda(5*NL+1:6*NL),'r+');
plot(alphas_lambda0(6*NL+1:7*NL),alphas_lambda(6*NL+1:7*NL),'r.');

plot(alphas_lambda0(3*NL+1:4*NL),alphas_lambda(3*NL+1:4*NL),'g+');


% ylim([0 2*pi])
% xlim([0 2*pi])
xlabel('\lambda_{ini}')
ylabel('\lambda_{final}')


titre=strcat(['Redistribution in pitch angle of ' num2str(Ekin_avg)],' keV helium ions');
title(titre)




figure(5);
grid on;
hold on;
set(gca,'FontSize',16);

 pphi_recalc_outputG=(mHe/eV)*(Xpos_outputG+R0).*(vphi_outputG)-ZHe*(psi_value_outputG);
alphas_pphi_ini=pphi_recalc_outputG(:,1);
alphas_pphi=pphi_recalc_outputG(:,end_ts);

plot(alphas_pphi_ini(1:NL),alphas_pphi(1:NL),'b.');
plot(alphas_pphi_ini(7*NL+1:8*NL),alphas_pphi(7*NL+1:8*NL),'r.');
plot(alphas_pphi_ini(4*NL+1:5*NL),alphas_pphi(4*NL+1:5*NL),'g+')
legend('counter passing','co passing','trapped')

plot(alphas_pphi_ini(NL+1:2*NL),alphas_pphi(NL+1:2*NL),'b.')
plot(alphas_pphi_ini(2*NL+1:3*NL),alphas_pphi(2*NL+1:3*NL),'b+');
plot(alphas_pphi_ini(5*NL+1:6*NL),alphas_pphi(5*NL+1:6*NL),'r+');
plot(alphas_pphi_ini(6*NL+1:7*NL),alphas_pphi(6*NL+1:7*NL),'r.');

plot(alphas_pphi_ini(3*NL+1:4*NL),alphas_pphi(3*NL+1:4*NL),'g+');


% ylim([0 2*pi])
% xlim([0 2*pi])
xlabel('p\phi_{ini}')
ylabel('p\phi_{final}')
plot([1 3.2],[1 3.2],'k--','LineWidth',2)
xlim([0 4])
ylim([0 4])

titre=strcat(['Redistribution in pphi of ' num2str(Ekin_avg)],' keV helium ions');
title(titre)




