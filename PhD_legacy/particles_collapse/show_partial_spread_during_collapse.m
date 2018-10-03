% load('particles_omega_spread_corr_Glisa140213_fc1h2.mat')

phipos_outputG_wrap=wrap2pi(phipos_outputG);
omega_outputG=theta_outputG-phipos_outputG_wrap;
omega_outputG=wrap2pi(omega_outputG);

alphas_psi0=psipos_outputG(:,1);
alphas_psi=psipos_outputG(:,end);

% ALPHAS_POPULATION=find(abs(alphas_lambda0-alphas_lambda)>0.2);

end_ts=1287;

trapped_particles_list=zeros(Nalphas_simulated,1);
for n=1:Nalphas_simulated
    if ~isempty(find(abs(vparallel_outputG(n,1:end_ts))<1e5))
        trapped_particles_list(n)=1;
    end
end
ALPHAS_POPULATION=find(~trapped_particles_list);
disp('Displaying data for sub-population of #alphas:')
length(ALPHAS_POPULATION)
alphas_pop=zeros(length(alphas_psi),1);
alphas_pop(ALPHAS_POPULATION)=1;

clear alphas_psi0 alphas_psi


Ekin_avg=1.6
NApsi=28
NL=16*NApsi;
close all

alphas_omega0=omega_outputG(:,1).*alphas_pop;
alphas_omega=omega_outputG(:,end_ts).*alphas_pop;
% alphas_omega0=omega_output(:,1);
% alphas_omega=omega_output(:,7000);


figure(1);
grid on;
hold on;
set(gca,'FontSize',16);

plot(alphas_omega0(1:NL),alphas_omega(1:NL),'b.');
plot(alphas_omega0(5*NL+1:6*NL),alphas_omega(5*NL+1:6*NL),'r.');
plot(alphas_omega0(2*NL+1:3*NL),alphas_omega(2*NL+1:3*NL),'g.')
legend('counter passing','co passing','trapped')

plot(alphas_omega0(NL+1:2*NL),alphas_omega(NL+1:2*NL),'b.')
plot(alphas_omega0(4*NL+1:5*NL),alphas_omega(4*NL+1:5*NL),'r.')

plot(alphas_omega0(3*NL+1:4*NL),alphas_omega(3*NL+1:4*NL),'g.')


ylim([0 2*pi])
xlim([0 2*pi])
xlabel('\omega_{ini}')
ylabel('\omega_{final}')

plot([0 2*pi],0.5*[pi pi],'k--','LineWidth',2)
plot([0 2*pi],1.5*[pi pi],'k--','LineWidth',2)

titre=strcat(['Redistribution in \omega of ' num2str(Ekin_avg)],' keV helium ions');
title(titre)



% alphas_psi0=interp1(psi_scale,1:257,psipos_outputG(:,1));
% alphas_psi=interp1(psi_scale,1:257,psipos_outputG(:,7));

figure(2);
grid on;
hold on;
set(gca,'FontSize',16);

alphas_psi0=psipos_outputG(:,1).*alphas_pop;
alphas_psi=psipos_outputG(:,end_ts).*alphas_pop;

plot(alphas_psi0(1:NL),alphas_psi(1:NL),'b.');
plot(alphas_psi0(5*NL+1:6*NL),alphas_psi(5*NL+1:6*NL),'r.');
plot(alphas_psi0(2*NL+1:3*NL),alphas_psi(2*NL+1:3*NL),'g.')
legend('counter passing','co passing','trapped')

plot(alphas_psi0(NL+1:2*NL),alphas_psi(NL+1:2*NL),'b.')
plot(alphas_psi0(4*NL+1:5*NL),alphas_psi(4*NL+1:5*NL),'r.')
% 
plot(alphas_psi0(3*NL+1:4*NL),alphas_psi(3*NL+1:4*NL),'g.')
% 

ylim([1 150])
xlim([1 150])
xlabel('\psi_{ini}')
ylabel('\psi_{final}')

plot([0 150],[133 133],'k--','LineWidth',2)


titre=strcat(['Redistribution in \psi of ' num2str(Ekin_avg)],' keV helium ions');
title(titre)

B0=3.5
mm_outputG=Eperp_outputG./Bfield_outputG;
lambda_outputG=B0*mm_outputG./Ekin_outputG;

alphas_mm0=mm_outputG(:,1).*alphas_pop;
alphas_mm=mm_outputG(:,end_ts).*alphas_pop;




figure(3);
grid on;
hold on;
set(gca,'FontSize',16);


plot(alphas_mm0(1:NL),alphas_mm(1:NL),'b.');
plot(alphas_mm0(5*NL+1:6*NL),alphas_mm(5*NL+1:6*NL),'r.');
plot(alphas_mm0(2*NL+1:3*NL),alphas_mm(2*NL+1:3*NL),'g.')
legend('counter passing','co passing','trapped')

plot(alphas_mm0(NL+1:2*NL),alphas_mm(NL+1:2*NL),'b.')
plot(alphas_mm0(3*NL+1:4*NL),alphas_mm(3*NL+1:4*NL),'g.')

plot(alphas_mm0(4*NL+1:5*NL),alphas_mm(4*NL+1:5*NL),'r.')

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

alphas_lambda0=lambda_outputG(:,1).*alphas_pop;
alphas_lambda=lambda_outputG(:,end_ts).*alphas_pop;

plot(alphas_lambda0(1:NL),alphas_lambda(1:NL),'b.');
plot(alphas_lambda0(5*NL+1:6*NL),alphas_lambda(5*NL+1:6*NL),'r.');
plot(alphas_lambda0(2*NL+1:3*NL),alphas_lambda(2*NL+1:3*NL),'g.')
legend('counter passing','co passing','trapped')

plot(alphas_lambda0(NL+1:2*NL),alphas_lambda(NL+1:2*NL),'b.')
plot(alphas_lambda0(3*NL+1:4*NL),alphas_lambda(3*NL+1:4*NL),'g.')

plot(alphas_lambda0(4*NL+1:5*NL),alphas_lambda(4*NL+1:5*NL),'r.')

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

end_ts=1287;
pphi_recalc_outputG=(mHe/eV)*(Xpos_outputG+R0).*(vphi_outputG)-ZHe*(psi_value_outputG);
alphas_pphi_ini=pphi_recalc_outputG(:,1).*alphas_pop;
alphas_pphi=pphi_recalc_outputG(:,end_ts).*alphas_pop;

plot(alphas_pphi_ini(1:NL),alphas_pphi(1:NL),'b.');
plot(alphas_pphi_ini(5*NL+1:6*NL),alphas_pphi(5*NL+1:6*NL),'r.');
plot(alphas_pphi_ini(2*NL+1:3*NL),alphas_pphi(2*NL+1:3*NL),'g.')
legend('counter passing','co passing','trapped')

plot(alphas_pphi_ini(NL+1:2*NL),alphas_pphi(NL+1:2*NL),'b.')
plot(alphas_pphi_ini(3*NL+1:4*NL),alphas_pphi(3*NL+1:4*NL),'g.')

plot(alphas_pphi_ini(4*NL+1:5*NL),alphas_pphi(4*NL+1:5*NL),'r.')

% ylim([0 2*pi])
% xlim([0 2*pi])
xlabel('p\phi_{ini}')
ylabel('p\phi_{final}')


titre=strcat(['Redistribution in pphi of ' num2str(Ekin_avg)],' keV helium ions');
title(titre)





