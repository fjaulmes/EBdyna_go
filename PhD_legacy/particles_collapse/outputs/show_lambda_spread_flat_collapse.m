close all

alphas_psi0=psi_pos0;
Ekin_avg=round(mean(alphas_Ekin)*1e-3)

CO_PASSING_POP=(alphas_kappa>1).*(vpll0>0);
COUNTER_PASSING_POP=(alphas_kappa>1).*(vpll0<0);
ALL_TRAPPED_POP=(alphas_kappa<=1);

CO_PASSING=find(CO_PASSING_POP.*(~alphas_ejected));
COUNTER_PASSING=find(COUNTER_PASSING_POP.*(~alphas_ejected));
ALL_TRAPPED=find(ALL_TRAPPED_POP.*(~alphas_ejected));


% figure(2);
% grid on;
% hold on;
% set(gca,'FontSize',16);
% 
% plot(alphas_psi0(COUNTER_PASSING),alphas_psi(COUNTER_PASSING),'b.');
% plot(alphas_psi0(CO_PASSING),alphas_psi(CO_PASSING),'r.');
% plot(alphas_psi0(ALL_TRAPPED),alphas_psi(ALL_TRAPPED),'g+')
% legend('counter passing','co passing','trapped')
% 
% ylim([1 160])
% xlim([1 160])
% xlabel('\psi_{ini}')
% ylabel('\psi_{final}')
% 
% plot([0 160],[size_r-2 size_r-2],'k--','LineWidth',2)
% plot([0 160],[psi_rank_q1 psi_rank_q1],'k--','LineWidth',2)
% plot([psi_rank_q1 psi_rank_q1],[0 160],'k--','LineWidth',2)
% plot([psi_core psi_core],[0 160],'r--','LineWidth',2)
% plot([psi_outer psi_outer],[0 160],'y--','LineWidth',2)
% 
% 
% titre=strcat(['Redistribution in \psi of ' num2str(Ekin_avg)],' keV helium ions');
% title(titre)


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


plot(alphas_mm0(COUNTER_PASSING),alphas_mm(COUNTER_PASSING),'b.');
plot(alphas_mm0(CO_PASSING),alphas_mm(CO_PASSING),'r.');
plot(alphas_mm0(ALL_TRAPPED),alphas_mm(ALL_TRAPPED),'g+')
legend('counter passing','co passing','trapped')



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


plot(alphas_lambda0(COUNTER_PASSING),alphas_lambda(COUNTER_PASSING),'b.');
plot(alphas_lambda0(CO_PASSING),alphas_lambda(CO_PASSING),'r.');
plot(alphas_lambda0(ALL_TRAPPED),alphas_lambda(ALL_TRAPPED),'g+')
legend('counter passing','co passing','trapped')


% ylim([0 2*pi])
% xlim([0 2*pi])
xlabel('\lambda_{ini}')
ylabel('\lambda_{final}')


titre=strcat(['Redistribution in pitch angle of ' num2str(Ekin_avg)],' keV helium ions');
title(titre)




figure(5);
grid on;
hold on;
set(gca,'FontSize',26);

alphas_pphi_ini=alphas_pphi0;

plot(alphas_pphi_ini(ALL_TRAPPED),alphas_pphi(ALL_TRAPPED),'g+')
plot(alphas_pphi_ini(COUNTER_PASSING),alphas_pphi(COUNTER_PASSING),'b.');
plot(alphas_pphi_ini(CO_PASSING),alphas_pphi(CO_PASSING),'r.');
legend('trapped','counter passing','co passing')




% ylim([0 2*pi])
% xlim([0 2*pi])
xl=xlabel('p$$_\varphi$$ini','interpreter','latex')
set(xl,'Interpreter','latex');
yl=xlabel('p$$_\varphi$$fin','interpreter','latex')
set(yl,'Interpreter','latex');
plot([-0.5 4],[-0.5 4],'k--','LineWidth',2)
axis equal
xlim([-0.5 4.0])
ylim([-0.5 4.0])

% titre=strcat(['Redistribution in pphi of ' num2str(Ekin_avg)],' keV helium ions');
% title(titre)




