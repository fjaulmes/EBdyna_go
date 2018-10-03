Ecrit=60

figure(2)
subplot(2,2,1)
set(gca,'FontSize',26);
% title('q_0=0.94')
hold on
grid on
% plot(Ekin_bins,Dpphi_Ekin_co_pos(1:end)./Dpphi_Ekin_counter_pos(1:end),'r','LineWidth',4)
% plot(Ekin_bins,Dpphi_Ekin_co_neg(1:end)./Dpphi_Ekin_counter_neg(1:end),'b--','LineWidth',4)
% plot(Ekin_bins,Dpphi_Ekin_counter_pos(1:end),'b--','LineWidth',4)
% plot(Ekin_bins,Dpphi_Ekin_counter_neg(1:end),'b--','LineWidth',4)
% plot(Ekin_bins,Dpphi_Ekin_co_pos(1:end),'r','LineWidth',4)
% plot(Ekin_bins,Dpphi_Ekin_co_neg(1:end),'r','LineWidth',4)
% plot(Ekin_bins,Dpphi_Ekin_counter_pos(1:end),'b--','LineWidth',4)
% plot(Ekin_bins,Dpphi_Ekin_counter_neg(1:end),'b--','LineWidth',4)
% plot([Ecrit Ecrit],[-1 1],'--','color',[0.1 0.7 0.1],'LineWidth',6);
plot([Ecrit Ecrit],[-1 1],'--','color',[0.1 0.7 0.1],'LineWidth',6);
% plot(Ekin_bins,Dpphi_Ekin_co(1:end)/Dpphi_Ekin_co(1),'r','LineWidth',4)
% plot(Ekin_bins,Dpphi_Ekin_counter(1:end)/Dpphi_Ekin_counter(1),'b--','LineWidth',4)
plot(Ekin_bins,Dpphi_Ekin_co(1:end),'r','LineWidth',4)
plot(Ekin_bins,Dpphi_Ekin_counter(1:end),'b--','LineWidth',4)

% plot(Ekin_bins,Dpphi_Ekin_co_pos(1:end)/Dpphi_Ekin_co_pos(1),'r','LineWidth',4)
% plot(Ekin_bins,-Dpphi_Ekin_co_neg(1:end)/Dpphi_Ekin_co_neg(1),'r','LineWidth',4)
% plot(Ekin_bins,Dpphi_Ekin_counter_pos(1:end)/Dpphi_Ekin_counter_pos(1),'b--','LineWidth',4)
% plot(Ekin_bins,-Dpphi_Ekin_counter_neg(1:end)/Dpphi_Ekin_counter_neg(1),'b--','LineWidth',4)
% plot([Ebin1value 400 800 2200 3200],Dpphi_Ekin_counter_neg(1:end),'b-.','LineWidth',4)
% plot([Ebin1value 400 800 2200 3200],Dpphi_Ekin_co_neg(1:end),'r','LineWidth',4)
% set(gca,'XTick',Ekin_v_bins); % Change x-axis ticks
% set(gca,'XTickLabel',Ekin_bins);
% xlim([200 2300])
% xlabel('Ekin');
% yl=ylabel('$$|\tilde{\Delta p_\varphi }|$$','interpreter','latex');
yl=ylabel('$$|\Delta p_\varphi |$$','interpreter','latex');
set(yl,'interpreter','latex')
% ylim([0 0.23])
% ylim([-0.23 0.23])

% figure(2)
% subplot(4,1,2)
% set(gca,'FontSize',26);
% hold on
% grid on
% plot([200 400 800 2200 3200],dr_Ekin_counter(2:end),'b-.','LineWidth',4)
% plot([200 400 800 2200 3200],dr_Ekin_co(2:end),'r','LineWidth',4)
% xlim([200 3200])
% % xlabel('Ekin');
% ylabel('\deltar (m)');

figure(2)
subplot(2,2,3)
set(gca,'FontSize',26);
hold on
grid on

% plot(Ekin_bins,g_Ekin_counter1(1:end),'b-.','LineWidth',4)
% plot(Ekin_bins,g_Ekin_co1(1:end),'r','LineWidth',4)
plot(Ekin_bins,g_Ekin_counter(1:end),'b-.','LineWidth',4)
plot(Ekin_bins,g_Ekin_co(1:end),'r','LineWidth',4)

plot([Ecrit Ecrit],[-1 2.0],'--','color',[0.1 0.7 0.1],'LineWidth',6);

plot(Ekin_bins,g_Ekin_counter(1:end),'b-.','LineWidth',4)
plot(Ekin_bins,g_Ekin_co(1:end),'r','LineWidth',4)

% plot(Ekin_bins,g_Ekin_counter1(1:end),'b-.','LineWidth',4)
% plot(Ekin_bins,g_Ekin_co1(1:end),'r','LineWidth',4)
% plot(Ekin_bins,g_Ekin_counter2(1:end),'b-.','LineWidth',4)
% plot(Ekin_bins,g_Ekin_co2(1:end),'r','LineWidth',4)

% set(gca,'XTick',Ekin_v_bins); % Change x-axis ticks
% set(gca,'XTickLabel',Ekin_bins);

% xlim([200 2300])
% ylim([0 0.9])
% xlabel('Ekin (keV)');
ylabel('|\omega_{v_D}|/|\omega_{\psi}|');
hl=legend('counter-passing','co-passing')
set(hl,'fontsize',20)
xlabel('E_{kin} (keV)');

annotation('textbox',...
    [0.245 0.41 0.0483 0.02],...
    'EdgeColor','none',...
    'String',{'E_{cc}'},'fontsize',16,...
    'FitBoxToText','off','color',[0 0.4 0.1]);

annotation('textbox',...
    [0.245 0.86 0.0483 0.02],...
    'EdgeColor','none',...
    'String',{'E_{cc}'},'fontsize',16,...
    'FitBoxToText','off','color',[0 0.4 0.1]);

% annotation('textbox',...
%     [0.22 0.61 0.0483 0.02],...
%     'EdgeColor','none',...
%     'String',{'(1)'},'fontsize',18,...
%     'FitBoxToText','off','color',[0 0 0]);
% 
% annotation('textbox',...
%     [0.22 0.87 0.0483 0.02],...
%     'EdgeColor','none',...
%     'String',{'(2)'},'fontsize',18,...
%     'FitBoxToText','off','color',[0 0 0]);
% 
% 
