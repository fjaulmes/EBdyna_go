close all

tmax=time_stamp

pphi_output_recalc=(mHe/eV)*(R0+Xpos_part_output(:,:)).*vphi_output(:,:)-ZHe*psi_value_output(:,:);

algo_perf=(pphi_output(end,:)-pphi_output_recalc(end,:))./mean(pphi_output(end,:));
POP_LARGE=find((~alphas_ejected(1:length(algo_perf))).*(abs(algo_perf)>2e-4)');
length(POP_LARGE)

P1=ALL_TRAPPED(20);
P2=200;
PVPLL=find(vpll_avg<-1e7);
P3=PVPLL(5);
PVPLL=find(vpll_avg>1e7);
P4=PVPLL(5)


P1=POP_LARGE(1);
P2=POP_LARGE(2);
P3=POP_LARGE(3);
P4=POP_LARGE(4);

figure(2)
subplot(3,1,1)
grid on
hold on 
plot(Ekin_gc_output(:,P1))
plot(Ekin_output(:,P1),'r--')
plot(Etot_output(:,P1),'g--')
legend('v^2','Ekin output')
legend('gc','vtot')
% plot((mHe/eV)*(R0+Xpos_output(:,P1)).*vphi_output(:,P1)-ZHe*psi_value_output(:,P1),'r--')
xlim([0 tmax])

subplot(3,1,2)
grid on
hold on 
plot(pphi_output(:,P1))
plot((mHe/eV)*(R0+Xpos_part_output(:,P1)).*vphi_output(:,P1)-ZHe*psi_value_output(:,P1),'r--')
xlim([0 tmax])

subplot(3,1,3)
plot(vparallel_output(:,P1),'b')

grid on
hold on 
xlim([0 tmax])

%%

figure(3)
subplot(2,1,1)
grid on
hold on 
plot(Ekin_gc_output(:,P2))
plot(Ekin_output(:,P2),'r--')
legend('v^2','Ekin output')
% plot((mHe/eV)*(R0+Xpos_output(:,P2)).*vphi_output(:,P2)-ZHe*psi_value_output(:,P2),'r--')
xlim([0 tmax])

subplot(2,1,2)
grid on
hold on 
plot(pphi_output(:,P2))
plot((mHe/eV)*(R0+Xpos_part_output(:,P2)).*vphi_output(:,P2)-ZHe*psi_value_output(:,P2),'r--')
xlim([0 tmax])


figure(5)
subplot(2,1,1)
grid on
hold on 
plot(Ekin_gc_output(:,P3))
plot(Ekin_output(:,P3),'r--')
legend('v^2','Ekin output')
% plot((mHe/eV)*(R0+Xpos_output(:,P3)).*vphi_output(:,P3)-ZHe*psi_value_output(:,P3),'r--')
xlim([0 tmax])


subplot(2,1,2)
grid on
hold on 
plot(pphi_output(:,P3))
plot((mHe/eV)*(R0+Xpos_part_output(:,P3)).*vphi_output(:,P3)-ZHe*psi_value_output(:,P3),'r--')
xlim([0 tmax])


figure(6)
subplot(2,1,1)
grid on
hold on 
plot(Ekin_gc_output(:,P4))
plot(Ekin_output(:,P4),'r--')
legend('v^2','Ekin output')
% plot((mHe/eV)*(R0+Xpos_output(:,P4)).*vphi_output(:,P4)-ZHe*psi_value_output(:,P4),'r--')
xlim([0 tmax])

subplot(2,1,2)
grid on
hold on 
plot(pphi_output(:,P4))
plot((mHe/eV)*(R0+Xpos_part_output(:,P4)).*vphi_output(:,P4)-ZHe*psi_value_output(:,P4),'r--')
xlim([0 tmax])



%%
time_scale_norm=time_scale/(2*pi/(omega_TAE));

figure(8)
subplot(3,1,1)
set(gca,'fontsize',20)

grid on
hold on 
plot(time_scale_norm,(vparallel_output(:,P1)),'g')
plot(time_scale_norm,(vparallel_output(:,P2)),'r--')
plot(time_scale_norm,(vparallel_output(:,P3)),'b')
plot(time_scale_norm,(vparallel_output(:,P4)),'k--')

% plot(time_scale_norm,(Ekin_output(:,P1)-Ekin_gc_output(:,P1))./(Ekin_output(:,P4)),'b')
% plot(time_scale_norm,(Ekin_output(:,P2)-Ekin_gc_output(:,P2))./(Ekin_output(:,P4)),'r--')
% plot(time_scale_norm,(Ekin_output(:,P3)-Ekin_gc_output(:,P3))./(Ekin_output(:,P4)),'g')
% plot(time_scale_norm,(Ekin_output(:,P4)-Ekin_gc_output(:,P4))./(Ekin_output(:,P4)),'k--')

% plot((mHe/eV)*(R0+Xpos_output(:,P4)).*vphi_output(:,P4)-ZHe*psi_value_output(:,P4),'r--')
xlim([0 time_scale_norm(tmax)])

yl=ylabel('$v_\parallel$')
set(yl,'Interpreter','latex')


subplot(3,1,2)
set(gca,'fontsize',20)
grid on
hold on 

plot(time_scale_norm,(Ekin_output(:,P1))/(Ekin_output(1,P1)),'g')
plot(time_scale_norm,(Ekin_output(:,P2))/(Ekin_output(1,P2)),'r--')
plot(time_scale_norm,(Ekin_output(:,P3))/(Ekin_output(1,P3)),'b')
plot(time_scale_norm,(Ekin_output(:,P4))/(Ekin_output(1,P4)),'k--')
xlim([0 time_scale_norm(tmax)])

yl=ylabel('$\mathcal{E}_{\textrm{kin}} / \mathcal{E}_{\textrm{ini}}$')
set(yl,'Interpreter','latex')

subplot(3,1,3)
set(gca,'fontsize',20)

grid on
hold on 
plot(time_scale_norm,(pphi_output(:,P1)-pphi_output_recalc(:,P1))./(pphi_output(:,P4)),'g')
plot(time_scale_norm,(pphi_output(:,P2)-pphi_output_recalc(:,P2))./(pphi_output(:,P4)),'r--')
plot(time_scale_norm,(pphi_output(:,P3)-pphi_output_recalc(:,P3))./(pphi_output(:,P4)),'b')
plot(time_scale_norm,(pphi_output(:,P4)-pphi_output_recalc(:,P4))./(pphi_output(:,P4)),'k--')
xlim([0 time_scale_norm(tmax)])

yl=ylabel('$\delta p_\varphi / p_\varphi$')
set(yl,'Interpreter','latex')


xl=xlabel('t (2\pi / \omega_{TAE})')


%%

figure(9)
grid on
hold on 
set(gca,'fontsize',20)

plot(time_scale_norm,(Ekin_output(:,P1)-Ekin_gc_output(:,P1))./(Ekin_output(:,P4)),'b')
plot(time_scale_norm,(Ekin_output(:,P2)-Ekin_gc_output(:,P2))./(Ekin_output(:,P4)),'r--')
plot(time_scale_norm,(Ekin_output(:,P3)-Ekin_gc_output(:,P3))./(Ekin_output(:,P4)),'g')
plot(time_scale_norm,(Ekin_output(:,P4)-Ekin_gc_output(:,P4))./(Ekin_output(:,P4)),'k--')


yl=xlabel('$\delta \mathcal{E}_{\textrm{kin}} / \mathcal{E}_{\textrm{kin}}$')
set(yl,'Interpreter','latex')

xl=xlabel('$t (2 \pi / \omega_{TAE})$')
set(xl,'Interpreter','latex')


%%
figure(10)

subplot(2,1,1)
% grid on
% hold on 

set(gca,'fontsize',20)

algo_perf=(pphi_output(end,:)-pphi_output_recalc(end,:))./mean(pphi_output(end,:));
algo_perf=(pphi_output(end,:)-pphi_output_recalc(end,:));
algo_perf=min(algo_perf,0.1);
algo_perf=max(algo_perf,-0.1);
hist(algo_perf,(-0.000182:4e-6:0.000182));
xlim([-0.000184 0.000184])

% xl=xlabel('$\delta p_\varphi / \left< p_\varphi \right>$');
xl=xlabel('$\delta p_\varphi $');
set(xl,'Interpreter','latex')


subplot(2,1,2)
% grid on
% hold on 

set(gca,'fontsize',20)

algo_perf2=(Ekin_output(end,:)-Ekin_gc_output(end,:))./(Ekin_output(end,:));
algo_perf2(isnan(algo_perf2))=0;
algo_perf2(isinf(algo_perf2))=0;
algo_perf2=min(algo_perf2,0.1);
algo_perf2=max(algo_perf2,-0.1);
hist(algo_perf2,(-0.0002:16e-6:0.0002));
xlim([-0.0002 0.0002])

xl=xlabel('$\delta \mathcal{E}_{\textrm{kin}} / \mathcal{E}_{\textrm{kin}}$');
set(xl,'Interpreter','latex')



POP_LARGE=find((~alphas_ejected(1:length(algo_perf))).*(abs(algo_perf)>2e-4)');
length(POP_LARGE)