mHe=mD
ZHe=1
tmax=time_stamp
TINI=2

pphi_recalc_output=(mHe/eV)*(R0+Xpos_part_output(:,:)).*vphi_output(:,:)-ZHe*psi_value_output(:,:);
Dpphi_recalc=pphi_recalc_output(tmax,:)-pphi_recalc_output(TINI,:);
Dpphi=pphi_output(tmax,:)-pphi_output(TINI,:);
DE=Etot_output(tmax,:)-Etot_output(TINI,:);
DE_gc=Ekin_gc_output(tmax,:)-Ekin_gc_output(TINI,:);

% DE_th=Etot_th_output(tmax,:)-Etot_th_output(TINI,:);

DE_straight=(omega_TAE/nTAE)*Dpphi;

P1=1;
P2=10;
P3=20;
P4=30;

P5=100;
P6=150;
P7=200;
P8=300;

figure(1)
subplot(2,1,1)
grid on
hold on 
plot(pphi_output(:,P1))
plot(pphi_recalc_output(:,P1),'r--')
xlim([0 tmax])

subplot(2,1,2)
grid on
hold on 
plot(pphi_output(:,P2))
plot(pphi_recalc_output(:,P2),'r--')
xlim([0 tmax])



figure(2)
subplot(2,1,1)
grid on
hold on 
plot(pphi_output(:,P3))
plot(pphi_recalc_output(:,P3),'r--')
xlim([0 tmax])

subplot(2,1,2)
grid on
hold on 
plot(pphi_output(:,P4))
plot(pphi_recalc_output(:,P4),'r--')
xlim([0 tmax])


figure(3)
subplot(2,1,1)
grid on
hold on 
plot(pphi_output(:,P5))
plot(pphi_recalc_output(:,P5),'r--')
xlim([0 tmax])

subplot(2,1,2)
grid on
hold on 
plot(pphi_output(:,P6))
plot(pphi_recalc_output(:,P6),'r--')
xlim([0 tmax])


figure(4)
subplot(2,1,1)
grid on
hold on 
plot(pphi_output(:,P7))
plot(pphi_recalc_output(:,P7),'r--')
xlim([0 tmax])

subplot(2,1,2)
grid on
hold on 
plot(pphi_output(:,P8))
plot(pphi_recalc_output(:,P8),'r--')
xlim([0 tmax])

%%
POP_NON_EJ=find(~(alphas_ejected));
figure(8)
hold on
plot(Dpphi(POP_NON_EJ),DE_gc(POP_NON_EJ),'g.')
plot(Dpphi_recalc(POP_NON_EJ),DE(POP_NON_EJ),'b.')
plot(Dpphi(POP_NON_EJ),DE(POP_NON_EJ),'r.')
plot([-0.1 0.1],(omega_TAE/nTAE)*[-0.1 0.1],'k--')


%%
figure(7)

hold on

plot(mean(vparallel_output(1:tmax,POP_NON_EJ),1),DE(POP_NON_EJ),'b.')
Epll_output=0.5*(mHe/eV)*vparallel_output.^2;
Eperp_output=max(Ekin_output-Epll_output,0);
vperp_output=sqrt(2*(eV/mHe)*Eperp_output);

plot([vA_TAE/5 vA_TAE/5],[-4.5 4.5]*1e4,'g--','linewidth',3)
plot([vA3_TAE vA3_TAE],[-4.5 4.5]*1e4,'g--','linewidth',3)
plot([vA_TAE vA_TAE],[-4.5 4.5]*1e4,'g--','linewidth',3)
plot(-[vA_TAE/5 vA_TAE/5],[-4.5 4.5]*1e4,'g--','linewidth',3)
plot(-[vA3_TAE vA3_TAE],[-4.5 4.5]*1e4,'g--','linewidth',3)
plot(-[vA_TAE vA_TAE],[-4.5 4.5]*1e4,'g--','linewidth',3)



figure(6)
plot(mean(Ekin_output(1:tmax,POP_NON_EJ),1),DE(POP_NON_EJ),'b.')

%%
figure(5)
hold on

plot(mean(vperp_output(1:tmax,POP_NON_EJ),1),DE(POP_NON_EJ),'b.')
plot([vA_TAE/5 vA_TAE/5],[-4.5 4.5]*1e4,'g--','linewidth',3)
plot([vA3_TAE vA3_TAE],[-4.5 4.5]*1e4,'g--','linewidth',3)
plot([vA_TAE vA_TAE],[-4.5 4.5]*1e4,'g--','linewidth',3)
