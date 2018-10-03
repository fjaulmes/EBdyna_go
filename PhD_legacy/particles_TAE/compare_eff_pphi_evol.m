

tmax=time_stamp
close all
% PARTS_EFF=find((mean(power_exchange_evol(:,:),1)>=1.2e8).*(alphas_Ekin'<=2.8e6));
PARTS_EFF=find(((mean(power_exchange_evol(1:tmax,:),1))>=1.2e6));

figure(2)
subplot(2,1,1)
grid on
hold on 
plot(pphi_output(:,PARTS_EFF),'b')
plot((mHe/eV)*(R0+Xpos_part_output(:,PARTS_EFF)).*vphi_output(:,PARTS_EFF)-ZHe*psi_value_output(:,PARTS_EFF),'r--')
xlim([0 tmax])

subplot(2,1,2)
grid on
hold on 
plot(Ekin_output(:,PARTS_EFF))
plot(Ekin_gc_output(:,PARTS_EFF),'r--')
xlim([0 tmax])


figure(3)
subplot(2,1,1)
grid on
hold on 
plot(vphi_output(:,PARTS_EFF))
xlim([0 tmax])

subplot(2,1,2)
grid on
hold on 
plot(power_exchange_evol(:,PARTS_EFF))
xlim([0 tmax])

%%
figure(4)
Etot_rel_output=Etot_output*0;
Ekin_rel_output=Etot_output*0;
mm_rel_output=Etot_output*0;
Epot_rel_output=Etot_output*0;

for n=1:length(PARTS_EFF)
    Etot_rel_output(:,PARTS_EFF(n))=Etot_output(:,PARTS_EFF(n))-Etot_output(3,PARTS_EFF(n));
    Ekin_rel_output(:,PARTS_EFF(n))=Ekin_output(:,PARTS_EFF(n))-Ekin_output(3,PARTS_EFF(n));
%     mm_rel_output(:,PARTS_EFF(n))=mm_output(:,PARTS_EFF(n))-mm_output(3,PARTS_EFF(n));
    Epot_rel_output(:,PARTS_EFF(n))=Epot_output(:,PARTS_EFF(n))-Epot_output(3,PARTS_EFF(n));
end
subplot(2,1,1)
grid on
hold on 
plot(Ekin_rel_output(:,PARTS_EFF))
xlim([0 tmax])

subplot(2,1,2)
grid on
hold on 

% plot(mm_rel_output(:,PARTS_EFF))
% 
% plot(Etot_rel_output(:,PARTS_EFF))
% plot(Ekin_rel_output(:,PARTS_EFF))
% plot(ZHe*Epot_rel_output(:,PARTS_EFF))
plot(Etot_output(:,PARTS_EFF)-(ZHe*Epot_output(:,PARTS_EFF)+Ekin_output(:,PARTS_EFF)))
xlim([0 tmax])


% figure(5)
% subplot(2,1,1)
% grid on
% hold on 
% plot(theta_output(:,PARTS_EFF))
% xlim([0 tmax])
% 
% subplot(2,1,2)
% grid on
% hold on 
% plot(Etot_output(:,PARTS_EFF))
% xlim([0 tmax])
