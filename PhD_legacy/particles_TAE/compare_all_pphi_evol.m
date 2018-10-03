% tmax=600
tmax=time_stamp

pphi_recalc_output=(mHe/eV)*(R0+Xpos_part_output(:,:)).*vphi_output(:,:)-ZHe*psi_value_output(:,:);


figure(2)
grid on
hold on 
plot(sum(pphi_output(:,:),2))
plot(sum(pphi_recalc_output,2),'r--')
xlim([0 tmax])


LARGE_DEVIATION=find(sum(abs(pphi_recalc_output-pphi_output),1)>1e-3);
PART_NUMBER=4;

figure(3)
subplot(2,1,1)
plot(Ekin_output(1:tmax,LARGE_DEVIATION(PART_NUMBER)))

subplot(2,1,2)
hold on
plot(pphi_output(1:tmax,LARGE_DEVIATION(PART_NUMBER)))
plot(pphi_recalc_output(1:tmax,LARGE_DEVIATION(PART_NUMBER)),'r--')



figure(4)
subplot(2,1,1)
plot(vparallel_output(1:tmax,LARGE_DEVIATION(PART_NUMBER)))

subplot(2,1,2)
hold on
plot(Eperp_output(1:tmax,LARGE_DEVIATION(PART_NUMBER)))


figure(5)
subplot(2,1,1)
hold on
plot(Etot_output(1:tmax,LARGE_DEVIATION(PART_NUMBER)))
plot(Ekin_output(1:tmax,LARGE_DEVIATION(PART_NUMBER)),'r--')

subplot(2,1,2)
hold on
plot(psi_star_output(1:tmax,LARGE_DEVIATION(PART_NUMBER)))


% figure(5)
% plot(Xpos_output(1:tmax,LARGE_DEVIATION(2)),Zpos_output(1:tmax,LARGE_DEVIATION(2)))
% 
% 
