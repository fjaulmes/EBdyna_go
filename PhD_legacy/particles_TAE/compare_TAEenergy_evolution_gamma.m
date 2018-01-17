% load('../B_maps/B0001.mat')
% deltaB_B0=max(max(max(sqrt(BsX_map_phi.^2+BsZ_map_phi.^2+Bsphi_map_phi.^2))))/Bavg
time_stamp_go_recalc=length(WTAE_evol)*(time_stamp-1)/length(time_scale)

time_scale_go=(1:length(WTAE_evol))*time_scale(end)/length(WTAE_evol);
tmax=time_stamp
tinit=15

t_tau_max=time_scale_go(tmax)*omega_TAE/(2*pi)

figure(1)
subplot(3,1,1)
set(gca,'fontsize',16)
hold on
grid on

plot(time_scale_go*omega_TAE/(2*pi),WTAE_AVG*WTAE_evol-WTAE_AVG*WTAE_evol(1),'r','linewidth',3);
% plot(time_scale_go*omega_TAE/(2*pi),Epart_tot_vD_rel_evol,'b--')
plot(time_scale_go*omega_TAE/(2*pi),Epart_tot_rel_evol,'k--','linewidth',3)
% plot((NB_PART_RESCALE*eV)*(sum(Ekin_output,2))-(NB_PART_RESCALE*eV)*(sum(Ekin_output(1,:),2)),'b')
% plot(Etot_part_evol-Etot_part_evol(tinit))
% plot(time_scale_go*omega_TAE/(2*pi),Etot_th_part_evol-Etot_th_part_evol(tinit),'g--')
xlim([0 tmax])
hl=legend('$$W_{TAE}$$','$$\Sigma \mathcal{E}$$');
set(hl,'Interpreter','latex')
ylabel('J')
% xlabel('\tau _{TAE}')

subplot(3,1,2)
set(gca,'fontsize',16)
hold on
grid on
plot(time_scale_go*omega_TAE/(2*pi),100*sqrt(WTAE_AVG*WTAE_evol),'r','linewidth',3);hold on
xlim([0 tmax])
ylabel('\Phi_0 (V)')
% xlabel('\tau _{TAE}')

subplot(3,1,3)
set(gca,'fontsize',16)
hold on
grid on
plot(time_scale*omega_TAE/(2*pi),(eV*1e14)*(sum(pphi_output(:,:),2)-sum(pphi_output(1,:),2)),'b','linewidth',3)
xlim([0 tmax])
% ylim([-5 0.1]*1e-3)
ylabel('$$\Sigma p_\varphi$$ (Js)','Interpreter','latex')


[minvalue inipos]=min(abs(100*sqrt(WTAE_AVG*WTAE_evol)-120));
[minvalue endpos]=min(abs(100*sqrt(WTAE_AVG*WTAE_evol)-2850));
gamma_value=mean(gamma_TAE_evol(inipos:endpos))

subplot(3,1,2)
plot(time_scale_go(inipos:endpos)*omega_TAE/(2*pi),110*exp(gamma_value*(time_scale_go(inipos:endpos)-time_scale_go(inipos))),'k--','linewidth',3);hold on

% subplot(3,1,3)
% set(gca,'fontsize',16)
% grid on
% hold on
% 
% plot(time_scale_go*omega_TAE/(2*pi),gamma_TAE_evol,'r','linewidth',3);hold on
% plot(time_scale_go*omega_TAE/(2*pi),gamma_TAE_vD_evol,'b--','linewidth',3);hold on
% xlim([1 t_tau_max])

tmax=time_stamp;
xlabel('\tau _{TAE}')


save(SAVENAME,'-append','gamma_value')
