% load('../B_maps/B0001.mat')
% deltaB_B0=max(max(max(sqrt(BsX_map_phi.^2+BsZ_map_phi.^2+Bsphi_map_phi.^2))))/Bavg
time_stamp_go=length(WTAE_evol);
% list_pos =find(WTAE_evol==0);
% time_stamp_go=list_pos(1)-1

% [wmax time_stamp_go ]=max(WTAE_evol)
% time_stamp_go=length(WTAE_evol)*(time_stamp-1)/length(time_scale)
% time_stamp_go=time_stamp

figure(1)
hold on

tmax=time_stamp_go
tinit=2

plot(WTAE_AVG*(WTAE_evol-WTAE_evol(1)),'r');
% plot(Epart_tot_vD_rel_evol,'b--')
plot(Epart_tot_rel_evol,'k')
% plot((NB_PART_RESCALE*eV)*(sum(Ekin_output,2))-(NB_PART_RESCALE*eV)*(sum(Ekin_output(1,:),2)),'b')
% plot(Etot_part_evol-Etot_part_evol(tinit))
plot(Etot_part_evol-Etot_part_evol(tinit),'g--')
xlim([2 tmax])

%%
figure(2)
tinit=2


subplot(2,1,1)
hold on
% plot(Ekin_part_evol-Ekin_part_evol(1))
% plot(Epart_tot_rel_evol,'k--')
% plot(Etot_part_evol-Etot_part_evol(tinit))
plot(Etot_part_evol-Etot_part_evol(tinit),'b')
plot(Ekin_part_evol-Ekin_part_evol(tinit),'k--')
plot(Epart_tot_rel_evol,'r-.')
% plot(Epart_tot_vD_rel_evol,'r')
% plot(Epart_tot_rel_evol,'b--')
% plot(Ekin_part_evol-Ekin_part_evol(2))
% plot((NB_PART_RESCALE*eV)*(sum(Ekin_output,2)-sum(Ekin_output(4,:))),'r--')
% plot(Ekin_part_evol-Ekin_part_evol(tinit))
xlim([2 tmax])

subplot(2,1,2)
hold on
plot(Etot_part_evol-Etot_part_evol(tinit),'b')
plot(Etot_th_part_evol-Etot_th_part_evol(tinit),'g--')
xlim([2 tmax])

% plot(Epart_tot_rel_evol,'r')
% plot(WTAE_evol-WTAE_evol(10),'r');hold on
% plot((NB_PART_RESCALE*eV)*(sum(Ekin_output,2)-sum(Ekin_output(10,:),2)),'r')

%%
figure(3)
subplot(2,1,1)
hold on
plot((WTAE_evol)-(WTAE_evol(tinit)),'r');hold on
plot(Etot_part_evol-Etot_part_evol(tinit),'b')
plot(Etot_th_part_evol-Etot_th_part_evol(tinit),'g--')
xlim([2 tmax])

subplot(2,1,2)
hold on
plot(gamma_TAE_evol,'r');hold on
plot(gamma_TAE_vD_evol,'b--');hold on
xlim([2 tmax])

tmax=time_stamp

%%
figure(4)
time_scale_end=NB_OSCILLATIONS*2*pi/omega_TAE;
tinit=round(1*length(WTAE_evol)/NB_OSCILLATIONS)

t_lin_eval=time_stamp_go-100

ref_rough=(0:length(time_scale)-1)/(length(time_scale)-1);
ref_precise=(0:length(WTAE_evol)-1)/(length(WTAE_evol)-1);

time_scale_go=((1:length(WTAE_evol))-1)*time_scale_end/length(WTAE_evol);

% time_scale_go=interp1(ref_rough,time_scale,ref_precise());
plot(time_scale_go(tinit:t_lin_eval),log(100*sqrt(WTAE_evol(tinit:t_lin_eval))),'r');hold on




%%
figure(5)


tinit=1
time_scale_end=NB_OSCILLATIONS*2*pi/omega_TAE;

time_scale_go=((1:length(WTAE_evol))-1)*time_scale_end/length(WTAE_evol);
time_scale_go_norm=time_scale_go/(2*pi/omega_TAE);
subplot(3,1,1)
set(gca,'fontsize',22)

grid on
hold on
plot(time_scale_go_norm(tinit:end),WTAE_AVG*WTAE_evol(tinit:end),'r','linewidth',3);hold on
plot(time_scale_go_norm(tinit:end),Etot_part_evol(tinit:end)-Etot_part_evol(tinit),'b--','linewidth',2)
hl=legend('$W_{TAE} (\textrm{J})$','$\sum \mathcal{E}_{ions} (\textrm{J})$')
set(hl,'Interpreter','latex')
% xl=xlabel('$t \ (2\pi / \omega_{\textrm{TAE}})$')
% set(xl,'Interpreter','latex')
ylabel('J')

% plot(time_scale_go_norm(tinit:end),Etot_th_part_evol(tinit:end)-Etot_th_part_evol(tinit),'g--')
xlim([0 NB_OSCILLATIONS])

subplot(3,1,2)
grid on
hold on
set(gca,'fontsize',22)

plot(time_scale_go_norm(tinit:end),gamma_TAE_evol(tinit:end),'r','linewidth',2);hold on
% plot(time_scale_go(tinit:end),gamma_TAE_vD_evol(tinit:end),'b--');hold on
% xl=xlabel('$t \ (2\pi / \omega_{\textrm{TAE}})$')
% set(xl,'Interpreter','latex')
ylabel('\gamma_{TAE} (s^{-1})')

xlim([0 NB_OSCILLATIONS])
ylim([-2.6 4.6]*1e4)


subplot(3,1,3)
%
grid on
hold on
set(gca,'fontsize',22)

% lin_fit=0.37*time_scale_go_norm+3.1;
% lin_fit=0.31*time_scale_go_norm+3.3;
% lin_fit=0.25*time_scale_go_norm+2.9;
% lin_fit=0.19*time_scale_go_norm+3.7;
% lin_fit=0.051*time_scale_go_norm+3.6;

lin_fit=0.32*time_scale_go_norm+4.58;

plot(time_scale_go_norm(tinit:end),log(100*sqrt(WTAE_evol(tinit:end))),'r','linewidth',2);hold on
plot(time_scale_go_norm(tinit:end),lin_fit(tinit:end),'b--','linewidth',2);hold on
% plot(time_scale_go(tinit:end),gamma_TAE_vD_evol(tinit:end),'b--');hold on
xl=xlabel('$t \ (2\pi / \omega_{\textrm{TAE}})$')
set(xl,'Interpreter','latex')
yl=ylabel('$\textrm{ln} (\Phi _0)$')
set(yl,'Interpreter','latex')
xlim([0 NB_OSCILLATIONS])
ylim([4 12])

hl=legend('$\textrm{ln} (\Phi _0)$','linear fit')
set(hl,'Interpreter','latex')
