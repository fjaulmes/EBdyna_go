rhoTAE=interp1(psi_scale,radial_r_value_flux,psiTAE)
rhoTAE=rhoTAE/max(radial_r_value_flux)

tmax=time_stamp-1
close all
% PARTS_EFF=find((mean(power_exchange_evol(:,:),1)>=1.2e8).*(alphas_Ekin'<=2.8e6));
% PARTS_EFF=find(((mean(power_exchange_evol(1:tmax,:),1))>=2.05e7));
PARTS_EFF=find((mean(power_exchange_evol(1:tmax,:),1)>=0.939e7).*(mean(vparallel_output(1:tmax,:),1)<8*1e6));

figure(2)
subplot(5,1,1)
set(gca,'fontsize',16)
grid on
hold on 
% plot(time_scale*omega_TAE/(2*pi),pphi_output(:,PARTS_EFF),'linewidth',2)
plot(time_scale*omega_TAE/(2*pi),vparallel_output(:,PARTS_EFF),'linewidth',2)
% plot((mHe/eV)*(R0+Xpos_output(:,PARTS_EFF)).*vphi_output(:,PARTS_EFF)-ZHe*psi_value_output(:,PARTS_EFF),'r--')
xlim([0 2.2])
ylabel('$$v_\parallel$$','Interpreter','latex')
% ylabel('$$p_\varphi$$','Interpreter','latex')

% subplot(4,1,2)
% grid on
% hold on 
% plot(Ekin_output(:,PARTS_EFF))
% xlim([0 tmax])
% 
% 
% figure(3)
radial_output=interp1(1:Nradial,radial_r_value_flux,psipos_output);
radial_output=radial_output/max(radial_r_value_flux);


subplot(5,1,2)
set(gca,'fontsize',16)
grid on
hold on 
plot(time_scale*omega_TAE/(2*pi),radial_output(:,PARTS_EFF),'linewidth',2)
xlim([0 2.2])
ylabel('\rho')
plot([0 2.2],[rhoTAE rhoTAE],'g--')
ylim([0.35 0.42])

subplot(5,1,3)
set(gca,'fontsize',16)
grid on
hold on 
plot(time_scale*omega_TAE/(2*pi),power_exchange_evol(:,PARTS_EFF),'linewidth',2)
xlim([0 2.2])
ylabel('P(W)')
ylim([-0.5 2]*1e7)

%%
% figure(4)
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
subplot(5,1,4)
set(gca,'fontsize',16)
grid on
hold on 
plot(time_scale*omega_TAE/(2*pi),Ekin_rel_output(:,PARTS_EFF),'linewidth',2)
xlim([0 2.2])
ylabel('$$\Delta \mathcal{E}_kin$$','Interpreter','latex')
ylim([-1 0.1]*1e3)

subplot(5,1,5)
set(gca,'fontsize',16)
grid on
hold on 

% plot(mm_rel_output(:,PARTS_EFF))
% 
% plot(Etot_rel_output(:,PARTS_EFF))
% plot(Ekin_rel_output(:,PARTS_EFF))
plot(time_scale*omega_TAE/(2*pi),Epot_output(:,PARTS_EFF),'linewidth',2)
% plot(Etot_output(:,PARTS_EFF)-(ZHe*Epot_output(:,PARTS_EFF)+Ekin_output(:,PARTS_EFF)))
xlim([0 2.2])

xlabel('\tau_{TAE}')
ylabel('\Phi (V)')

% ylim([-300 150])
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
