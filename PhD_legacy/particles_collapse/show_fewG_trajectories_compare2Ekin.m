close all

load('fewG_Ekin5_collapse_local280113_fc2h2.mat')
P1=1356
mm_outputG1=Eperp_outputG./Bfield_outputG;
Xpos_outputG1=Xpos_outputG;
Zpos_outputG1=Zpos_outputG;
Etot_outputG1=Etot_outputG;
Epot_outputG1=Epot_outputG;
Eperp_outputG1=Eperp_outputG;
vphi_outputG1=vphi_outputG;
vparallel_outputG1=vparallel_outputG;
psi_value_outputG1=psi_value_outputG;
pphi_outputG1=pphi_outputG;
psi_star_outputG1=psi_star_outputG;

% time_stamp_end = 1231
% find(Ekin_outputG(:,1)<60*1e3)
load('../data_tokamak/q_profile.mat', 'psi_rank_q1')
end_ts=size(Xpos_outputG,2)
Nalphas_simulated=size(Xpos_outputG,1)
omega_cr=pi/(1.44*1e-4)
 
NAPsi=28;

pphi_recalc1=(mHe/eV)*(Xpos_outputG(:,1:end)+R0).*(vphi_outputG(:,1:end))-ZHe*(psi_value_outputG(:,1:end));

% for (n=1:Nalphas_simulated)
% pphi_outputG(n,:)=pphi_outputG(n,:)-pphi_outputG(n,1)+pphi_recalc(n,1);
% end
calculate_vD_XZ_map;
calculate_gradZ_theta_map;
calculate_gradX_theta_map;

vDZ_output_corr1=interp2(scale_X,scale_Z,vD_Z_XZ_map',Xpos_outputG,Zpos_outputG,'*linear');
vDX_output_corr1=interp2(scale_X,scale_Z,vD_X_XZ_map',Xpos_outputG,Zpos_outputG,'*linear');
vDphi_output_corr1=interp2(scale_X,scale_Z,vD_phi_XZ_map',Xpos_outputG,Zpos_outputG,'*linear');
GZ_theta_output_pop1=interp2(scale_X,scale_Z,gradZ_theta_map',Xpos_outputG,Zpos_outputG,'*linear');
GX_theta_output_pop1=interp2(scale_X,scale_Z,gradX_theta_map',Xpos_outputG,Zpos_outputG,'*linear');
B_output_pop1=interp2(scale_X,scale_Z,Btot_XZ_map',Xpos_outputG,Zpos_outputG,'*linear');


for ts=1:end_ts
    vDZ_output_corr1(:,ts)=(vDZ_output_corr1(:,ts)).*(2*Ekin_outputG(:,ts)-B_output_pop1(:,ts).*mm_outputG1(:,ts))/(ZHe);
    vDX_output_corr1(:,ts)=(vDX_output_corr1(:,ts)).*(2*Ekin_outputG(:,ts)-B_output_pop1(:,ts).*mm_outputG1(:,ts))/(ZHe);
    vDphi_output_corr1(:,ts)=(vDphi_output_corr1(:,ts)).*(2*Ekin_outputG(:,ts)-B_output_pop1(:,ts).*mm_outputG1(:,ts))/(ZHe);
end


bX_output_pop1=interp2(scale_X,scale_Z,bX_XZ_map',Xpos_outputG,Zpos_outputG,'*linear');
bZ_output_pop1=interp2(scale_X,scale_Z,bZ_XZ_map',Xpos_outputG,Zpos_outputG,'*linear');
bphi_output_pop1=sqrt(1-(bX_output_pop1.^2+bZ_output_pop1.^2));

for number=1:Nalphas_simulated
%     number=ALL_TRAPPED(n);
    omega_precess_avg1(number,:)=(GX_theta_output_pop1(number,:).*vDX_output_corr1(number,:)+GZ_theta_output_pop1(number,:).*vDZ_output_corr1(number,:)-...
       vDphi_output_corr1(number,:)./(R0+Xpos_outputG(number,:)));
    omega_psi_avg1(number,:)=(GZ_theta_output_pop1(number,:).*bZ_output_pop1(number,:).*vparallel_outputG(number,:)-...
       bphi_output_pop1(number,:).*vparallel_outputG(number,:)./(R0+Xpos_outputG(number,:)));
end

% vperp_sq_outputG=vDsq_outputG;
% Bfield_part_outputG=interp2(scale_X,scale_Z,Btot_XZ_map',Xpos_outputG,Zpos_outputG,'*linear');
mm_outputG1=Eperp_outputG./Bfield_outputG;
lambda_outputG1=Bavg*mm_outputG1./Ekin_outputG;
% psi_value_outputG=interp2(scale_X,scale_Z,psi_XZsmall_map',Xpos_outputG,Zpos_outputG,'*linear');
% lap_psi_outputG=interp2(scale_X,scale_Z,lap_psi_XZsmall_map',Xpos_outputG,Zpos_outputG,'*linear');

% for (n=1:24)
% psi_value_corr_outputG(n,:)=psipos_outputG(n,:)-alphas_psi_star_ini(n)+psi_star_outputG(n,:);
% end
lambda_P1=lambda_outputG1(P1,1);


Ekin_P1=round(Ekin_outputG(P1,1));


lP1=strcat('\lambda _{ini} = ',num2str(lambda_P1));


EP1=strcat('Ekin _{ini} = ',num2str(Ekin_P1));




phipos_outputG_wrap1=wrap2pi(phipos_outputG);
omega_outputG1=theta_outputG-phipos_outputG_wrap1;
omega_outputG1=wrap2pi(omega_outputG1);










load('fewG_Ekin90_collapse_local290113_fc2h2.mat')
P2=1356
mm_outputG2=Eperp_outputG./Bfield_outputG;
Xpos_outputG2=Xpos_outputG;
Zpos_outputG2=Zpos_outputG;
Etot_outputG2=Etot_outputG;
Epot_outputG2=Epot_outputG;
Eperp_outputG2=Eperp_outputG;
vphi_outputG2=vphi_outputG;
vparallel_outputG2=vparallel_outputG;
psi_value_outputG2=psi_value_outputG;
pphi_outputG2=pphi_outputG;
psi_star_outputG2=psi_star_outputG;

% time_stamp_end = 1231
% find(Ekin_outputG(:,1)<60*1e3)

end_ts=size(Xpos_outputG,2)
Nalphas_simulated=size(Xpos_outputG,1)
omega_cr=pi/(1.44*1e-4)
 
NAPsi=28;

pphi_recalc1=(mHe/eV)*(Xpos_outputG(:,1:end)+R0).*(vphi_outputG(:,1:end))-ZHe*(psi_value_outputG(:,1:end));

% for (n=1:Nalphas_simulated)
% pphi_outputG(n,:)=pphi_outputG(n,:)-pphi_outputG(n,1)+pphi_recalc(n,1);
% end

vDZ_output_corr2=interp2(scale_X,scale_Z,vD_Z_XZ_map',Xpos_outputG,Zpos_outputG,'*linear');
vDX_output_corr2=interp2(scale_X,scale_Z,vD_X_XZ_map',Xpos_outputG,Zpos_outputG,'*linear');
vDphi_output_corr2=interp2(scale_X,scale_Z,vD_phi_XZ_map',Xpos_outputG,Zpos_outputG,'*linear');
GZ_theta_output_pop2=interp2(scale_X,scale_Z,gradZ_theta_map',Xpos_outputG,Zpos_outputG,'*linear');
GX_theta_output_pop2=interp2(scale_X,scale_Z,gradX_theta_map',Xpos_outputG,Zpos_outputG,'*linear');
B_output_pop2=interp2(scale_X,scale_Z,Btot_XZ_map',Xpos_outputG,Zpos_outputG,'*linear');


for ts=1:end_ts
    vDZ_output_corr2(:,ts)=(vDZ_output_corr2(:,ts)).*(2*Ekin_outputG(:,ts)-B_output_pop2(:,ts).*mm_outputG2(:,ts))/(ZHe);
    vDX_output_corr2(:,ts)=(vDX_output_corr2(:,ts)).*(2*Ekin_outputG(:,ts)-B_output_pop1(:,ts).*mm_outputG1(:,ts))/(ZHe);
    vDphi_output_corr2(:,ts)=(vDphi_output_corr2(:,ts)).*(2*Ekin_outputG(:,ts)-B_output_pop2(:,ts).*mm_outputG2(:,ts))/(ZHe);
end


bX_output_pop2=interp2(scale_X,scale_Z,bX_XZ_map',Xpos_outputG,Zpos_outputG,'*linear');
bZ_output_pop2=interp2(scale_X,scale_Z,bZ_XZ_map',Xpos_outputG,Zpos_outputG,'*linear');
bphi_output_pop2=sqrt(1-(bX_output_pop2.^2+bZ_output_pop2.^2));

for number=1:Nalphas_simulated
%     number=ALL_TRAPPED(n);
    omega_precess_avg2(number,:)=(GX_theta_output_pop2(number,:).*vDX_output_corr2(number,:)+GZ_theta_output_pop2(number,:).*vDZ_output_corr2(number,:)-...
       vDphi_output_corr2(number,:)./(R0+Xpos_outputG(number,:)));
    omega_psi_avg2(number,:)=(GZ_theta_output_pop2(number,:).*bZ_output_pop2(number,:).*vparallel_outputG(number,:)-...
       bphi_output_pop1(number,:).*vparallel_outputG(number,:)./(R0+Xpos_outputG(number,:)));
end

% vperp_sq_outputG=vDsq_outputG;
% Bfield_part_outputG=interp2(scale_X,scale_Z,Btot_XZ_map',Xpos_outputG,Zpos_outputG,'*linear');
mm_outputG2=Eperp_outputG./Bfield_outputG;
lambda_outputG2=Bavg*mm_outputG2./Ekin_outputG;
% psi_value_outputG=interp2(scale_X,scale_Z,psi_XZsmall_map',Xpos_outputG,Zpos_outputG,'*linear');
% lap_psi_outputG=interp2(scale_X,scale_Z,lap_psi_XZsmall_map',Xpos_outputG,Zpos_outputG,'*linear');

% for (n=1:24)
% psi_value_corr_outputG(n,:)=psipos_outputG(n,:)-alphas_psi_star_ini(n)+psi_star_outputG(n,:);
% end

lambda_P2=lambda_outputG2(P2,1);


Ekin_P2=round(Ekin_outputG(P2,1));

lP2=strcat('\lambda _{ini} = ',num2str(lambda_P2));

EP2=strcat('Ekin _{ini} = ',num2str(Ekin_P2));



phipos_outputG_wrap2=wrap2pi(phipos_outputG);
omega_outputG2=theta_outputG-phipos_outputG_wrap2;
omega_outputG2=wrap2pi(omega_outputG2);












size_phi=257;
scale_tor=40*pi*((0:size_phi-1)/(size_phi-1));
Btot_X_vector=mean(Btot_XZ_map,2);
[Btot_Xphi_map tmp]=meshgrid(Btot_X_vector,ones(size_phi,1));
Btot_Xphi_map=Btot_Xphi_map';

end_time=size(time_scale_G,2);
ejection_time=size(time_scale_G,2);
% end_time=2920;
time_end=time_scale_G(end_time);
time_stamp_end=end_time

%time0=round(0.75*end_time);
time0=1;








%%

close all

figure(1);
set(gca,'FontSize',16);
hold on;

plot(Xpos_outputG1(P1,time0:time_stamp_end),Zpos_outputG1(P1,time0:time_stamp_end),'b--','LineWidth',3);
plot(Xpos_outputG2(P2,time0:time_stamp_end),Zpos_outputG2(P2,time0:time_stamp_end),'g','LineWidth',3);
legend(EP1,EP2);

contour(scale_X,scale_Z,(psi_XZsmall_map'),12,'y','linewidth',2);
contour(scale_X,scale_Z,(psi_XZsmall_map'),[0 0],'k','linewidth',3);
xlim([-0.2 0.5]);
ylim([-0.6 0.6]);

contour(scale_X,scale_Z,Btot_XZ_map',(2:0.2:4.4));axis xy square
colormap('summer')
colorbar;

plot(Xpos_outputG1(P1,time0:time_stamp_end),Zpos_outputG1(P1,time0:time_stamp_end),'b--','LineWidth',3);
plot(Xpos_outputG2(P2,time0:time_stamp_end),Zpos_outputG2(P2,time0:time_stamp_end),'g','LineWidth',3);

xlabel('X (m)');
ylabel('Z (m)');







%%

figure(5);

subplot(3,1,1);
set(gca,'FontSize',18);
grid on;
hold on;

plot(time_scale_G(1:end_time),(mHe/eV)*(Xpos_outputG1(P1,1:end_time)+R0).*(vphi_outputG1(P1,1:end_time))-ZHe*(psi_value_outputG1(P1,1:end_time)),'b--','LineWidth',3);
plot(time_scale_G(1:end_time),(mHe/eV)*(Xpos_outputG2(P2,1:end_time)+R0).*(vphi_outputG2(P2,1:end_time))-ZHe*(psi_value_outputG2(P2,1:end_time)),'g','LineWidth',5);

set(gca,'FontSize',16);
legend(lP1,lP2);
set(gca,'FontSize',18);

xlim([time_scale_G(1) time_scale_G(end-1)])
plot(time_scale_G(1:end_time),(mHe/eV)*(Xpos_outputG1(P1,1:end_time)+R0).*(vphi_outputG1(P1,1:end_time))-ZHe*(psi_value_outputG1(P1,1:end_time)),'b--','LineWidth',4);
plot(time_scale_G(1:end_time),pphi_outputG1(P1,1:end_time),'r--','LineWidth',2);
plot(time_scale_G(1:end_time),(mHe/eV)*(Xpos_outputG2(P2,1:end_time)+R0).*(vphi_outputG2(P2,1:end_time))-ZHe*(psi_value_outputG2(P2,1:end_time)),'g','LineWidth',5);
plot(time_scale_G(1:end_time),pphi_outputG2(P2,1:end_time),'r--','LineWidth',2);

xlabel('time (s)');
yl=ylabel('p$$_\varphi$$','interpreter','latex')
set(yl,'Interpreter','latex');
xlim([0 time_end]);
% ylim([1.5 2.7])

omega_precess_avg1(P1,1:10)=0;
omega_precess_avg1(P1,end-5:end)=0;
subplot(3,1,3);
set(gca,'FontSize',18);
grid on;
hold on
% plot(time_scale_G(1:end_time),omega_precess_avg(P1,1:end_time)/mean(omega_psi_avg(P1,1:end_time)),'b--','LineWidth',4);
% plot(time_scale_G(1:end_time),omega_precess_avg(P2,1:end_time)/mean(omega_psi_avg(P2,1:end_time)),'g','LineWidth',5);
% plot(time_scale_G(1:end_time),omega_precess_avg(P1,1:end_time)/omega_cr,'b--','LineWidth',4);
% plot(time_scale_G(1:end_time),omega_precess_avg(P2,1:end_time)/omega_cr,'g','LineWidth',5);
plot(time_scale_G(1:end_time),time_scale_G(1:end_time)*0+omega_cr,'k-.','LineWidth',3);
plot(time_scale_G(1:end_time),abs(slidingavg(omega_precess_avg1(P1,1:end_time),50)),'b--','LineWidth',2);
plot(time_scale_G(1:end_time),abs(slidingavg(omega_precess_avg2(P2,1:end_time),50)),'g','LineWidth',4);
plot(time_scale_G(1:end_time),time_scale_G(1:end_time)*0+omega_cr,'k-.','LineWidth',3);
xlabel('time (s)');
ylabel('|\omega_v_D| (rad/s)');
xlim([0 time_end]);
ylim([0 2]*1e5)
legend('\omega_{cr}')


subplot(3,1,2);
set(gca,'FontSize',18);
grid on;
hold on
plot(time_scale_G(1:end_time),psi_star_outputG1(P1,1:end_time),'b--','LineWidth',2);
plot(time_scale_G(1:end_time),psi_star_outputG2(P2,1:end_time),'g','LineWidth',2);
xlabel('time (s)');
ylabel('\psi_*');
xlim([0 time_end]);

% subplot(3,1,2);
% set(gca,'FontSize',18);
% grid on;
% hold on
% plot(time_scale_G(1:end_time),ZHe*Epot_outputG1(P1,1:end_time)+0.5*(mHe/eV)*(vparallel_outputG1(P1,1:end_time).^2)+Eperp_outputG1(P1,1:end_time),'b--','LineWidth',4);
% plot(time_scale_G(1:end_time),Etot_outputG1(P1,1:end_time),'r--','LineWidth',2);
% plot(time_scale_G(1:end_time),ZHe*Epot_outputG2(P2,1:end_time)+0.5*(mHe/eV)*(vparallel_outputG2(P2,1:end_time).^2)+Eperp_outputG2(P2,1:end_time),'g','LineWidth',5);
% plot(time_scale_G(1:end_time),Etot_outputG2(P2,1:end_time),'r--','LineWidth',2);
% xlabel('time (s)');
% ylabel('Etot (eV)');
% xlim([0 time_end]);
% ylim([3.9 4.2]*1e5)



% figure(3)
% 
% subplot(2,1,1);
% set(gca,'FontSize',26);
% grid on;
% hold on
% plot(time_scale_G(1:end_time),(vparallel_outputG(P1,1:end_time)),'b--','LineWidth',4);
% plot(time_scale_G(1:end_time),(vparallel_outputG(P2,1:end_time)),'g','LineWidth',5);
% xlabel('time (s)');
% ylabel('\omega_\psi (rad/s)');
% xlim([0 time_end]);
% 
% 
% subplot(2,1,2);
% set(gca,'FontSize',26);
% grid on;
% hold on
% plot(time_scale_G(1:end_time),abs(slidingavg(omega_precess_avg(P1,1:end_time),30)),'b--','LineWidth',4);
% plot(time_scale_G(1:end_time),abs(slidingavg(omega_precess_avg(P2,1:end_time),20)),'g','LineWidth',5);
% xlabel('time (s)');
% ylabel('\omega_\psi (rad/s)');
% xlim([0 time_end]);
% 
% 
% 
% 

