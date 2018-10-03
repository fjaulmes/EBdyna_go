% time_stamp_end = 1231
% find(Ekin_outputG(:,1)<60*1e3)


% % 
P1=5;
P2=6;
P3=P1+8;
P4=P2+8;
P5=P1+8*2;
P6=P2+8*2;
P7=P1+8*3;
P8=P2+8*3;
P9=P1+8*4;
P10=P2+8*4;
P11=P1+16*3;
P12=P2+16*3;
P13=P1+16*6+8;
P14=P2+16*6+8;
P15=P1+16*7;
P16=P2+16*7;

PARTLIST=[P1 P2 P3 P4 P5 P6 P7 P8 P9 P10 P11 P12 P13 P14 P15 P16]

load('../data_tokamak/B_fields.mat', 'Btor_PR_map')
Bavg=mean(Btor_PR_map(:,1),1)

% vperp_sq_outputG=vDsq_outputG;
% Bfield_part_outputG=interp2(scale_X,scale_Z,Btot_XZ_map',Xpos_outputG,Zpos_outputG,'*linear');
mm_outputG=Eperp_outputG./Bfield_outputG;
lambda_outputG=Bavg*mm_outputG./Ekin_outputG;
% psi_value_outputG=interp2(scale_X,scale_Z,psi_XZsmall_map',Xpos_outputG,Zpos_outputG,'*linear');
% lap_psi_outputG=interp2(scale_X,scale_Z,lap_psi_XZsmall_map',Xpos_outputG,Zpos_outputG,'*linear');

lambda_P_PART=PARTLIST*0;
Ekin_PART=PARTLIST*0;
color_PART=zeros(length(PARTLIST),3);

for (n=1:8)
    color_PART(n,:)=[n/8 0 0];
end
for (n=1:8)
    color_PART(n+8,:)=[0 0 n/8];
end
for (n=1:16)
    lambda_P_PART(n)=lambda_outputG(PARTLIST(n),1);
    Ekin_PART(n)=Ekin_outputG(PARTLIST(n),1);
end

%%

phipos_outputG_wrap=wrap2pi(phipos_outputG);
omega_outputG=theta_outputG-phipos_outputG_wrap;
omega_outputG=wrap2pi(omega_outputG);

size_phi=257;
scale_tor=40*pi*((0:size_phi-1)/(size_phi-1));
Btot_X_vector=mean(Btot_XZ_map,2);
[Btot_Xphi_map tmp]=meshgrid(Btot_X_vector,ones(size_phi,1));
Btot_Xphi_map=Btot_Xphi_map';

end_time=size(time_scale,2);
ejection_time=size(time_scale,2);
% end_time=250;
%end_time=150000;
time_end=time_scale(end_time);
time_stamp_end=end_time

%time0=round(0.75*end_time);
time0=1;


close all

figure(1);
set(gca,'FontSize',16);
hold on;

for (n=1:16)
    plot(Xpos_outputG(PARTLIST(n),time0:time_stamp_end),Zpos_outputG(PARTLIST(n),time0:time_stamp_end),'color',color_PART(n,:),'LineWidth',2);
end

contour(scale_X,scale_Z,(psi_norm_XZsmall_map'-257)*psi_global/257,12,'y');
xlim([-0.5 0.7]);
ylim([-0.7 0.7]);

contour(scale_X,scale_Z,Btot_XZ_map',(2:0.2:4.4));axis xy square
colormap('summer')
colorbar;


xlabel('X (m)');
ylabel('Z (m)');






figure(3);

subplot(2,1,1);
set(gca,'FontSize',16);
grid on;
hold on;
xlim([time_scale(1) time_scale(end-1)])

for (n=1:16)
    plot(time_scale(time0:end_time),phipos_outputG(PARTLIST(n),time0:time_stamp_end),'color',color_PART(n,:),'LineWidth',2);
end

xlabel('time (s)');
ylabel('\phi (rad)');
xlim([0 time_end]);

subplot(2,1,2);
set(gca,'FontSize',16);
grid on;
hold on
xlim([time_scale(1) time_scale(end-1)])
for (n=1:16)
    plot(time_scale(time0:end_time),vparallel_outputG(PARTLIST(n),time0:time_stamp_end),'color',color_PART(n,:),'LineWidth',2);
end

xlabel('time (s)');
ylabel('v_{||} (m/s)');
xlim([0 time_end]);





% 


%%
figure(5);

subplot(2,1,1);
set(gca,'FontSize',16);
grid on;
hold on;

pphi_recalc_outputG=(mHe/eV)*(Xpos_outputG(:,1:end_time)+R0).*(vphi_outputG(:,1:end_time))-ZHe*(psi_value_outputG(:,1:end_time));

pphi_outputG=squeeze(pphi_recalc_outputG(:,1));
pphi_error=pphi_recalc_outputG*0;

for n=1:size(pphi_recalc_outputG,1)
    pphi_error(n,:)=(squeeze(pphi_outputG(n))-pphi_recalc_outputG(n,:))./squeeze(pphi_outputG(n));
end


for (n=1:16)
    plot([time_scale(1) time_scale(end_time)],squeeze(pphi_outputG(PARTLIST(n)))*[1 1],'k--','LineWidth',1);
    plot(time_scale(time0:end_time),pphi_recalc_outputG(PARTLIST(n),time0:time_stamp_end),'color',color_PART(n,:),'LineWidth',2);
end




xlabel('time (s)');
ylabel('p_\phi');
xlim([0 time_end]);


subplot(2,1,2);
set(gca,'FontSize',16);
grid on;
hold on
for (n=1:16)
    plot(time_scale(time0:end_time),pphi_error(PARTLIST(n),time0:time_stamp_end),'color',color_PART(n,:),'LineWidth',2);
end

xlabel('time (s)')
ylabel('\delta p_\phi / p_\phi');
% xlim([time_scale(1) time_scale(end-1)])
xlim([0 time_end]);

%%






figure(7)

subplot(2,1,1);
set(gca,'FontSize',16);
grid on;
hold on

for n=1:16
    plot(time_scale(time0:end_time),mm_outputG(PARTLIST(n),time0:time_stamp_end),'color',color_PART(n,:),'LineWidth',2);
end
xlabel('time (s)');
ylabel('\mu (eV T^{-1})');
% xlim([time_scale(1) time_scale(end-1)])
xlim([0 time_end]);


subplot(2,1,2);
set(gca,'FontSize',16);
grid on;
hold on
for n=1:16
    plot(time_scale(time0:end_time),rhoL_outputG(PARTLIST(n),time0:time_stamp_end),'color',color_PART(n,:),'LineWidth',2);
end

xlabel('time (s)')
ylabel('\rho_L (m)');
% xlim([time_scale(1) time_scale(end-1)])
xlim([0 time_end]);


%%


figure(9)

subplot(2,1,1);
set(gca,'FontSize',16);
grid on;
hold on

for n=1:16
    plot(time_scale(time0:end_time),psi_value_outputG(PARTLIST(n),time0:time_stamp_end),'color',color_PART(n,:),'LineWidth',2);
end

xlabel('time (s)');
ylabel('\psi (T.m^{-2})');
% xlim([time_scale(1) time_scale(end-1)])
xlim([0 time_end]);

subplot(2,1,2);
set(gca,'FontSize',16);
grid on;
hold on
for n=1:16
    plot(time_scale(time0:end_time),vphi_outputG(PARTLIST(n),time0:time_stamp_end),'color',color_PART(n,:),'LineWidth',2);
end

xlabel('time (s)')
ylabel('v_\phi');
% xlim([time_scale(1) time_scale(end-1)])
xlim([0 time_end]);


%%

figure(10)

subplot(2,1,1);
set(gca,'FontSize',16);
grid on;
hold on

Ekin_error=Ekin_half_outputG*0;

alphas_Ekin_recalc=alphas_Ekin*0;

for n=1:size(alphas_Ekin_recalc,1)
    alphas_Ekin_recalc(n)=mean(Ekin_half_outputG(n,1:10));
end
for n=1:size(pphi_recalc_outputG,1)
    Ekin_error(n,:)=(squeeze(alphas_Ekin_recalc(n))-Ekin_half_outputG(n,:))./squeeze(alphas_Ekin(n));
end

for n=1:16
    plot(time_scale(time0:end_time),Ekin_error(PARTLIST(n),time0:time_stamp_end),'color',color_PART(n,:),'LineWidth',2);
end

ylim([-1 1]*1e-11)

xlabel('time (s)');
ylabel('Ekin error');
% xlim([time_scale(1) time_scale(end-1)])
xlim([0 time_end]);

subplot(2,1,2);
set(gca,'FontSize',16);
grid on;
hold on
for n=1:16
    plot(time_scale(time0:end_time),Eperp_outputG(PARTLIST(n),time0:time_stamp_end),'color',color_PART(n,:),'LineWidth',2);
end

xlabel('time (s)')
ylabel('Eperp');
% xlim([time_scale(1) time_scale(end-1)])
xlim([0 time_end]);





