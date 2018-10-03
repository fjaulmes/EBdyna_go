
close all

REINIT_ALL_TOKAMAK_DATA=1;

if REINIT_ALL_TOKAMAK_DATA==1
    clear all
    initialize_folder_names;
    initialize_collapse_map_calculation_context
    rescaling_to_XZsmall_maps
    DT_INTERPOLATION_METHOD='quadratic'     % by default
    CALCULATE_VD_DATA_FILE=0;
end
CALCULATE_VD_DATA_FILE=1;

% load('../B_maps/B0311.mat');
% load('../E_maps/E0311.mat');
%
rmix=radial_r_value_flux(size_r-4)

frame_rank=32

load('vExB_avg_evol.mat')

vExB_radial_PR_map_avg=squeeze(vExB_radial_evol(frame_rank,:,:));
vExB_helical_PR_map_avg=squeeze(vExB_helical_evol(frame_rank,:,:));


%%

figure(1)
subplot(2,2,1)
set(gca,'FontSize',26);
imagesc(radial_r_value_flux(1:size_r),((0:NP-1)/(NP-1))*2*pi , vExB_radial_PR_map_avg);
xlim([0 0.45])
xlabel('r (m)')
ylabel('\omega (rad)')
colorbar('West')

subplot(2,2,3)
set(gca,'FontSize',26);
hold on
grid on
plot(radial_r_value_flux(1:size_r),mean(vExB_radial_PR_map_avg(1:NP-1,:),1),'linewidth',3)
xlim([0 0.45])
xlabel('r (m)')
ylabel('<v_E_r> (m/s)')


subplot(2,2,2)
set(gca,'FontSize',26);
hold on
grid on
dr_hel1=round(0.1*size_r);
dr_hel2=round(0.3*size_r);
dr_hel4=round(0.6*size_r);
% plot(((0:NP-1)/(NP-1))*2*pi,mean(vExB_radial_PR_map_avg(:,end-dr_hel1:end),2),'k--','linewidth',2)
plot(mean(vExB_radial_PR_map_avg(2:NP,round((size_r-0.5*dr_hel4)-0.5*dr_hel1):round((size_r-0.5*dr_hel4)+0.5*dr_hel1)),2),((1:NP-1)/(NP-1))*2*pi,'k--','linewidth',2)
% plot(((0:NP-1)/(NP-1))*2*pi,mean(vExB_radial_PR_map_avg(:,end-dr_hel2:end),2),'b','linewidth',3)
plot(mean(vExB_radial_PR_map_avg(2:NP,round((size_r-0.5*dr_hel4)-0.5*dr_hel2):round((size_r-0.5*dr_hel4)+0.5*dr_hel2)),2),((1:NP-1)/(NP-1))*2*pi,'b','linewidth',3)
% dr_hel3=round(0.5*size_r);
% plot(((0:NP-1)/(NP-1))*2*pi,mean(vExB_radial_PR_map_avg(:,end-dr_hel3:end),2),'k--','linewidth',3)
% plot(((0:NP-1)/(NP-1))*2*pi,mean(vExB_radial_PR_map_avg(:,end-dr_hel4:end),2),'r-.','linewidth',4)
plot(mean(vExB_radial_PR_map_avg(2:NP,round((size_r-0.5*dr_hel4)-0.5*dr_hel4):round((size_r-0.5*dr_hel4)+0.5*dr_hel4)),2),((1:NP-1)/(NP-1))*2*pi,'r-.','linewidth',4)
width1=radial_r_value_flux(size_r)-radial_r_value_flux(size_r-dr_hel1)
width2=radial_r_value_flux(size_r)-radial_r_value_flux(size_r-dr_hel2)
% width3=radial_r_value_flux(size_r)-radial_r_value_flux(size_r-dr_hel3)
width4=radial_r_value_flux(size_r)-radial_r_value_flux(size_r-dr_hel4)
% legend(strcat('\delta_r=',num2str(width1,3)),strcat('\delta_r=',num2str(width2,3)),strcat('\delta_r=',num2str(width3,3)),strcat('\delta_r=',num2str(width4,3)));
h_legend=legend(strcat('\delta_r=',num2str(width1,3)),strcat('\delta_r=',num2str(width2,3)),strcat('\delta_r=',num2str(width4,3)));
set(h_legend,'fontsize',18)
ylim([0 2*pi])
ylabel('\omega (rad)')
xlabel('<v_E_r> (m/s)')
xlim([-2000 7000])

subplot(2,2,4)
set(gca,'FontSize',26);
hold on
grid on
plot(radial_r_value_flux(1:size_r),slidingavg(vExB_radial_PR_map_avg(round(0.5*NP),:),dr_hel1),'k--','linewidth',2)
plot(radial_r_value_flux(1:size_r),slidingavg(vExB_radial_PR_map_avg(round(0.5*NP),:),dr_hel2),'b','linewidth',3)
plot(radial_r_value_flux(1:size_r),slidingavg(vExB_radial_PR_map_avg(round(0.5*NP),:),dr_hel4),'r-.','linewidth',4)
ylim([0 7000])
% title('<v_E_r> (m/s)')
xlim([0 0.45])
xlabel('r (m)')
h_legend=legend(strcat('\delta_r=',num2str(width1,3)),strcat('\delta_r=',num2str(width2,3)),strcat('\delta_r=',num2str(width4,3)));
set(h_legend,'fontsize',18)

%%

% subplot(3,2,4)
% set(gca,'FontSize',26);
% imagesc(((0:NP-1)/(NP-1))*2*pi, radial_r_value_flux(1:size_r), (vExB_helical_PR_map_avg)');
% ylim([0 0.44])
% xlim([0 2*pi])
% xlabel('\omega')
% ylabel('r')
% colorbar('East')
% 
% 
% subplot(3,2,6)
% set(gca,'FontSize',26);
% hold on
% grid on
% dr_hel1=round(0.1*size_r);
% plot(((0:NP-1)/(NP-1))*2*pi,mean(vExB_helical_PR_map_avg(:,end-dr_hel1:end),2),'k--','linewidth',2)
% dr_hel2=round(0.3*size_r);
% plot(((0:NP-1)/(NP-1))*2*pi,mean(vExB_helical_PR_map_avg(:,end-dr_hel2:end),2),'b','linewidth',3)
% % dr_hel3=round(0.5*size_r);
% % plot(((0:NP-1)/(NP-1))*2*pi,mean(vExB_radial_PR_map_avg(:,end-dr_hel3:end),2),'k--','linewidth',3)
% dr_hel4=round(0.6*size_r);
% plot(((0:NP-1)/(NP-1))*2*pi,mean(vExB_helical_PR_map_avg(:,end-dr_hel4:end),2),'r-.','linewidth',4)
% width1=radial_r_value_flux(size_r)-radial_r_value_flux(size_r-dr_hel1)
% width2=radial_r_value_flux(size_r)-radial_r_value_flux(size_r-dr_hel2)
% % width3=radial_r_value_flux(size_r)-radial_r_value_flux(size_r-dr_hel3)
% width4=radial_r_value_flux(size_r)-radial_r_value_flux(size_r-dr_hel4)
% % legend(strcat('\delta_r=',num2str(width1,3)),strcat('\delta_r=',num2str(width2,3)),strcat('\delta_r=',num2str(width3,3)),strcat('\delta_r=',num2str(width4,3)));
% legend(strcat('\delta_r=',num2str(width1,3)),strcat('\delta_r=',num2str(width2,3)),strcat('\delta_r=',num2str(width4,3)));
% xlim([0 2*pi])
% xlabel('\omega')
% ylabel('\omega_{v_E} (rad/s)')


figure(2)
set(gca,'FontSize',26);
imagesc(radial_r_value_flux(1:size_r),((0:NP-1)/(NP-1))*2*pi , abs(vExB_helical_PR_map_avg));
xlim([0 0.44])
xlabel('r (m)')
ylabel('\omega (rad)')
colorbar


figure(3)

subplot(2,1,1)
set(gca,'FontSize',26);
hold on
imagesc(((0:NP-1)/(NP-1))*2*pi, radial_r_value_flux(size_r:-1:1), (vExB_helical_PR_map_avg(:,size_r:-1:1))');
plot((0:NP-1)*2*pi/(NP-1),(0:NP-1)*0+radial_r_value_flux(round((size_r-0.5*dr_hel2))),'k--','linewidth',2)
axis ij tight
% ylim([0 0.44])
% xlim([0 2*pi])
xlabel('\omega (rad)')
ylabel('r (m)')
colorbar('west')
dr_hel1=round(0.1*size_r);
dr_hel2=round(0.3*size_r);
dr_hel3=round(0.6*size_r);
dr_hel4=round(0.8*size_r);

for (r_pos=1:size_r)
    vExB_helical_PR_map_sliding_avg1(:,r_pos)=slidingavg(vExB_helical_PR_map_avg(:,r_pos),dr_hel1);
    vExB_helical_PR_map_sliding_avg2(:,r_pos)=slidingavg(vExB_helical_PR_map_avg(:,r_pos),dr_hel2);
    vExB_helical_PR_map_sliding_avg3(:,r_pos)=slidingavg(vExB_helical_PR_map_avg(:,r_pos),dr_hel3);
    vExB_helical_PR_map_sliding_avg4(:,r_pos)=slidingavg(vExB_helical_PR_map_avg(:,r_pos),dr_hel4);
end

subplot(2,1,2)
set(gca,'FontSize',26);
hold on
grid on
plot(((0:NP-1)/(NP-1))*2*pi,vExB_helical_PR_map_sliding_avg1(:,round(size_r-0.5*dr_hel2)),'r-.','linewidth',2)
plot(((0:NP-1)/(NP-1))*2*pi,vExB_helical_PR_map_sliding_avg2(:,round(size_r-0.5*dr_hel2)),'b','linewidth',3)
plot(((0:NP-1)/(NP-1))*2*pi,vExB_helical_PR_map_sliding_avg3(:,round(size_r-0.5*dr_hel2)),'k--','linewidth',4)
plot(((0:NP-1)/(NP-1))*2*pi,vExB_helical_PR_map_sliding_avg4(:,round(size_r-0.5*dr_hel2)),'g-.','linewidth',5)
width1=radial_r_value_flux(size_r)-radial_r_value_flux(size_r-dr_hel1)
width2=radial_r_value_flux(size_r)-radial_r_value_flux(size_r-dr_hel2)
width3=radial_r_value_flux(size_r)-radial_r_value_flux(size_r-dr_hel3)
width4=radial_r_value_flux(size_r)-radial_r_value_flux(size_r-dr_hel4)
h_legend=legend(strcat('\delta_r=',num2str(width1,3)),strcat('\delta_r=',num2str(width2,3)),strcat('\delta_r=',num2str(width3,3)),strcat('\delta_r=',num2str(width4,3)));
set(h_legend,'fontsize',18)
xlim([0 2*pi])
xlabel('\omega (rad)')
ylabel('\omega_{v_E} (rad/s)')
% colorbar





% figure(4)
% subplot(2,1,1)
% set(gca,'FontSize',26);
% imagesc(((0:NP-1)/(NP-1))*2*pi, radial_r_value_flux(1:size_r), (vExB_radial_PR_map_avg)');
% ylim([0 0.44])
% xlim([0 2*pi])
% xlabel('\omega')
% ylabel('r')
% colorbar
% 
% 
% subplot(2,1,2)
% set(gca,'FontSize',26);
% hold on
% grid on
% dr_hel1=round(0.1*size_r);
% plot(((0:NP-1)/(NP-1))*2*pi,mean(vExB_radial_PR_map_avg(:,end-dr_hel1:end),2),'k--','linewidth',2)
% dr_hel2=round(0.3*size_r);
% plot(((0:NP-1)/(NP-1))*2*pi,mean(vExB_radial_PR_map_avg(:,end-dr_hel2:end),2),'b','linewidth',3)
% % dr_hel3=round(0.5*size_r);
% % plot(((0:NP-1)/(NP-1))*2*pi,mean(vExB_radial_PR_map_avg(:,end-dr_hel3:end),2),'k--','linewidth',3)
% dr_hel4=round(0.6*size_r);
% plot(((0:NP-1)/(NP-1))*2*pi,mean(vExB_radial_PR_map_avg(:,end-dr_hel4:end),2),'r-.','linewidth',4)
% width1=radial_r_value_flux(size_r)-radial_r_value_flux(size_r-dr_hel1)
% width2=radial_r_value_flux(size_r)-radial_r_value_flux(size_r-dr_hel2)
% % width3=radial_r_value_flux(size_r)-radial_r_value_flux(size_r-dr_hel3)
% width4=radial_r_value_flux(size_r)-radial_r_value_flux(size_r-dr_hel4)
% % legend(strcat('\delta_r=',num2str(width1,3)),strcat('\delta_r=',num2str(width2,3)),strcat('\delta_r=',num2str(width3,3)),strcat('\delta_r=',num2str(width4,3)));
% legend(strcat('\delta_r=',num2str(width1,3)),strcat('\delta_r=',num2str(width2,3)),strcat('\delta_r=',num2str(width4,3)));
% xlim([0 2*pi])
% xlabel('\omega')
% ylabel('\omega_{v_E} (rad/s)')
% colorbar

%%

figure(6)
subplot(2,2,1)
set(gca,'FontSize',26);
imagesc(radial_r_value_flux(1:size_r),((0:NP-1)/(NP-1))*2*pi , vExB_radial_PR_map_avg);
xlim([0.01 0.45])
xlabel('r (m)')
ylabel('\omega (rad)')
colorbar('West')
title('<v_E_r> (m/s)');

subplot(2,2,3)
set(gca,'FontSize',26);
hold on
grid on
plot(radial_r_value_flux(1:size_r),mean(vExB_radial_PR_map_avg(1:NP-1,:),1),'linewidth',3)
xlim([0.01 0.45])
xlabel('r (m)')
ylabel('<v_E_r> (m/s)')



subplot(2,2,2)
set(gca,'FontSize',26);
hold on
imagesc(((0:NP-1)/(NP-1))*2*pi, radial_r_value_flux(size_r:-1:1), (vExB_helical_PR_map_avg(:,size_r:-1:1))');
plot((0:NP-1)*2*pi/(NP-1),(0:NP-1)*0+radial_r_value_flux(round((size_r-0.5*dr_hel2))),'k--','linewidth',2)
axis ij tight
% ylim([0 0.44])
% xlim([0 2*pi])
xlabel('\omega (rad)')
ylabel('r (m)')
colorbar('west')
title('<\omega _{vE}> (rad/s)');

subplot(2,2,4)
set(gca,'FontSize',26);
hold on
grid on
plot(((0:NP-1)/(NP-1))*2*pi,vExB_helical_PR_map_sliding_avg1(:,round(size_r-0.5*dr_hel2)),'r-.','linewidth',2)
plot(((0:NP-1)/(NP-1))*2*pi,vExB_helical_PR_map_sliding_avg2(:,round(size_r-0.5*dr_hel2)),'b','linewidth',3)
plot(((0:NP-1)/(NP-1))*2*pi,vExB_helical_PR_map_sliding_avg3(:,round(size_r-0.5*dr_hel2)),'k--','linewidth',4)
plot(((0:NP-1)/(NP-1))*2*pi,vExB_helical_PR_map_sliding_avg4(:,round(size_r-0.5*dr_hel2)),'g-.','linewidth',5)
width1=radial_r_value_flux(size_r)-radial_r_value_flux(size_r-dr_hel1)
width2=radial_r_value_flux(size_r)-radial_r_value_flux(size_r-dr_hel2)
width3=radial_r_value_flux(size_r)-radial_r_value_flux(size_r-dr_hel3)
width4=radial_r_value_flux(size_r)-radial_r_value_flux(size_r-dr_hel4)
h_legend=legend(strcat('\delta_r=',num2str(width1,3)),strcat('\delta_r=',num2str(width2,3)),strcat('\delta_r=',num2str(width3,3)),strcat('\delta_r=',num2str(width4,3)));
set(h_legend,'fontsize',18)
xlim([0 2*pi])
xlabel('\omega (rad)')
ylabel('\omega_{v_E} (rad/s)')

%% end section
figure(8)
subplot(2,1,1)
set(gca,'FontSize',26);
imagesc(radial_r_value_flux(1:size_r),((0:NP-1)/(NP-1))*2*pi , vExB_radial_PR_map_avg);
xlim([0 0.45])
xlabel('r (m)')
ylabel('\omega (rad)')
colorbar('West')

subplot(2,1,2)
set(gca,'FontSize',26);
hold on
grid on
plot(radial_r_value_flux(1:size_r),mean(vExB_radial_PR_map_avg(1:NP-1,:),1),'linewidth',3)
xlim([0 0.45])
xlabel('r (m)')
ylabel('<v_E_r> (m/s)')

