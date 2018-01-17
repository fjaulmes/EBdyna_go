% close all
reset_data_analysis_environment;
load('../data_tokamak/q_profile.mat', 'psi_rank_q1')
psi_mix=size_r-4
delta_psi_q1=20
psi_core=round(psi_rank_q1-delta_psi_q1)
psi_outer=round(1.1*psi_mix)
r_mix=interp1(1:257,radial_r_value_flux,psi_mix);
r_q1=interp1(1:257,radial_r_value_flux,psi_rank_q1);
r_core=interp1(1:257,radial_r_value_flux,psi_core);
r_outer=interp1(1:257,radial_r_value_flux,psi_outer);


%%
load('initialG_alphas_vA_all_pre_collapse.mat');
load('initialG_alphas_vA_all_precession_stats.mat');
X0=alphas_pos_x;
Z0=alphas_pos_z;

alphas_theta_ini=interp2(scale_X,scale_Z,theta_XZsmall_map',alphas_pos_x,alphas_pos_z,'*linear');
psi_ini=interp2(scale_X,scale_Z,psi_norm_XZsmall_map',alphas_pos_x,alphas_pos_z,'*linear');
r_ini=interp1(1:257,radial_r_value_flux,psi_ini);

% ALL_TRAPPED_POP=(alphas_kappa>1);
% CO_PASSING_POP=(alphas_kappa<=1).*(alphas_vpll>=0);
% COUNTER_PASSING_POP=(alphas_kappa<=1).*(alphas_vpll<0);

figure(3);
grid on;
hold on;
set(gca,'FontSize',26);

alphas_pphi_ini=alphas_pphi0;

load('alphas_vA_all_collapse_Glisa_fc1h2_G110414.mat')
alphas_theta_end=interp2(scale_X,scale_Z,theta_XZsmall_map',alphas_pos_x,alphas_pos_z,'*linear');
psi_end=interp2(scale_X,scale_Z,psi_norm_XZsmall_map',alphas_pos_x,alphas_pos_z,'*linear');
r_end=interp1(1:257,radial_r_value_flux,psi_end);
delta_r=r_end-r_ini;
% PARASITE=((delta_r<0).*(r_ini<0.7*r_q1));
% find(PARASITE)
ALPHAS_POP=(~alphas_ejected);

PART_POP=find(COUNTER_PASSING_POP.*ALPHAS_POP);
plot(r_ini(PART_POP(1)),r_end(PART_POP(1)),'b.');
PART_POP=find(CO_PASSING_POP.*ALPHAS_POP);
plot(r_ini(PART_POP(1)),r_end(PART_POP(1)),'r.');
PART_POP=find(ALL_TRAPPED_POP.*ALPHAS_POP);
plot(r_ini(PART_POP(1)),r_end(PART_POP(1)),'+','color',[0 0.9 0])

plot([r_core r_core],[0 2],'r--','LineWidth',5)
% plot([0 0.8],[r_outer r_outer],'--','Color',[1.0 0.8 0.4],'LineWidth',2)
plot([r_outer r_outer],[0 2],'--','Color',[1.0 0.8 0.4],'LineWidth',5)





PART_POP=find(COUNTER_PASSING_POP.*ALPHAS_POP);
plot(r_ini(PART_POP(1:2:end)),r_end(PART_POP(1:2:end)),'b.');
PART_POP=find(CO_PASSING_POP.*ALPHAS_POP);
plot(r_ini(PART_POP(1:2:end)),r_end(PART_POP(1:2:end)),'r.');
PART_POP=find(ALL_TRAPPED_POP.*ALPHAS_POP);
plot(r_ini(PART_POP(1:2:end)),r_end(PART_POP(1:2:end)),'+','color',[0 0.9 0])
hl=legend('counter passing','co passing','trapped')
set(hl,'FontSize',22);


plot([0 2],[r_q1 r_q1],'k','LineWidth',2)
plot([r_q1 r_q1],[0 2],'k','LineWidth',2)
plot([0 2],[r_mix r_mix],'b','LineWidth',2)
plot([r_mix r_mix],[0 2],'b','LineWidth',2)






xlabel('r_{ini} (m)')
ylabel('r_{final} (m)')
plot([0 2],[0 2],'k--','LineWidth',3)
axis equal
xlim([0 1.5])
ylim([0 1.5])


% annotation('textbox',...
%     [0.5 0.8 0.0483 0.02],...
%     'EdgeColor','none',...
%     'String',{'q=1'},'fontsize',16,...
%     'FitBoxToText','off','color',[0 0 0]);
% 
% annotation('textbox',...
%     [0.8 0.392 0.0483 0.02],...
%     'EdgeColor','none',...
%     'String',{'r_{mix}'},'fontsize',16,...
%     'FitBoxToText','off','color',[0 0.0 0.9]);


annotation('textbox',...
    [0.46 0.85 0.0483 0.02],...
    'EdgeColor','none',...
    'String',{'q=1'},'fontsize',16,...
    'FitBoxToText','off','color',[0 0 0]);

annotation('textbox',...
    [0.6 0.85 0.0483 0.02],...
    'EdgeColor','none',...
    'String',{'r_{mix}'},'fontsize',16,...
    'FitBoxToText','off','color',[0 0.0 0.9]);

% annotation('textbox',...
%     [0.1 0.9 0.0483 0.02],...
%     'EdgeColor','none',...
%     'String',{'(1)'},'fontsize',22,...
%     'FitBoxToText','off','color',[1 0.0 0.0]);


%%
figure(7)

subplot(2,1,1);
set(gca,'FontSize',22);

binXlims=(-1.4:0.15:1.9);
binZlims=(-1.8:0.15:1.8);

mHist = hist2d ([Z0, X0], binZlims,binXlims);
Hist_norm=max(max(mHist));
% plot2Dhist(mHist, -0.7:0.05:0.7, -0.7:0.05:0.7,  'X', 'Z', 'Nb part ini'); 
contourf(binXlims(1:end-1)+0.025,binZlims(1:end-1)+0.025,mHist(:,1:end)/Hist_norm,(0:0.06:1))
axis xy
xlabel('X (m)')
ylabel('Z (m)')
plot_q1_psi_surface
colorbar

subplot(2,1,2);
set(gca,'FontSize',22);

mHistpost = hist2d ([alphas_pos_z, alphas_pos_x], binZlims,binXlims);
mHistpost(1,1)=Hist_norm;
% plot2Dhist(mHist, -0.7:0.05:0.7, -0.7:0.05:0.7,  'X', 'Z', 'Nb part'); 
contourf(binXlims(1:end-1)+0.5*0.15,binZlims(1:end-1)+0.5*0.15,mHistpost(:,1:end)/Hist_norm,(0:0.06:1))
axis xy
plot_q1_psi_surface
colorbar

xlabel('X (m)')
ylabel('Z (m)')

