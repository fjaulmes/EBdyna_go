% close all
load('../../data_tokamak/q_profile.mat', 'psi_rank_q1')
psi_mix=size_r-4
delta_psi_q1=20
psi_core=round(psi_rank_q1-delta_psi_q1)
psi_outer=round(1.1*psi_mix)
r_mix=interp1(1:257,radial_r_value_flux,psi_mix);
r_q1=interp1(1:257,radial_r_value_flux,psi_rank_q1);
r_core=interp1(1:257,radial_r_value_flux,psi_core);
r_outer=interp1(1:257,radial_r_value_flux,psi_outer);

load('initialG_flatD_5keV_pre_collapse.mat');
load('initialG_flatD_5keV_precession_stats.mat');
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

load('flatD_5keV_collapse_Glisa_fc2h2_G290114.mat')
psi_end=interp2(scale_X,scale_Z,psi_norm_XZsmall_map',alphas_pos_x,alphas_pos_z,'*linear');
r_end=interp1(1:257,radial_r_value_flux,psi_end);
delta_r=r_end-r_ini;
PARASITE=((delta_r<0).*(r_ini<0.7*r_q1));
find(PARASITE)
ALPHAS_POP=(~alphas_ejected).*(~PARASITE);

PART_POP=find(COUNTER_PASSING_POP.*ALPHAS_POP);
plot(r_ini(PART_POP(1)),r_end(PART_POP(1)),'b.');
PART_POP=find(CO_PASSING_POP.*ALPHAS_POP);
plot(r_ini(PART_POP(1)),r_end(PART_POP(1)),'r.');
PART_POP=find(ALL_TRAPPED_POP.*ALPHAS_POP);
plot(r_ini(PART_POP(1)),r_end(PART_POP(1)),'+','color',[0 0.9 0])

plot([0 0.8],[r_q1 r_q1],'k','LineWidth',2)
plot([r_q1 r_q1],[0 0.8],'k','LineWidth',2)
plot([0 0.8],[r_mix r_mix],'b','LineWidth',2)
plot([r_mix r_mix],[0 0.8],'b','LineWidth',2)




PART_POP=find(COUNTER_PASSING_POP.*ALPHAS_POP);
plot(r_ini(PART_POP(1:2:end)),r_end(PART_POP(1:2:end)),'b.');
PART_POP=find(CO_PASSING_POP.*ALPHAS_POP);
plot(r_ini(PART_POP(1:2:end)),r_end(PART_POP(1:2:end)),'r.');
PART_POP=find(ALL_TRAPPED_POP.*ALPHAS_POP);
plot(r_ini(PART_POP(1:2:end)),r_end(PART_POP(1:2:end)),'+','color',[0 0.9 0])
hl=legend('counter passing','co passing','trapped')
set(hl,'FontSize',22);






plot([r_core r_core],[0 0.8],'r--','LineWidth',5)
% plot([0 0.8],[r_outer r_outer],'--','Color',[1.0 0.8 0.4],'LineWidth',2)
plot([r_outer r_outer],[0 0.8],'--','Color',[1.0 0.8 0.4],'LineWidth',5)


xlabel('r_{ini} (m)')
ylabel('r_{final} (m)')
plot([0 0.68],[0 0.68],'k--','LineWidth',3)
axis equal
xlim([0 0.48])
ylim([0 0.48])


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
