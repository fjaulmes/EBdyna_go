run('reset_data_analysis_environment')
close all

load('initialG_alphas_vA_all_precession_stats.mat');

%%
load('initialG_alphas_vA_all_pre_collapse.mat');



EKIN_INF=20*1e3

PART_POP=find(alphas_Ekin>EKIN_INF);
% PART_POP=find((alphas_Ekin>EKIN_INF));

PICH_BIN_SIZE=1*1e6;
PITCH_BINS=(1:PICH_BIN_SIZE:1.5*1e7);
pitch_values=PITCH_BINS(1:end-1)+0.5*PICH_BIN_SIZE;

R_BIN_SIZE=1*1e6;
R_BINS=(-1.5*1e7:R_BIN_SIZE:1.5*1e7);
r_values=R_BINS(1:end-1)+0.5*R_BIN_SIZE;

alphas_vpll_ini=alphas_vpll;

alphas_B_ini=interp2(scale_X,scale_Z,Btot_XZ_map',alphas_pos_x,alphas_pos_z,'*linear');
alphas_Eperp=max(alphas_Ekin-0.5*(mHe/eV)*alphas_vpll.^2,0);
alphas_vperp=sqrt(2*(eV/mHe)*alphas_Eperp);
alphas_pitch=alphas_vpll./alphas_vperp;
alphas_lambda_ini=atan(1./alphas_pitch);
alphas_lambda0_ini=Bavg*alphas_mm./alphas_Ekin;
alphas_r=interp1(1:257,radial_r_value_flux,alphas_psi);
alphas_Ekin_ini=alphas_Ekin;

% here can try alphas_r or r_avg or even pphi [if small range of Ekin]

mHist = hist2d ([alphas_vpll(PART_POP) alphas_vperp(PART_POP)], R_BINS, PITCH_BINS);
% pcolor (r_values, pitch_values, mHist);
% shading faceted; 


figure(3)
subplot(2,1,1);hold on

set(gca,'fontsize',22)
contourf ( pitch_values, r_values,mHist/max(max(mHist)),(0 :0.05:1));
plot([0 4*1e6],[X_axis X_axis],'k--','Linewidth',4)
%plot2Dhist(mHist,PITCH_BINS,R_BINS,   'pitch',  'r (m)','Nb part ini'); 
axis  ij
xlabel('pitch')
ylabel('X (m)')
colorbar

load('alphas_vA_all_collapse_Glisa_fc0p8h2_G160414.mat')
PART_POP=find(alphas_Ekin>EKIN_INF);
% PART_POP=find((alphas_Ekin>EKIN_INF));


EKIN_INF=1*1e3

% PART_POP=find(ALL_PASSING_POP.*(alphas_Ekin>EKIN_INF));

% PICH_BIN_SIZE=0.08;
% PITCH_BINS=(-0.8:PICH_BIN_SIZE:0.8);
% pitch_values=PITCH_BINS(1:end-1)+0.5*PICH_BIN_SIZE;
% 
% R_BIN_SIZE=0.06;
% R_BINS=(0:R_BIN_SIZE:0.6);
% r_values=R_BINS(1:end-1)+0.5*R_BIN_SIZE;

alphas_B_end=interp2(scale_X,scale_Z,Btot_XZ_map',pos_X_gc,pos_Z_gc,'*linear');
alphas_theta_end=interp2(scale_X,scale_Z,theta_XZsmall_map',pos_X_gc,pos_Z_gc,'*linear');

alphas_Eperp=max(alphas_Ekin-0.5*(mHe/eV)*alphas_vpll.^2,0);
alphas_vperp=sqrt(2*(eV/mHe)*alphas_Eperp);
alphas_pitch=alphas_vpll./alphas_vperp;
alphas_r=interp1(1:257,radial_r_value_flux,alphas_psi);
alphas_lambda_end=atan(1./alphas_pitch);
alphas_lambda0_end=Bavg*alphas_mm./alphas_Ekin;
alphas_Ekin_end=alphas_Ekin;

% here can try alphas_r or r_avg or even pphi [if small range of Ekin]

mHistpost = hist2d ([alphas_vpll(PART_POP) alphas_vperp(PART_POP)], R_BINS, PITCH_BINS);
% pcolor (r_values, pitch_values, mHist);
% shading faceted; 


figure(3)
subplot(2,1,2);
hold on
set(gca,'fontsize',22)
% mHistpost(1,1)=max(max(mHist));
contourf ( pitch_values, r_values,mHistpost/max(max(mHistpost)),(0 :0.05:1));
plot([0 4*1e6],[X_axis X_axis],'k--','Linewidth',4)
%plot2Dhist(mHist,PITCH_BINS,R_BINS,   'pitch',  'r (m)','Nb part ini'); 
axis  ij
xlabel('Ekin (eV)')
ylabel('X (m)')
colorbar


figure(6)
hold on
mHistpost(1,1)=0;
set(gca,'fontsize',22)
mHistpost(1,1)=0;
contourf ( pitch_values, r_values,mHistpost-mHist,(-400:20:400));
plot([0 4*1e6],[X_axis X_axis],'k--','Linewidth',4)
%plot2Dhist(mHist,PITCH_BINS,R_BINS,   'pitch',  'r (m)','Nb part ini'); 
axis  ij
xlabel('Ekin (eV)')
ylabel('X (m)')


%%

figure(7)
hold on
set(gca,'fontsize',22)
contourf ( pitch_values, r_values,(mHistpost-mHist)./mHist,(-2:0.1:2));
plot([0 4*1e6],[X_axis X_axis],'k--','Linewidth',4)
%plot2Dhist(mHist,PITCH_BINS,R_BINS,   'pitch',  'r (m)','Nb part ini'); 
axis  ij
xlabel('Ekin (eV)')
ylabel('X (m)')