run('reset_data_analysis_environment')
close all

load('initial_NBI60keV_precession_stats_all.mat');

%%
load('initial_NBI60keV_pre_collapse_all.mat');



EKIN_INF=1*1e3

% PART_POP=find((alphas_Ekin>EKIN_INF).*(alphas_pos_z<0.2).*(alphas_pos_z>-0.2));
PART_POP=find((alphas_Ekin>EKIN_INF));

PICH_BIN_SIZE=0.08;
PITCH_BINS=(-1.15:PICH_BIN_SIZE:1.15);
pitch_values=PITCH_BINS(1:end-1)+0.5*PICH_BIN_SIZE;

R_BIN_SIZE=0.15;
R_BINS=(0:R_BIN_SIZE:0.6);
r_values=R_BINS(1:end-1)+0.5*R_BIN_SIZE;
volume_r_values=2*pi*R0*(pi*r_values.^2);
volume_r_values(2:end)=volume_r_values(2:end)-volume_r_values(1:end-1);

alphas_vpll_ini=alphas_vpll;

alphas_B_ini=interp2(scale_X,scale_Z,Btot_XZ_map',alphas_pos_x,alphas_pos_z,'*linear');
alphas_Eperp=max(alphas_Ekin-0.5*(mHe/eV)*alphas_vpll.^2,0);
alphas_vperp=sqrt(2*(eV/mHe)*alphas_Eperp);
alphas_vtot=sqrt(2*(eV/mHe)*alphas_Ekin);
alphas_pitch=alphas_vpll./alphas_vtot;
alphas_lambda_ini=atan(alphas_vperp./alphas_vpll);
alphas_lambda0_ini=Bavg*alphas_mm./alphas_Ekin;
alphas_r=interp1(1:257,radial_r_value_flux,alphas_psi);
alphas_Ekin_ini=alphas_Ekin;

% here can try alphas_r or r_avg or even pphi [if small range of Ekin]

mHist = hist2d ([alphas_r(PART_POP) alphas_pitch(PART_POP)], R_BINS, PITCH_BINS);
for r=1:length(r_values)
    mHist_density(r,:)=mHist(r,:)/volume_r_values(r);
end
% pcolor (r_values, pitch_values, mHist);
% shading faceted; 


figure(3)
subplot(2,1,1);hold on
set(gca,'fontsize',22)
title('before sawtooth')
mHist_contour=mHist_density/max(max(mHist_density));
mHist_contour(1,1)=1.2;
contourf ( pitch_values, r_values,mHist_contour,(0 :0.05:1.0));
plot([-1 1],[r_value_q1_mean r_value_q1_mean],'k--','Linewidth',4)
%plot2Dhist(mHist,PITCH_BINS,R_BINS,   'pitch',  'r (m)','Nb part ini'); 
axis  ij
xlabel('v_{||}/v_{tot}')
ylabel('r (m)')
colorbar
xlim([-1 1])

load('NBI60keV_fc1h1p6_all.mat')
% PART_POP=find((alphas_Ekin>EKIN_INF).*(alphas_pos_z<0.2).*(alphas_pos_z>-0.2));
% PART_POP=find((alphas_Ekin>EKIN_INF));



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
alphas_vtot=sqrt(2*(eV/mHe)*alphas_Ekin);
alphas_pitch=alphas_vpll./alphas_vtot;
alphas_lambda_end=atan(alphas_vperp./alphas_vpll);
alphas_r=interp1(1:257,radial_r_value_flux,alphas_psi);
alphas_lambda0_end=Bavg*alphas_mm./alphas_Ekin;
alphas_Ekin_end=alphas_Ekin;

% here can try alphas_r or r_avg or even pphi [if small range of Ekin]

mHistpost = hist2d ([alphas_r(PART_POP) alphas_pitch(PART_POP)], R_BINS, PITCH_BINS);
for r=1:length(r_values)
    mHistpost_density(r,:)=mHistpost(r,:)/volume_r_values(r);
end% pcolor (r_values, pitch_values, mHist);
% shading faceted; 


figure(3)
subplot(2,1,2);
set(gca,'fontsize',22)
title('after sawtooth (200 \mus)')
hold on
mHistpost_density(1,1)=max(max(mHist_density));
contourf ( pitch_values, r_values,mHistpost_density/max(max(mHist_density)),(0 :0.025:1.0));
plot([-1 1],[r_value_q1_mean r_value_q1_mean],'k--','Linewidth',4)
%plot2Dhist(mHist,PITCH_BINS,R_BINS,   'pitch',  'r (m)','Nb part ini'); 
axis  ij
xlabel('v_{||}/v_{tot}')
ylabel('r (m)')
colorbar
xlim([-1 1])


figure(6)
hold on
mHistpost(1,1)=0;
set(gca,'fontsize',22)
mHistpost(1,1)=0;
contourf ( pitch_values, r_values,mHistpost-mHist,(-4000:400:4000));
plot([-1 1],[r_value_q1_mean r_value_q1_mean],'k--','Linewidth',4)
%plot2Dhist(mHist,PITCH_BINS,R_BINS,   'pitch',  'r (m)','Nb part ini'); 
axis  ij
xlabel('pitch')
ylabel('r (m)')
xlim([-1 1])


%%
% LAMBDA_BIN_SIZE=0.04;
% LAMBDA_BINS=(0:LAMBDA_BIN_SIZE:2.2);
% l_values=LAMBDA_BINS(1:end-1)+0.5*LAMBDA_BIN_SIZE;
% 
% figure(8)
% hold on
% set(gca,'fontsize',22)
% dist_lambda_ini=histc(alphas_lambda_ini,LAMBDA_BINS)
% plot(l_values,dist_lambda_ini(1:end-1))
% xlabel('\lambda')
% 
% dist_lambda_end=histc(alphas_lambda_end,LAMBDA_BINS)
% plot(l_values,dist_lambda_end(1:end-1),'r')
% xlabel('\lambda')
% 
% xlim([0.0 2.2])
% 
% figure(9)
% hold on
% set(gca,'fontsize',22)
% dist_lambda_ini=histc(alphas_lambda0_ini,LAMBDA_BINS)
% plot(l_values,dist_lambda_ini(1:end-1))
% xlabel('\lambda')
% 
% dist_lambda_end=histc(alphas_lambda0_end,LAMBDA_BINS)
% plot(l_values,dist_lambda_end(1:end-1),'r')
% xlabel('\lambda')
% 
% xlim([0.0 1.2])
