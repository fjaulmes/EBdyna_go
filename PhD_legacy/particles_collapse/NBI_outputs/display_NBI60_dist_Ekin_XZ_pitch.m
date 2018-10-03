reset_data_analysis_environment
close all

PRE_COLLAPSE_FILENAME='initial_NBI60keV_Rlab_pre_collapse_all.mat'
STATS_FILENAME='initial_NBI60keV_R_precession_stats_all.mat'


load(PRE_COLLAPSE_FILENAME)
load(STATS_FILENAME)
alphas_psi=interp2(scale_X,scale_Z,psi_norm_XZsmall_map',pos_X_gc,pos_Z_gc,'*linear');
% r_avg=interp1(1:Nradial,radial_r_value_flux,alphas_psi);

psi_pos_q1=interp1(psi_scale,1:Nradial,psi_q1)

[eps_min inversion_radius_pos ]=min(abs(psi_star_initial(1:round(psi_pos_q1))-psi_star_final(1:round(psi_pos_q1))));
inversion_radius_pos

psi_core=1.2*inversion_radius_pos
psi_core_max=1.25*psi_core

r_core=interp1(1:Nradial,radial_r_value_flux,psi_core)

r_core_max=interp1(1:Nradial,radial_r_value_flux,psi_core_max)

EKIN_INF=1*1e3

% PART_POP=find(((alphas_pos_x-X_axis)>0.15).*((alphas_pos_x-X_axis)<0.6).*(alphas_pos_z<0.25).*(alphas_pos_z>-0.25));
% PART_POP=find((alphas_pos_z<2).*(alphas_pos_z>-2).*(alphas_psi<psi_pos_q1));
% PART_POP=find((r_avg<r_core));
% PART_POP=find((alphas_psi>r_core).*(alphas_psi<=r_core_max));
PART_POP=find((r_avg>r_core).*(r_avg<=r_core_max).*(pos_X_gc>0.1));

PICH_BIN_SIZE=0.08;
PITCH_BINS=(-1.16:PICH_BIN_SIZE:1.16);
pitch_values=PITCH_BINS(1:end-1)+0.5*PICH_BIN_SIZE;

E_BIN_SIZE=6*1e3;
E_BINS=(-E_BIN_SIZE:E_BIN_SIZE:112*1e3);
E_values=E_BINS(1:end-1)+0.5*E_BIN_SIZE;

alphas_vpll_ini=alphas_vpll;

alphas_B_ini=interp2(scale_X,scale_Z,Btot_XZ_map',pos_X_gc,pos_Z_gc,'*linear');
alphas_Eperp=max(alphas_Ekin-0.5*(mHe/eV)*alphas_vpll.^2,0);
alphas_vperp=sqrt(2*(eV/mHe)*alphas_Eperp);
alphas_vtot=sqrt(2*(eV/mHe)*alphas_Ekin);
% alphas_pitch=alphas_vpll./alphas_vtot;
% changed for MST presentation
alphas_pitch=-alphas_vpll./alphas_vtot;
alphas_lambda_ini=atan(1./alphas_pitch);
% alphas_lambda0_ini=Bavg*alphas_mm./alphas_Ekin;
alphas_r=interp1(1:257,radial_r_value_flux,alphas_psi);
alphas_Ekin_ini=alphas_Ekin;

% here can try alphas_r or r_avg or even pphi [if small range of Ekin]

mHist = hist2d ([alphas_Ekin(PART_POP) alphas_pitch(PART_POP)], E_BINS, PITCH_BINS);
% pcolor (E_values, pitch_values, mHist);
% shading faceted; 

%%
figure(3)
title('initial NBI 60 keV (inside inv. radius)')

% subplot(2,1,1);
hold on;
set(gca,'fontsize',22)
contourf ( E_values,pitch_values, mHist'/max(max(mHist)),(0 :0.05:1));
% plot([-1 1],[X_axis X_axis],'k--','Linewidth',4)
colorbar
ylim([-1 1])
xlim([0 80]*1e3)

%plot2Dhist(mHist,PITCH_BINS,E_BINS,   'pitch',  'r (m)','Nb part ini'); 
% axis  ij
% xlabel('pitch')
% ylabel('Ekin')
ylabel('pitch')
xlabel('Ekin')


load('NBI60keV_Rlab_fc1p6h1p6_all.mat');

alphas_psi_pos=interp2(scale_X,scale_Z,psi_norm_XZsmall_map',pos_X_gc,pos_Z_gc);
% PART_POP=find(((alphas_pos_x-X_axis)>0.15).*((alphas_pos_x-X_axis)<0.6).*(alphas_pos_z<0.2).*(alphas_pos_z>-0.2));
% PART_POP=find((alphas_pos_z<2).*(alphas_pos_z>-2).*(alphas_psi<psi_pos_q1));
PART_POP=find((alphas_psi_pos>psi_core).*(alphas_psi_pos<=psi_core_max).*(pos_X_gc>0.1));


% EKIN_INF=1*1e3

% PART_POP=find(ALL_PASSING_POP.*(alphas_Ekin>EKIN_INF));

% PICH_BIN_SIZE=0.08;
% PITCH_BINS=(-0.8:PICH_BIN_SIZE:0.8);
% pitch_values=PITCH_BINS(1:end-1)+0.5*PICH_BIN_SIZE;
% 
% E_BIN_SIZE=0.06;
% E_BINS=(0:E_BIN_SIZE:0.6);
% E_values=E_BINS(1:end-1)+0.5*E_BIN_SIZE;

alphas_B_end=interp2(scale_X,scale_Z,Btot_XZ_map',pos_X_gc,pos_Z_gc,'*linear');
alphas_theta_end=interp2(scale_X,scale_Z,theta_XZsmall_map',pos_X_gc,pos_Z_gc,'*linear');

alphas_Eperp=max(alphas_Ekin-0.5*(mHe/eV)*alphas_vpll.^2,0);
alphas_vperp=sqrt(2*(eV/mHe)*alphas_Eperp);
alphas_vtot=sqrt(2*(eV/mHe)*alphas_Ekin);
% alphas_pitch=alphas_vpll./alphas_vtot;
% changed for MST presentation
alphas_pitch=-alphas_vpll./alphas_vtot;
alphas_r=interp1(1:257,radial_r_value_flux,alphas_psi);
alphas_lambda_end=atan(1./alphas_pitch);
alphas_lambda0_end=Bavg*alphas_mm./alphas_Ekin;
alphas_Ekin_end=alphas_Ekin;

% here can try alphas_r or r_avg or even pphi [if small range of Ekin]

mHistpost = hist2d ([alphas_Ekin(PART_POP) alphas_pitch(PART_POP)], E_BINS, PITCH_BINS);
% pcolor (E_values, pitch_values, mHist);
% shading faceted; 


%%

figure(4)
% subplot(2,1,2);
hold on;
set(gca,'fontsize',22)
mHistpost(1,1)=max(max(mHist));
contourf ( E_values,pitch_values, mHistpost'/max(max(mHistpost)),(0 :0.05:1));
% plot([-1 1],[X_axis X_axis],'k--','Linewidth',4)
colorbar
r_core_max
%plot2Dhist(mHist,PITCH_BINS,E_BINS,   'pitch',  'r (m)','Nb part ini'); 
% axis  ij
% xlabel('pitch')
% ylabel('Ekin')
ylabel('pitch')
xlabel('Ekin (keV)')
ylim([-1 1])
xlim([0 80]*1e3)

title('after sawtooth')


%%
figure(6)
set(gca,'fontsize',22)
mHistpost(1,1)=0;

contourf (E_values, pitch_values, (mHistpost'-mHist')/max(max(mHist)),(-0.5 :0.05:0.5));
%plot2Dhist(mHist,PITCH_BINS,E_BINS,   'pitch',  'r (m)','Nb part ini'); 

ylabel('pitch')
xlabel('Ekin')
colorbar


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
