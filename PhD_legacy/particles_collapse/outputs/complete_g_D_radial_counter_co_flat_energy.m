run('reset_data_analysis_environment.m')

% close all;

Ekin_bins=[5 16 32 48 64 90 ]
Ekin_v_bins=sqrt(Ekin_bins);
Ecrit=60

NEbins=length(Ekin_bins)

rmix=radial_r_value_flux(size_r-4)
rmix=radial_r_value_flux(size_r-4)
psi_mix=size_r-4
delta_psi_q1=20;
psi_core=round(psi_rank_q1-delta_psi_q1)
psi_outer=max(round(1.1*psi_mix))

PRECESS_R_BIN_SIZE=0.04;
precession_r_bins=(0:PRECESS_R_BIN_SIZE:20*PRECESS_R_BIN_SIZE);
r_bins_values=(0.5*PRECESS_R_BIN_SIZE:PRECESS_R_BIN_SIZE:19.5*PRECESS_R_BIN_SIZE);

pphi_dr_value=zeros(length(precession_r_bins)-1,5);
pphi_g_value_counter=zeros(length(precession_r_bins)-1,5);
pphi_g_value_co=zeros(length(precession_r_bins)-1,5);
pphi_Dpphi_value_counter=zeros(length(precession_r_bins)-1,5);
pphi_Dpphi_value_co=zeros(length(precession_r_bins)-1,5);
pphi_Dpphi_value_pos=zeros(length(precession_r_bins)-1,5);
pphi_Dpphi_value_neg=zeros(length(precession_r_bins)-1,5);

g_Ekin_counter=zeros(NEbins,1);
g_Ekin_co=zeros(NEbins,1);

g_Ekin_counter1=zeros(NEbins,1);
g_Ekin_co1=zeros(NEbins,1);
g_Ekin_counter2=zeros(NEbins,1);
g_Ekin_co2=zeros(NEbins,1);


dr_Ekin_counter=zeros(NEbins,1);
dr_Ekin_co=zeros(NEbins,1);
omega_vD_Ekin_counter=zeros(NEbins,1);
omega_vD_Ekin_co=zeros(NEbins,1);
omega_psi_Ekin_counter=zeros(NEbins,1);
omega_psi_Ekin_co=zeros(NEbins,1);

omega_vD_Ekin_counter1=zeros(NEbins,1);
omega_vD_Ekin_co1=zeros(NEbins,1);
omega_psi_Ekin_counter1=zeros(NEbins,1);
omega_psi_Ekin_co1=zeros(NEbins,1);

omega_vD_Ekin_counter2=zeros(NEbins,1);
omega_vD_Ekin_co2=zeros(NEbins,1);
omega_psi_Ekin_counter2=zeros(NEbins,1);
omega_psi_Ekin_co2=zeros(NEbins,1);

Dpphi_Ekin_counter=zeros(NEbins,1);
Dpphi_Ekin_co=zeros(NEbins,1);
Dpphi_Ekin_counter_pos=zeros(NEbins,1);
Dpphi_Ekin_co_pos=zeros(NEbins,1);
Dpphi_Ekin_counter_neg=zeros(NEbins,1);
Dpphi_Ekin_co_neg=zeros(NEbins,1);

ORBIT_POP_name='ALL_PASSING_POP'

Ebin=1
load('initialG_flatD_5keV_pre_collapse.mat')
load('initialG_flatD_5keV_precession_stats.mat')
ORBIT_POP=evalin('base',ORBIT_POP_name);
REGION_POP1=(alphas_psi<=psi_core);
REGION_POP2=(alphas_psi>=psi_core).*(alphas_psi<=psi_outer);

pphi_ini=alphas_pphi0;
% r_ini=interp1(1:Nradial,radial_r_value_flux,alphas_psi);
r_ini=r_avg';
lambda_ini=Bavg*alphas_mm./alphas_Ekin;
for bin=1:(length(precession_r_bins)-1)
    PRECESS_LAMBDA_POP_CO=find((r_ini>precession_r_bins(bin)).*(r_ini<=precession_r_bins(bin+1)).*CO_PASSING_POP);
    PRECESS_LAMBDA_POP_COUNTER=find((r_ini>precession_r_bins(bin)).*(r_ini<=precession_r_bins(bin+1)).*COUNTER_PASSING_POP);
    pphi_g_value_counter(bin,Ebin)=(mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)))));
    pphi_g_value_co(bin,Ebin)=(mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_CO)))));
%     pphi_g_value_counter(bin,Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)));
%     pphi_g_value_co(bin,Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)));
    pphi_dr_value(bin,Ebin)=mean(delta_r_avg(PRECESS_LAMBDA_POP_COUNTER))/mean(delta_r_avg(PRECESS_LAMBDA_POP_CO));
end
PRECESS_LAMBDA_POP_COUNTER=find(COUNTER_PASSING_POP);
PRECESS_LAMBDA_POP_CO=find(CO_PASSING_POP);
PRECESS_LAMBDA_POP_COUNTER=find(REGION_POP1.*COUNTER_PASSING_POP);
PRECESS_LAMBDA_POP_CO=find(REGION_POP1.*CO_PASSING_POP);
g_Ekin_counter1(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)));
g_Ekin_co1(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_CO)));
PRECESS_LAMBDA_POP_COUNTER=find(REGION_POP2.*COUNTER_PASSING_POP);
PRECESS_LAMBDA_POP_CO=find(REGION_POP2.*CO_PASSING_POP);
g_Ekin_counter2(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)));
g_Ekin_co2(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_CO)));
dr_Ekin_counter(Ebin)=mean(delta_r_avg(PRECESS_LAMBDA_POP_COUNTER));
dr_Ekin_co(Ebin)=mean(delta_r_avg(PRECESS_LAMBDA_POP_CO));
omega_vD_Ekin_counter(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)));
omega_vD_Ekin_co(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)));
omega_psi_Ekin_counter(Ebin)=mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)));
omega_psi_Ekin_co(Ebin)=mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_CO)));

PRECESS_LAMBDA_POP_COUNTER=find(COUNTER_PASSING_POP.*REGION_POP1);
PRECESS_LAMBDA_POP_CO=find(CO_PASSING_POP.*REGION_POP1);
omega_vD_Ekin_counter1(Ebin)=mean((omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)));
omega_vD_Ekin_co1(Ebin)=mean((omega_precess_avg(PRECESS_LAMBDA_POP_CO)));
omega_psi_Ekin_counter1(Ebin)=mean((omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)));
omega_psi_Ekin_co1(Ebin)=mean((omega_psi_avg(PRECESS_LAMBDA_POP_CO)));
PRECESS_LAMBDA_POP_COUNTER=find(COUNTER_PASSING_POP.*REGION_POP2);
PRECESS_LAMBDA_POP_CO=find(CO_PASSING_POP.*REGION_POP2);
omega_vD_Ekin_counter2(Ebin)=mean((omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)));
omega_vD_Ekin_co2(Ebin)=mean((omega_precess_avg(PRECESS_LAMBDA_POP_CO)));
omega_psi_Ekin_counter2(Ebin)=mean((omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)));
omega_psi_Ekin_co2(Ebin)=mean((omega_psi_avg(PRECESS_LAMBDA_POP_CO)));

load('flatD_5keV_collapse_Glisa_fc2h2_G290114.mat')
lambda_end=Bavg*alphas_mm./alphas_Ekin;
pphi_end=alphas_pphi0;
delta_pphi=pphi_end-pphi_ini; 
delta_lambda=lambda_end-lambda_ini;
ALPHAS_POP=(~alphas_ejected).*(alphas_psi<=psi_outer);
for bin=1:(length(precession_r_bins)-1)
    PRECESS_LAMBDA_POP_CO=find((r_ini>precession_r_bins(bin)).*(r_ini<=precession_r_bins(bin+1)).*CO_PASSING_POP);
    PRECESS_LAMBDA_POP_COUNTER=find((r_ini>precession_r_bins(bin)).*(r_ini<=precession_r_bins(bin+1)).*COUNTER_PASSING_POP);
    pphi_Dpphi_value_counter(bin,Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP_COUNTER)));
    pphi_Dpphi_value_co(bin,Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP_CO)));
end
PRECESS_LAMBDA_POP_COUNTER=find(COUNTER_PASSING_POP.*ALPHAS_POP);
PRECESS_LAMBDA_POP_CO=find(CO_PASSING_POP.*ALPHAS_POP);
PRECESS_LAMBDA_POP_CO_POS=find(REGION_POP1.*CO_PASSING_POP.*ALPHAS_POP);
PRECESS_LAMBDA_POP_COUNTER_POS=find(REGION_POP1.*COUNTER_PASSING_POP.*ALPHAS_POP);
PRECESS_LAMBDA_POP_CO_NEG=find(REGION_POP2.*CO_PASSING_POP.*ALPHAS_POP);
PRECESS_LAMBDA_POP_COUNTER_NEG=find(REGION_POP2.*COUNTER_PASSING_POP.*ALPHAS_POP);
g_Ekin_counter(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)));
g_Ekin_co(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_CO)));
Dpphi_Ekin_counter(Ebin)=mean(abs(delta_pphi(PRECESS_LAMBDA_POP_COUNTER)));
Dpphi_Ekin_co(Ebin)=mean(abs(delta_pphi(PRECESS_LAMBDA_POP_CO)));
Dpphi_Ekin_counter_pos(Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP_COUNTER_POS)));
Dpphi_Ekin_co_pos(Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP_CO_POS)));
Dpphi_Ekin_counter_neg(Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP_COUNTER_NEG)));
Dpphi_Ekin_co_neg(Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP_CO_NEG)));

Ebin1value=round(mean(alphas_Ekin)*1e-3)

% Ebin=2
% load('initialG_200keV_flat_pre_collapse.mat')
% load('initialG_200keV_flat_precession_stats.mat')
% ORBIT_POP=evalin('base',ORBIT_POP_name);
% r_ini=alphas_pphi0;
% lambda_ini=Bavg*alphas_mm./alphas_Ekin;
% for bin=1:(length(precession_lambda_bins)-1)
%     PRECESS_LAMBDA_POP_CO=find((alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*CO_PASSING_POP);
%     PRECESS_LAMBDA_POP_COUNTER=find((alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*COUNTER_PASSING_POP);
%     lambda_g_value(bin,Ebin)=(mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER))/mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)))));
%     lambda_g_value(bin,Ebin)=lambda_g_value(bin,Ebin)*(mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_CO))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)))));
%     lambda_dr_value(bin,Ebin)=mean(delta_r_avg(PRECESS_LAMBDA_POP_COUNTER))/mean(delta_r_avg(PRECESS_LAMBDA_POP_CO));
% end
% PRECESS_LAMBDA_POP_COUNTER=find(COUNTER_PASSING_POP);
% PRECESS_LAMBDA_POP_CO=find(CO_PASSING_POP);
% g_Ekin_counter(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)));
% g_Ekin_co(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_CO)));
% dr_Ekin_counter(Ebin)=mean(delta_r_avg(PRECESS_LAMBDA_POP_COUNTER));
% dr_Ekin_co(Ebin)=mean(delta_r_avg(PRECESS_LAMBDA_POP_CO));
% omega_vD_Ekin_counter(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)));
% omega_vD_Ekin_co(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)));
% omega_psi_Ekin_counter(Ebin)=mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)));
% omega_psi_Ekin_co(Ebin)=mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_CO)));
% load('flat200keV_collapse_Glocal_fc1h2_G210813.mat')
% lambda_end=Bavg*alphas_mm./alphas_Ekin;
% pphi_end=alphas_pphi0;
% delta_pphi=pphi_end-pphi_ini; 
% delta_lambda=lambda_end-lambda_ini;
% ALPHAS_POP=(~alphas_ejected).*(alphas_psi<size_r+DELTA_SIZE);
% for bin=1:(length(precession_lambda_bins)-1)
%     PRECESS_LAMBDA_POP_CO=find((alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*CO_PASSING_POP.*ALPHAS_POP);
%     PRECESS_LAMBDA_POP_COUNTER=find((alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*COUNTER_PASSING_POP.*ALPHAS_POP);
%     lambda_Dpphi_value(bin,Ebin)=mean(abs(delta_pphi(PRECESS_LAMBDA_POP_COUNTER)))/mean(abs(delta_pphi(PRECESS_LAMBDA_POP_CO)));
% end
% PRECESS_LAMBDA_POP_COUNTER=find(COUNTER_PASSING_POP.*ALPHAS_POP);
% PRECESS_LAMBDA_POP_CO=find(CO_PASSING_POP.*ALPHAS_POP);
% Dpphi_Ekin_counter(Ebin)=mean(abs(delta_pphi(PRECESS_LAMBDA_POP_COUNTER)));
% Dpphi_Ekin_co(Ebin)=mean(abs(delta_pphi(PRECESS_LAMBDA_POP_CO)));




Ebin=Ebin+1
load('initialG_flatD_16keV_pre_collapse.mat')
load('initialG_flatD_16keV_precession_stats.mat')
ORBIT_POP=evalin('base',ORBIT_POP_name);
REGION_POP1=(alphas_psi<=psi_core);
REGION_POP2=(alphas_psi>=psi_core).*(alphas_psi<=psi_outer);

pphi_ini=alphas_pphi0;
% r_ini=interp1(1:Nradial,radial_r_value_flux,alphas_psi);
r_ini=r_avg';
lambda_ini=Bavg*alphas_mm./alphas_Ekin;
for bin=1:(length(precession_r_bins)-1)
    PRECESS_LAMBDA_POP_CO=find((r_ini>precession_r_bins(bin)).*(r_ini<=precession_r_bins(bin+1)).*CO_PASSING_POP);
    PRECESS_LAMBDA_POP_COUNTER=find((r_ini>precession_r_bins(bin)).*(r_ini<=precession_r_bins(bin+1)).*COUNTER_PASSING_POP);
    pphi_g_value_counter(bin,Ebin)=(mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)))));
    pphi_g_value_co(bin,Ebin)=(mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_CO)))));
%     pphi_g_value_counter(bin,Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)));
%     pphi_g_value_co(bin,Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)));
    pphi_dr_value(bin,Ebin)=mean(delta_r_avg(PRECESS_LAMBDA_POP_COUNTER))/mean(delta_r_avg(PRECESS_LAMBDA_POP_CO));
end
PRECESS_LAMBDA_POP_COUNTER=find(COUNTER_PASSING_POP);
PRECESS_LAMBDA_POP_CO=find(CO_PASSING_POP);
PRECESS_LAMBDA_POP_COUNTER=find(REGION_POP1.*COUNTER_PASSING_POP);
PRECESS_LAMBDA_POP_CO=find(REGION_POP1.*CO_PASSING_POP);
g_Ekin_counter1(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)));
g_Ekin_co1(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_CO)));
PRECESS_LAMBDA_POP_COUNTER=find(REGION_POP2.*COUNTER_PASSING_POP);
PRECESS_LAMBDA_POP_CO=find(REGION_POP2.*CO_PASSING_POP);
g_Ekin_counter2(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)));
g_Ekin_co2(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_CO)));
dr_Ekin_counter(Ebin)=mean(delta_r_avg(PRECESS_LAMBDA_POP_COUNTER));
dr_Ekin_co(Ebin)=mean(delta_r_avg(PRECESS_LAMBDA_POP_CO));
omega_vD_Ekin_counter(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)));
omega_vD_Ekin_co(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)));
omega_psi_Ekin_counter(Ebin)=mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)));
omega_psi_Ekin_co(Ebin)=mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_CO)));

PRECESS_LAMBDA_POP_COUNTER=find(COUNTER_PASSING_POP.*REGION_POP1);
PRECESS_LAMBDA_POP_CO=find(CO_PASSING_POP.*REGION_POP1);
omega_vD_Ekin_counter1(Ebin)=mean((omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)));
omega_vD_Ekin_co1(Ebin)=mean((omega_precess_avg(PRECESS_LAMBDA_POP_CO)));
omega_psi_Ekin_counter1(Ebin)=mean((omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)));
omega_psi_Ekin_co1(Ebin)=mean((omega_psi_avg(PRECESS_LAMBDA_POP_CO)));
PRECESS_LAMBDA_POP_COUNTER=find(COUNTER_PASSING_POP.*REGION_POP2);
PRECESS_LAMBDA_POP_CO=find(CO_PASSING_POP.*REGION_POP2);
omega_vD_Ekin_counter2(Ebin)=mean((omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)));
omega_vD_Ekin_co2(Ebin)=mean((omega_precess_avg(PRECESS_LAMBDA_POP_CO)));
omega_psi_Ekin_counter2(Ebin)=mean((omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)));
omega_psi_Ekin_co2(Ebin)=mean((omega_psi_avg(PRECESS_LAMBDA_POP_CO)));

load('flatD_16keV_collapse_Glisa_fc2h2_G290114.mat')
lambda_end=Bavg*alphas_mm./alphas_Ekin;
pphi_end=alphas_pphi0;
delta_pphi=pphi_end-pphi_ini; 
delta_lambda=lambda_end-lambda_ini;
ALPHAS_POP=(~alphas_ejected).*(alphas_psi<=psi_outer);
for bin=1:(length(precession_r_bins)-1)
    PRECESS_LAMBDA_POP_CO=find((r_ini>precession_r_bins(bin)).*(r_ini<=precession_r_bins(bin+1)).*CO_PASSING_POP.*ALPHAS_POP);
    PRECESS_LAMBDA_POP_COUNTER=find((r_ini>precession_r_bins(bin)).*(r_ini<=precession_r_bins(bin+1)).*COUNTER_PASSING_POP.*ALPHAS_POP);
    pphi_Dpphi_value_counter(bin,Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP_COUNTER)));
    pphi_Dpphi_value_co(bin,Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP_CO)));
end
PRECESS_LAMBDA_POP_COUNTER=find(COUNTER_PASSING_POP.*ALPHAS_POP);
PRECESS_LAMBDA_POP_CO=find(CO_PASSING_POP.*ALPHAS_POP);
PRECESS_LAMBDA_POP_CO_POS=find(REGION_POP1.*CO_PASSING_POP.*ALPHAS_POP);
PRECESS_LAMBDA_POP_COUNTER_POS=find(REGION_POP1.*COUNTER_PASSING_POP.*ALPHAS_POP);
PRECESS_LAMBDA_POP_CO_NEG=find(REGION_POP2.*CO_PASSING_POP.*ALPHAS_POP);
PRECESS_LAMBDA_POP_COUNTER_NEG=find(REGION_POP2.*COUNTER_PASSING_POP.*ALPHAS_POP);
g_Ekin_counter(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)));
g_Ekin_co(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_CO)));
Dpphi_Ekin_counter(Ebin)=mean(abs(delta_pphi(PRECESS_LAMBDA_POP_COUNTER)));
Dpphi_Ekin_co(Ebin)=mean(abs(delta_pphi(PRECESS_LAMBDA_POP_CO)));
Dpphi_Ekin_counter_pos(Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP_COUNTER_POS)));
Dpphi_Ekin_co_pos(Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP_CO_POS)));
Dpphi_Ekin_counter_neg(Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP_COUNTER_NEG)));
Dpphi_Ekin_co_neg(Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP_CO_NEG)));




Ebin=Ebin+1

load('initialG_flatD_32keV_pre_collapse.mat')
load('initialG_flatD_32keV_precession_stats.mat')
ORBIT_POP=evalin('base',ORBIT_POP_name);
REGION_POP1=(alphas_psi<=psi_core);
REGION_POP2=(alphas_psi>=psi_core).*(alphas_psi<=psi_outer);

pphi_ini=alphas_pphi0;
% r_ini=interp1(1:Nradial,radial_r_value_flux,alphas_psi);
r_ini=r_avg';
lambda_ini=Bavg*alphas_mm./alphas_Ekin;
for bin=1:(length(precession_r_bins)-1)
    PRECESS_LAMBDA_POP_CO=find((r_ini>precession_r_bins(bin)).*(r_ini<=precession_r_bins(bin+1)).*CO_PASSING_POP);
    PRECESS_LAMBDA_POP_COUNTER=find((r_ini>precession_r_bins(bin)).*(r_ini<=precession_r_bins(bin+1)).*COUNTER_PASSING_POP);
    pphi_g_value_counter(bin,Ebin)=(mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)))));
    pphi_g_value_co(bin,Ebin)=(mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_CO)))));
%     pphi_g_value_counter(bin,Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)));
%     pphi_g_value_co(bin,Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)));
    pphi_dr_value(bin,Ebin)=mean(delta_r_avg(PRECESS_LAMBDA_POP_COUNTER))/mean(delta_r_avg(PRECESS_LAMBDA_POP_CO));
end
PRECESS_LAMBDA_POP_COUNTER=find(COUNTER_PASSING_POP);
PRECESS_LAMBDA_POP_CO=find(CO_PASSING_POP);
PRECESS_LAMBDA_POP_COUNTER=find(REGION_POP1.*COUNTER_PASSING_POP);
PRECESS_LAMBDA_POP_CO=find(REGION_POP1.*CO_PASSING_POP);
g_Ekin_counter1(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)));
g_Ekin_co1(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_CO)));
PRECESS_LAMBDA_POP_COUNTER=find(REGION_POP2.*COUNTER_PASSING_POP);
PRECESS_LAMBDA_POP_CO=find(REGION_POP2.*CO_PASSING_POP);
g_Ekin_counter2(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)));
g_Ekin_co2(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_CO)));
dr_Ekin_counter(Ebin)=mean(delta_r_avg(PRECESS_LAMBDA_POP_COUNTER));
dr_Ekin_co(Ebin)=mean(delta_r_avg(PRECESS_LAMBDA_POP_CO));
omega_vD_Ekin_counter(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)));
omega_vD_Ekin_co(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)));
omega_psi_Ekin_counter(Ebin)=mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)));
omega_psi_Ekin_co(Ebin)=mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_CO)));

PRECESS_LAMBDA_POP_COUNTER=find(COUNTER_PASSING_POP.*REGION_POP1);
PRECESS_LAMBDA_POP_CO=find(CO_PASSING_POP.*REGION_POP1);
omega_vD_Ekin_counter1(Ebin)=mean((omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)));
omega_vD_Ekin_co1(Ebin)=mean((omega_precess_avg(PRECESS_LAMBDA_POP_CO)));
omega_psi_Ekin_counter1(Ebin)=mean((omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)));
omega_psi_Ekin_co1(Ebin)=mean((omega_psi_avg(PRECESS_LAMBDA_POP_CO)));
PRECESS_LAMBDA_POP_COUNTER=find(COUNTER_PASSING_POP.*REGION_POP2);
PRECESS_LAMBDA_POP_CO=find(CO_PASSING_POP.*REGION_POP2);
omega_vD_Ekin_counter2(Ebin)=mean((omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)));
omega_vD_Ekin_co2(Ebin)=mean((omega_precess_avg(PRECESS_LAMBDA_POP_CO)));
omega_psi_Ekin_counter2(Ebin)=mean((omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)));
omega_psi_Ekin_co2(Ebin)=mean((omega_psi_avg(PRECESS_LAMBDA_POP_CO)));

load('flatD_32keV_collapse_Glisa_fc2h2_G290114.mat')
lambda_end=Bavg*alphas_mm./alphas_Ekin;
pphi_end=alphas_pphi0;
delta_pphi=pphi_end-pphi_ini; 
delta_lambda=lambda_end-lambda_ini;
ALPHAS_POP=(~alphas_ejected).*(alphas_psi<=psi_outer);
for bin=1:(length(precession_r_bins)-1)
    PRECESS_LAMBDA_POP_CO=find((r_ini>precession_r_bins(bin)).*(r_ini<=precession_r_bins(bin+1)).*CO_PASSING_POP.*ALPHAS_POP);
    PRECESS_LAMBDA_POP_COUNTER=find((r_ini>precession_r_bins(bin)).*(r_ini<=precession_r_bins(bin+1)).*COUNTER_PASSING_POP.*ALPHAS_POP);
    pphi_Dpphi_value_counter(bin,Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP_COUNTER)));
    pphi_Dpphi_value_co(bin,Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP_CO)));
end
PRECESS_LAMBDA_POP_COUNTER=find(COUNTER_PASSING_POP.*ALPHAS_POP);
PRECESS_LAMBDA_POP_CO=find(CO_PASSING_POP.*ALPHAS_POP);
PRECESS_LAMBDA_POP_CO_POS=find(REGION_POP1.*CO_PASSING_POP.*ALPHAS_POP);
PRECESS_LAMBDA_POP_COUNTER_POS=find(REGION_POP1.*COUNTER_PASSING_POP.*ALPHAS_POP);
PRECESS_LAMBDA_POP_CO_NEG=find(REGION_POP2.*CO_PASSING_POP.*ALPHAS_POP);
PRECESS_LAMBDA_POP_COUNTER_NEG=find(REGION_POP2.*COUNTER_PASSING_POP.*ALPHAS_POP);
g_Ekin_counter(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)));
g_Ekin_co(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_CO)));
Dpphi_Ekin_counter(Ebin)=mean(abs(delta_pphi(PRECESS_LAMBDA_POP_COUNTER)));
Dpphi_Ekin_co(Ebin)=mean(abs(delta_pphi(PRECESS_LAMBDA_POP_CO)));
Dpphi_Ekin_counter_pos(Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP_COUNTER_POS)));
Dpphi_Ekin_co_pos(Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP_CO_POS)));
Dpphi_Ekin_counter_neg(Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP_COUNTER_NEG)));
Dpphi_Ekin_co_neg(Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP_CO_NEG)));


Ebin=Ebin+1

load('initialG_flatD_48keV_pre_collapse.mat')
load('initialG_flatD_48keV_precession_stats.mat')
ORBIT_POP=evalin('base',ORBIT_POP_name);
REGION_POP1=(alphas_psi<=psi_core);
REGION_POP2=(alphas_psi>=psi_core).*(alphas_psi<=psi_outer);

pphi_ini=alphas_pphi0;
% r_ini=interp1(1:Nradial,radial_r_value_flux,alphas_psi);
r_ini=r_avg';
lambda_ini=Bavg*alphas_mm./alphas_Ekin;
for bin=1:(length(precession_r_bins)-1)
    PRECESS_LAMBDA_POP_CO=find((r_ini>precession_r_bins(bin)).*(r_ini<=precession_r_bins(bin+1)).*CO_PASSING_POP);
    PRECESS_LAMBDA_POP_COUNTER=find((r_ini>precession_r_bins(bin)).*(r_ini<=precession_r_bins(bin+1)).*COUNTER_PASSING_POP);
    pphi_g_value_counter(bin,Ebin)=(mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)))));
    pphi_g_value_co(bin,Ebin)=(mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_CO)))));
%     pphi_g_value_counter(bin,Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)));
%     pphi_g_value_co(bin,Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)));
    pphi_dr_value(bin,Ebin)=mean(delta_r_avg(PRECESS_LAMBDA_POP_COUNTER))/mean(delta_r_avg(PRECESS_LAMBDA_POP_CO));
end
PRECESS_LAMBDA_POP_COUNTER=find(COUNTER_PASSING_POP);
PRECESS_LAMBDA_POP_CO=find(CO_PASSING_POP);
PRECESS_LAMBDA_POP_COUNTER=find(REGION_POP1.*COUNTER_PASSING_POP);
PRECESS_LAMBDA_POP_CO=find(REGION_POP1.*CO_PASSING_POP);
g_Ekin_counter1(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)));
g_Ekin_co1(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_CO)));
PRECESS_LAMBDA_POP_COUNTER=find(REGION_POP2.*COUNTER_PASSING_POP);
PRECESS_LAMBDA_POP_CO=find(REGION_POP2.*CO_PASSING_POP);
g_Ekin_counter2(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)));
g_Ekin_co2(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_CO)));
dr_Ekin_counter(Ebin)=mean(delta_r_avg(PRECESS_LAMBDA_POP_COUNTER));
dr_Ekin_co(Ebin)=mean(delta_r_avg(PRECESS_LAMBDA_POP_CO));
omega_vD_Ekin_counter(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)));
omega_vD_Ekin_co(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)));
omega_psi_Ekin_counter(Ebin)=mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)));
omega_psi_Ekin_co(Ebin)=mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_CO)));

PRECESS_LAMBDA_POP_COUNTER=find(COUNTER_PASSING_POP.*REGION_POP1);
PRECESS_LAMBDA_POP_CO=find(CO_PASSING_POP.*REGION_POP1);
omega_vD_Ekin_counter1(Ebin)=mean((omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)));
omega_vD_Ekin_co1(Ebin)=mean((omega_precess_avg(PRECESS_LAMBDA_POP_CO)));
omega_psi_Ekin_counter1(Ebin)=mean((omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)));
omega_psi_Ekin_co1(Ebin)=mean((omega_psi_avg(PRECESS_LAMBDA_POP_CO)));
PRECESS_LAMBDA_POP_COUNTER=find(COUNTER_PASSING_POP.*REGION_POP2);
PRECESS_LAMBDA_POP_CO=find(CO_PASSING_POP.*REGION_POP2);
omega_vD_Ekin_counter2(Ebin)=mean((omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)));
omega_vD_Ekin_co2(Ebin)=mean((omega_precess_avg(PRECESS_LAMBDA_POP_CO)));
omega_psi_Ekin_counter2(Ebin)=mean((omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)));
omega_psi_Ekin_co2(Ebin)=mean((omega_psi_avg(PRECESS_LAMBDA_POP_CO)));

load('flatD_48keV_collapse_Glisa_fc2h2_G290114.mat')
lambda_end=Bavg*alphas_mm./alphas_Ekin;
pphi_end=alphas_pphi0;
delta_pphi=pphi_end-pphi_ini; 
delta_lambda=lambda_end-lambda_ini;
ALPHAS_POP=(~alphas_ejected).*(alphas_psi<=psi_outer);
for bin=1:(length(precession_r_bins)-1)
    PRECESS_LAMBDA_POP_CO=find((r_ini>precession_r_bins(bin)).*(r_ini<=precession_r_bins(bin+1)).*CO_PASSING_POP.*ALPHAS_POP);
    PRECESS_LAMBDA_POP_COUNTER=find((r_ini>precession_r_bins(bin)).*(r_ini<=precession_r_bins(bin+1)).*COUNTER_PASSING_POP.*ALPHAS_POP);
    pphi_Dpphi_value_counter(bin,Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP_COUNTER)));
    pphi_Dpphi_value_co(bin,Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP_CO)));
end
PRECESS_LAMBDA_POP_COUNTER=find(COUNTER_PASSING_POP.*ALPHAS_POP);
PRECESS_LAMBDA_POP_CO=find(CO_PASSING_POP.*ALPHAS_POP);
PRECESS_LAMBDA_POP_CO_POS=find(REGION_POP1.*CO_PASSING_POP.*ALPHAS_POP);
PRECESS_LAMBDA_POP_COUNTER_POS=find(REGION_POP1.*COUNTER_PASSING_POP.*ALPHAS_POP);
PRECESS_LAMBDA_POP_CO_NEG=find(REGION_POP2.*CO_PASSING_POP.*ALPHAS_POP);
PRECESS_LAMBDA_POP_COUNTER_NEG=find(REGION_POP2.*COUNTER_PASSING_POP.*ALPHAS_POP);
g_Ekin_counter(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)));
g_Ekin_co(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_CO)));
Dpphi_Ekin_counter(Ebin)=mean(abs(delta_pphi(PRECESS_LAMBDA_POP_COUNTER)));
Dpphi_Ekin_co(Ebin)=mean(abs(delta_pphi(PRECESS_LAMBDA_POP_CO)));
Dpphi_Ekin_counter_pos(Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP_COUNTER_POS)));
Dpphi_Ekin_co_pos(Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP_CO_POS)));
Dpphi_Ekin_counter_neg(Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP_COUNTER_NEG)));
Dpphi_Ekin_co_neg(Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP_CO_NEG)));


Ebin=Ebin+1

load('initialG_flatD_64keV_pre_collapse.mat')
load('initialG_flatD_64keV_precession_stats.mat')
ORBIT_POP=evalin('base',ORBIT_POP_name);
REGION_POP1=(alphas_psi<=psi_core);
REGION_POP2=(alphas_psi>=psi_core).*(alphas_psi<=psi_outer);

pphi_ini=alphas_pphi0;
% r_ini=interp1(1:Nradial,radial_r_value_flux,alphas_psi);
r_ini=r_avg';
lambda_ini=Bavg*alphas_mm./alphas_Ekin;
for bin=1:(length(precession_r_bins)-1)
    PRECESS_LAMBDA_POP_CO=find((r_ini>precession_r_bins(bin)).*(r_ini<=precession_r_bins(bin+1)).*CO_PASSING_POP);
    PRECESS_LAMBDA_POP_COUNTER=find((r_ini>precession_r_bins(bin)).*(r_ini<=precession_r_bins(bin+1)).*COUNTER_PASSING_POP);
    pphi_g_value_counter(bin,Ebin)=(mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)))));
    pphi_g_value_co(bin,Ebin)=(mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_CO)))));
%     pphi_g_value_counter(bin,Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)));
%     pphi_g_value_co(bin,Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)));
    pphi_dr_value(bin,Ebin)=mean(delta_r_avg(PRECESS_LAMBDA_POP_COUNTER))/mean(delta_r_avg(PRECESS_LAMBDA_POP_CO));
end
PRECESS_LAMBDA_POP_COUNTER=find(COUNTER_PASSING_POP);
PRECESS_LAMBDA_POP_CO=find(CO_PASSING_POP);
PRECESS_LAMBDA_POP_COUNTER=find(REGION_POP1.*COUNTER_PASSING_POP);
PRECESS_LAMBDA_POP_CO=find(REGION_POP1.*CO_PASSING_POP);
g_Ekin_counter1(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)));
g_Ekin_co1(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_CO)));
PRECESS_LAMBDA_POP_COUNTER=find(REGION_POP2.*COUNTER_PASSING_POP);
PRECESS_LAMBDA_POP_CO=find(REGION_POP2.*CO_PASSING_POP);
g_Ekin_counter2(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)));
g_Ekin_co2(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_CO)));
dr_Ekin_counter(Ebin)=mean(delta_r_avg(PRECESS_LAMBDA_POP_COUNTER));
dr_Ekin_co(Ebin)=mean(delta_r_avg(PRECESS_LAMBDA_POP_CO));
omega_vD_Ekin_counter(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)));
omega_vD_Ekin_co(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)));
omega_psi_Ekin_counter(Ebin)=mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)));
omega_psi_Ekin_co(Ebin)=mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_CO)));

PRECESS_LAMBDA_POP_COUNTER=find(COUNTER_PASSING_POP.*REGION_POP1);
PRECESS_LAMBDA_POP_CO=find(CO_PASSING_POP.*REGION_POP1);
omega_vD_Ekin_counter1(Ebin)=mean((omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)));
omega_vD_Ekin_co1(Ebin)=mean((omega_precess_avg(PRECESS_LAMBDA_POP_CO)));
omega_psi_Ekin_counter1(Ebin)=mean((omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)));
omega_psi_Ekin_co1(Ebin)=mean((omega_psi_avg(PRECESS_LAMBDA_POP_CO)));
PRECESS_LAMBDA_POP_COUNTER=find(COUNTER_PASSING_POP.*REGION_POP2);
PRECESS_LAMBDA_POP_CO=find(CO_PASSING_POP.*REGION_POP2);
omega_vD_Ekin_counter2(Ebin)=mean((omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)));
omega_vD_Ekin_co2(Ebin)=mean((omega_precess_avg(PRECESS_LAMBDA_POP_CO)));
omega_psi_Ekin_counter2(Ebin)=mean((omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)));
omega_psi_Ekin_co2(Ebin)=mean((omega_psi_avg(PRECESS_LAMBDA_POP_CO)));

load('flatD_64keV_collapse_Glisa_fc2h2_G290114.mat')
lambda_end=Bavg*alphas_mm./alphas_Ekin;
pphi_end=alphas_pphi0;
delta_pphi=pphi_end-pphi_ini; 
delta_lambda=lambda_end-lambda_ini;
ALPHAS_POP=(~alphas_ejected).*(alphas_psi<=psi_outer);
for bin=1:(length(precession_r_bins)-1)
    PRECESS_LAMBDA_POP_CO=find((r_ini>precession_r_bins(bin)).*(r_ini<=precession_r_bins(bin+1)).*CO_PASSING_POP.*ALPHAS_POP);
    PRECESS_LAMBDA_POP_COUNTER=find((r_ini>precession_r_bins(bin)).*(r_ini<=precession_r_bins(bin+1)).*COUNTER_PASSING_POP.*ALPHAS_POP);
    pphi_Dpphi_value_counter(bin,Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP_COUNTER)));
    pphi_Dpphi_value_co(bin,Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP_CO)));
end
PRECESS_LAMBDA_POP_COUNTER=find(COUNTER_PASSING_POP.*ALPHAS_POP);
PRECESS_LAMBDA_POP_CO=find(CO_PASSING_POP.*ALPHAS_POP);
PRECESS_LAMBDA_POP_CO_POS=find(REGION_POP1.*CO_PASSING_POP.*ALPHAS_POP);
PRECESS_LAMBDA_POP_COUNTER_POS=find(REGION_POP1.*COUNTER_PASSING_POP.*ALPHAS_POP);
PRECESS_LAMBDA_POP_CO_NEG=find(REGION_POP2.*CO_PASSING_POP.*ALPHAS_POP);
PRECESS_LAMBDA_POP_COUNTER_NEG=find(REGION_POP2.*COUNTER_PASSING_POP.*ALPHAS_POP);
g_Ekin_counter(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)));
g_Ekin_co(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_CO)));
Dpphi_Ekin_counter(Ebin)=mean(abs(delta_pphi(PRECESS_LAMBDA_POP_COUNTER)));
Dpphi_Ekin_co(Ebin)=mean(abs(delta_pphi(PRECESS_LAMBDA_POP_CO)));
Dpphi_Ekin_counter_pos(Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP_COUNTER_POS)));
Dpphi_Ekin_co_pos(Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP_CO_POS)));
Dpphi_Ekin_counter_neg(Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP_COUNTER_NEG)));
Dpphi_Ekin_co_neg(Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP_CO_NEG)));



% 
% Ebin=6

Ebin=Ebin+1

load('initialG_flatD_90keV_pre_collapse.mat')
load('initialG_flatD_90keV_precession_stats.mat')
ORBIT_POP=evalin('base',ORBIT_POP_name);
REGION_POP1=(alphas_psi<=psi_core);
REGION_POP2=(alphas_psi>=psi_core).*(alphas_psi<=psi_outer);

pphi_ini=alphas_pphi0;
% r_ini=interp1(1:Nradial,radial_r_value_flux,alphas_psi);
r_ini=r_avg';
lambda_ini=Bavg*alphas_mm./alphas_Ekin;
for bin=1:(length(precession_r_bins)-1)
    PRECESS_LAMBDA_POP_CO=find((r_ini>precession_r_bins(bin)).*(r_ini<=precession_r_bins(bin+1)).*CO_PASSING_POP);
    PRECESS_LAMBDA_POP_COUNTER=find((r_ini>precession_r_bins(bin)).*(r_ini<=precession_r_bins(bin+1)).*COUNTER_PASSING_POP);
    pphi_g_value_counter(bin,Ebin)=(mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)))));
    pphi_g_value_co(bin,Ebin)=(mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_CO)))));
%     pphi_g_value_counter(bin,Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)));
%     pphi_g_value_co(bin,Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)));
    pphi_dr_value(bin,Ebin)=mean(delta_r_avg(PRECESS_LAMBDA_POP_COUNTER))/mean(delta_r_avg(PRECESS_LAMBDA_POP_CO));
end
PRECESS_LAMBDA_POP_COUNTER=find(COUNTER_PASSING_POP);
PRECESS_LAMBDA_POP_CO=find(CO_PASSING_POP);
PRECESS_LAMBDA_POP_COUNTER=find(REGION_POP1.*COUNTER_PASSING_POP);
PRECESS_LAMBDA_POP_CO=find(REGION_POP1.*CO_PASSING_POP);
g_Ekin_counter1(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)));
g_Ekin_co1(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_CO)));
PRECESS_LAMBDA_POP_COUNTER=find(REGION_POP2.*COUNTER_PASSING_POP);
PRECESS_LAMBDA_POP_CO=find(REGION_POP2.*CO_PASSING_POP);
g_Ekin_counter2(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)));
g_Ekin_co2(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_CO)));
dr_Ekin_counter(Ebin)=mean(delta_r_avg(PRECESS_LAMBDA_POP_COUNTER))
dr_Ekin_co(Ebin)=mean(delta_r_avg(PRECESS_LAMBDA_POP_CO))
omega_vD_Ekin_counter(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)));
omega_vD_Ekin_co(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)));
omega_psi_Ekin_counter(Ebin)=mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)));
omega_psi_Ekin_co(Ebin)=mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_CO)));

PRECESS_LAMBDA_POP_COUNTER=find(COUNTER_PASSING_POP.*REGION_POP1);
PRECESS_LAMBDA_POP_CO=find(CO_PASSING_POP.*REGION_POP1);
omega_vD_Ekin_counter1(Ebin)=mean((omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)));
omega_vD_Ekin_co1(Ebin)=mean((omega_precess_avg(PRECESS_LAMBDA_POP_CO)));
omega_psi_Ekin_counter1(Ebin)=mean((omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)));
omega_psi_Ekin_co1(Ebin)=mean((omega_psi_avg(PRECESS_LAMBDA_POP_CO)));
PRECESS_LAMBDA_POP_COUNTER=find(COUNTER_PASSING_POP.*REGION_POP2);
PRECESS_LAMBDA_POP_CO=find(CO_PASSING_POP.*REGION_POP2);
omega_vD_Ekin_counter2(Ebin)=mean((omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)));
omega_vD_Ekin_co2(Ebin)=mean((omega_precess_avg(PRECESS_LAMBDA_POP_CO)));
omega_psi_Ekin_counter2(Ebin)=mean((omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)));
omega_psi_Ekin_co2(Ebin)=mean((omega_psi_avg(PRECESS_LAMBDA_POP_CO)));

load('flatD_90keV_collapse_Glisa_fc2h2_G290114.mat')
lambda_end=Bavg*alphas_mm./alphas_Ekin;
pphi_end=alphas_pphi0;
delta_pphi=pphi_end-pphi_ini; 
delta_lambda=lambda_end-lambda_ini;
ALPHAS_POP=(~alphas_ejected).*(alphas_psi<=psi_outer);
for bin=1:(length(precession_r_bins)-1)
    PRECESS_LAMBDA_POP_CO=find((r_ini>precession_r_bins(bin)).*(r_ini<=precession_r_bins(bin+1)).*CO_PASSING_POP.*ALPHAS_POP);
    PRECESS_LAMBDA_POP_COUNTER=find((r_ini>precession_r_bins(bin)).*(r_ini<=precession_r_bins(bin+1)).*COUNTER_PASSING_POP.*ALPHAS_POP);
    pphi_Dpphi_value_counter(bin,Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP_COUNTER)));
    pphi_Dpphi_value_co(bin,Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP_CO)));
end
PRECESS_LAMBDA_POP_COUNTER=find(COUNTER_PASSING_POP.*ALPHAS_POP);
PRECESS_LAMBDA_POP_CO=find(CO_PASSING_POP.*ALPHAS_POP);
PRECESS_LAMBDA_POP_CO_POS=find(REGION_POP1.*CO_PASSING_POP.*ALPHAS_POP);
PRECESS_LAMBDA_POP_COUNTER_POS=find(REGION_POP1.*COUNTER_PASSING_POP.*ALPHAS_POP);
PRECESS_LAMBDA_POP_CO_NEG=find(REGION_POP2.*CO_PASSING_POP.*ALPHAS_POP);
PRECESS_LAMBDA_POP_COUNTER_NEG=find(REGION_POP2.*COUNTER_PASSING_POP.*ALPHAS_POP);
g_Ekin_counter(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)))
g_Ekin_co(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_CO)))
Dpphi_Ekin_counter(Ebin)=mean(abs(delta_pphi(PRECESS_LAMBDA_POP_COUNTER)));
Dpphi_Ekin_co(Ebin)=mean(abs(delta_pphi(PRECESS_LAMBDA_POP_CO)));
Dpphi_Ekin_counter_pos(Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP_COUNTER_POS)));
Dpphi_Ekin_co_pos(Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP_CO_POS)));
Dpphi_Ekin_counter_neg(Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP_COUNTER_NEG)));
Dpphi_Ekin_co_neg(Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP_CO_NEG)));



% 
% Ebin=Ebin+1
% 
% load('initialG_flat_He_2800keV_pre_collapse.mat')
% load('initialG_flat_He_2800keV_precession_stats.mat')
% ORBIT_POP=evalin('base',ORBIT_POP_name);
% REGION_POP1=(alphas_psi<=psi_core);
% REGION_POP2=(alphas_psi>=psi_core).*(alphas_psi<=psi_outer);
% 
% pphi_ini=alphas_pphi0;
% % r_ini=interp1(1:Nradial,radial_r_value_flux,alphas_psi);
% r_ini=r_avg';
% lambda_ini=Bavg*alphas_mm./alphas_Ekin;
% for bin=1:(length(precession_r_bins)-1)
%     PRECESS_LAMBDA_POP_CO=find((r_ini>precession_r_bins(bin)).*(r_ini<=precession_r_bins(bin+1)).*CO_PASSING_POP);
%     PRECESS_LAMBDA_POP_COUNTER=find((r_ini>precession_r_bins(bin)).*(r_ini<=precession_r_bins(bin+1)).*COUNTER_PASSING_POP);
% %     pphi_g_value_counter(bin,Ebin)=(mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)))));
% %     pphi_g_value_co(bin,Ebin)=(mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_CO)))));
%     pphi_g_value_counter(bin,Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)));
%     pphi_g_value_co(bin,Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)));
%     pphi_dr_value(bin,Ebin)=mean(delta_r_avg(PRECESS_LAMBDA_POP_COUNTER))/mean(delta_r_avg(PRECESS_LAMBDA_POP_CO));
% end
% PRECESS_LAMBDA_POP_COUNTER=find(COUNTER_PASSING_POP);
% PRECESS_LAMBDA_POP_CO=find(CO_PASSING_POP);
% PRECESS_LAMBDA_POP_COUNTER=find(REGION_POP1.*COUNTER_PASSING_POP);
% PRECESS_LAMBDA_POP_CO=find(REGION_POP1.*CO_PASSING_POP);
% g_Ekin_counter1(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)));
% g_Ekin_co1(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_CO)));
% PRECESS_LAMBDA_POP_COUNTER=find(REGION_POP2.*COUNTER_PASSING_POP);
% PRECESS_LAMBDA_POP_CO=find(REGION_POP2.*CO_PASSING_POP);
% g_Ekin_counter2(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)));
% g_Ekin_co2(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_CO)));
% dr_Ekin_counter(Ebin)=mean(delta_r_avg(PRECESS_LAMBDA_POP_COUNTER))
% dr_Ekin_co(Ebin)=mean(delta_r_avg(PRECESS_LAMBDA_POP_CO))
% omega_vD_Ekin_counter(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)));
% omega_vD_Ekin_co(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)));
% omega_psi_Ekin_counter(Ebin)=mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)));
% omega_psi_Ekin_co(Ebin)=mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_CO)));
% 
% PRECESS_LAMBDA_POP_COUNTER=find(COUNTER_PASSING_POP.*REGION_POP1);
% PRECESS_LAMBDA_POP_CO=find(CO_PASSING_POP.*REGION_POP1);
% omega_vD_Ekin_counter1(Ebin)=mean((omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)));
% omega_vD_Ekin_co1(Ebin)=mean((omega_precess_avg(PRECESS_LAMBDA_POP_CO)));
% omega_psi_Ekin_counter1(Ebin)=mean((omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)));
% omega_psi_Ekin_co1(Ebin)=mean((omega_psi_avg(PRECESS_LAMBDA_POP_CO)));
% PRECESS_LAMBDA_POP_COUNTER=find(COUNTER_PASSING_POP.*REGION_POP2);
% PRECESS_LAMBDA_POP_CO=find(CO_PASSING_POP.*REGION_POP2);
% omega_vD_Ekin_counter2(Ebin)=mean((omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)));
% omega_vD_Ekin_co2(Ebin)=mean((omega_precess_avg(PRECESS_LAMBDA_POP_CO)));
% omega_psi_Ekin_counter2(Ebin)=mean((omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)));
% omega_psi_Ekin_co2(Ebin)=mean((omega_psi_avg(PRECESS_LAMBDA_POP_CO)));
% 
% load('flat_He_2800keV_collapse_Glisa_fc1h2_G250314.mat')
% lambda_end=Bavg*alphas_mm./alphas_Ekin;
% pphi_end=alphas_pphi0;
% delta_pphi=pphi_end-pphi_ini;
% delta_lambda=lambda_end-lambda_ini;
% ALPHAS_POP=(~alphas_ejected).*(alphas_psi<=psi_outer);
% for bin=1:(length(precession_r_bins)-1)
%     PRECESS_LAMBDA_POP_CO=find((r_ini>precession_r_bins(bin)).*(r_ini<=precession_r_bins(bin+1)).*CO_PASSING_POP.*ALPHAS_POP);
%     PRECESS_LAMBDA_POP_COUNTER=find((r_ini>precession_r_bins(bin)).*(r_ini<=precession_r_bins(bin+1)).*COUNTER_PASSING_POP.*ALPHAS_POP);
%     pphi_Dpphi_value_counter(bin,Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP_COUNTER)));
%     pphi_Dpphi_value_co(bin,Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP_CO)));
% end
% PRECESS_LAMBDA_POP_COUNTER=find(COUNTER_PASSING_POP.*ALPHAS_POP);
% PRECESS_LAMBDA_POP_CO=find(CO_PASSING_POP.*ALPHAS_POP);
% PRECESS_LAMBDA_POP_CO_POS=find(REGION_POP1.*CO_PASSING_POP.*ALPHAS_POP);
% PRECESS_LAMBDA_POP_COUNTER_POS=find(REGION_POP1.*COUNTER_PASSING_POP.*ALPHAS_POP);
% PRECESS_LAMBDA_POP_CO_NEG=find(REGION_POP2.*CO_PASSING_POP.*ALPHAS_POP);
% PRECESS_LAMBDA_POP_COUNTER_NEG=find(REGION_POP2.*COUNTER_PASSING_POP.*ALPHAS_POP);
% g_Ekin_counter(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)))
% g_Ekin_co(Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_CO)))
% Dpphi_Ekin_counter(Ebin)=mean(abs(delta_pphi(PRECESS_LAMBDA_POP_COUNTER)));
% Dpphi_Ekin_co(Ebin)=mean(abs(delta_pphi(PRECESS_LAMBDA_POP_CO)));
% Dpphi_Ekin_counter_pos(Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP_COUNTER_POS)));
% Dpphi_Ekin_co_pos(Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP_CO_POS)));
% Dpphi_Ekin_counter_neg(Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP_COUNTER_NEG)));
% Dpphi_Ekin_co_neg(Ebin)=mean((delta_pphi(PRECESS_LAMBDA_POP_CO_NEG)));
% 

%%
figure(1)

% max_psi_pos=psi_outer;
% max_r_pos=radial_r_value_flux(max_psi_pos)
% r_Dpphi_value_co=zeros(Nradial,NEbins);
% r_Dpphi_value_counter=zeros(Nradial,NEbins);
% r_g_value_co=zeros(Nradial,NEbins);
% r_g_value_counter=zeros(Nradial,NEbins);
% 
% for Ebin=1:NEbins
%     pphi_range=find(~isnan(pphi_Dpphi_value_co(:,Ebin)));
%     min_pphi=pphi_range(1);
%     max_pphi=pphi_range(end);
%     nb_pphi=max_pphi-min_pphi+1;
%     half_r_bin=round(0.5*(psi_outer-1)/nb_pphi)+1;
%     psi_adapt_scale=psi_scale(nb_pphi*2*half_r_bin:-2*half_r_bin:half_r_bin+1);
% %     psi_adapt_scale=flipud(psi_adapt_scale);
%     r_Dpphi_value_co(:,Ebin)=interp1(psi_adapt_scale,pphi_Dpphi_value_co(min_pphi:max_pphi,Ebin),psi_scale);
%     r_g_value_co(:,Ebin)=interp1(psi_adapt_scale,pphi_g_value_co(min_pphi:max_pphi,Ebin),psi_scale);
%     
%     pphi_range=find(~isnan(pphi_Dpphi_value_counter(:,Ebin)));
%     min_pphi=pphi_range(1);
%     max_pphi=pphi_range(end);
%     nb_pphi=max_pphi-min_pphi+1;
%     half_r_bin=round(0.5*(psi_outer-1)/nb_pphi)+1;
% %     psi_adapt_scale=psi_scale(half_r_bin+1:2*half_r_bin:nb_pphi*2*half_r_bin);
%     psi_adapt_scale=psi_scale(nb_pphi*2*half_r_bin:-2*half_r_bin:half_r_bin+1);
% %     psi_adapt_scale=flipud(psi_adapt_scale);
%     r_Dpphi_value_counter(:,Ebin)=interp1(psi_adapt_scale,pphi_Dpphi_value_counter(min_pphi:max_pphi,Ebin),psi_scale);
%     r_g_value_counter(:,Ebin)=interp1(psi_adapt_scale,pphi_g_value_counter(min_pphi:max_pphi,Ebin),psi_scale);
%     
% end

subplot(2,2,1)
set(gca,'FontSize',22);
% title('Trapped particles');

hold on
grid on


% xlabel('\lambda_{0}')
ylabel('\Deltap_\phi');

title('co passing')

plot(r_bins_values,pphi_Dpphi_value_co(:,2),'g-+','LineWidth',2);
% plot(lambda_values,lambda_Dpphi_value(:,2),'k--','LineWidth',3);
plot(r_bins_values,pphi_Dpphi_value_co(:,3),'b-.','LineWidth',4);
plot(r_bins_values,pphi_Dpphi_value_co(:,4),'r','LineWidth',5);
plot(r_bins_values,pphi_Dpphi_value_co(:,5),'y--','LineWidth',6);
% plot(r_bins_values,r_bins_values.*0+1,'k','LineWidth',2);

% xlim([0.83 1.17])
% legend('\sim 65 keV','400 keV','800 keV','1600 keV','2800 keV')

% xlim([0.2 0.8])




subplot(2,2,3)
set(gca,'FontSize',22);
% title('Trapped particles');

hold on
grid on


% xlabel('\lambda_{0}')
ylabel('|\omega_{v_D}|/|\omega_{\psi}|');


plot(r_bins_values,pphi_g_value_co(:,2),'g-+','LineWidth',2);
% plot(lambda_values,pphi_g_value_counter(:,2),'k--','LineWidth',3);
plot(r_bins_values,pphi_g_value_co(:,3),'b-.','LineWidth',4);
plot(r_bins_values,pphi_g_value_co(:,4),'r','LineWidth',5);
plot(r_bins_values,pphi_g_value_co(:,5),'y--','LineWidth',6);
plot(r_bins_values,r_bins_values.*0+1,'k','LineWidth',2);




figure(1)
subplot(2,2,2)
set(gca,'FontSize',22);
% title('Trapped particles');
title('counter passing')

hold on
grid on


xlabel('\lambda_{0}')
ylabel('\Deltap_\phi');


plot(r_bins_values,pphi_Dpphi_value_counter(:,2),'g-+','LineWidth',2);
% plot(lambda_values,lambda_Dpphi_value(:,2),'k--','LineWidth',3);
plot(r_bins_values,pphi_Dpphi_value_counter(:,3),'b-.','LineWidth',4);
plot(r_bins_values,pphi_Dpphi_value_counter(:,4),'r','LineWidth',5);
plot(r_bins_values,pphi_Dpphi_value_counter(:,5),'y--','LineWidth',6);
% plot(r_bins_values,r_bins_values.*0+1,'k','LineWidth',2);



subplot(2,2,4)
set(gca,'FontSize',22);
% title('Trapped particles');

hold on
grid on


% xlabel('\lambda_{0}')
ylabel('|\omega_{v_D}|/|\omega_{\psi}|');


plot(r_bins_values,pphi_g_value_counter(:,2),'g-+','LineWidth',2);
% plot(lambda_values,lambda_dr_value(:,2),'k--','LineWidth',3);
plot(r_bins_values,pphi_g_value_counter(:,3),'b-.','LineWidth',4);
plot(r_bins_values,pphi_g_value_counter(:,4),'r','LineWidth',5);
plot(r_bins_values,pphi_g_value_counter(:,5),'y--','LineWidth',6);
plot(r_bins_values,r_bins_values.*0+1,'k','LineWidth',2);



%%
figure(2)
subplot(2,1,1)
set(gca,'FontSize',26);
hold on
grid on
% plot(Ekin_bins,Dpphi_Ekin_co_pos(1:end)./Dpphi_Ekin_counter_pos(1:end),'r','LineWidth',4)
% plot(Ekin_bins,Dpphi_Ekin_co_neg(1:end)./Dpphi_Ekin_counter_neg(1:end),'b--','LineWidth',4)
% plot(Ekin_bins,Dpphi_Ekin_counter_pos(1:end),'b--','LineWidth',4)
% plot(Ekin_bins,Dpphi_Ekin_counter_neg(1:end),'b--','LineWidth',4)
% plot(Ekin_bins,Dpphi_Ekin_co_pos(1:end),'r','LineWidth',4)
% plot(Ekin_bins,Dpphi_Ekin_co_neg(1:end),'r','LineWidth',4)
% plot(Ekin_bins,Dpphi_Ekin_counter_pos(1:end),'b--','LineWidth',4)
% plot(Ekin_bins,Dpphi_Ekin_counter_neg(1:end),'b--','LineWidth',4)
% plot([Ecrit Ecrit],[-1 1],'--','color',[0.1 0.7 0.1],'LineWidth',6);
plot([Ecrit Ecrit],[-1 1],'--','color',[0.1 0.7 0.1],'LineWidth',6);

plot(Ekin_bins,Dpphi_Ekin_co_pos(1:end)/Dpphi_Ekin_co_pos(1),'r','LineWidth',4)
plot(Ekin_bins,-Dpphi_Ekin_co_neg(1:end)/Dpphi_Ekin_co_neg(1),'r','LineWidth',4)
plot(Ekin_bins,Dpphi_Ekin_counter_pos(1:end)/Dpphi_Ekin_counter_pos(1),'b--','LineWidth',4)
plot(Ekin_bins,-Dpphi_Ekin_counter_neg(1:end)/Dpphi_Ekin_counter_neg(1),'b--','LineWidth',4)
% set(gca,'XTick',Ekin_v_bins); % Change x-axis ticks
% set(gca,'XTickLabel',Ekin_bins);
% xlim([200 2200])
% xlabel('Ekin');
yl=ylabel('$$\tilde{\Delta p_\varphi }$$','interpreter','latex');
set(yl,'interpreter','latex')
% ylim([-0.7 0.7])
% ylim([-0.23 0.23])



figure(2)
subplot(2,1,2)
set(gca,'FontSize',26);
hold on
grid on

% plot(Ekin_bins,g_Ekin_counter(1:end),'b-.','LineWidth',4)
% plot(Ekin_bins,g_Ekin_co(1:end),'r','LineWidth',4)
plot(Ekin_bins,abs(omega_vD_Ekin_counter1(1:end))./abs(omega_psi_Ekin_counter1(1:end)),'b-.','LineWidth',4)
plot(Ekin_bins,abs(omega_vD_Ekin_co1(1:end))./abs(omega_psi_Ekin_co1(1:end)),'r','LineWidth',4)

plot([Ecrit Ecrit],[-1 1],'--','color',[0.1 0.7 0.1],'LineWidth',6);

plot(Ekin_bins,abs(omega_vD_Ekin_counter2(1:end))./abs(omega_psi_Ekin_counter2(1:end)),'b-.','LineWidth',4)
plot(Ekin_bins,abs(omega_vD_Ekin_co2(1:end))./abs(omega_psi_Ekin_co2(1:end)),'r','LineWidth',4)

% plot(Ekin_bins,omega_vD_Ekin_counter2(1:end)./omega_psi_Ekin_counter2(1:end),'b-.','LineWidth',4)
% plot(Ekin_bins,omega_vD_Ekin_co2(1:end)./omega_psi_Ekin_co2(1:end),'r','LineWidth',4)
% plot(Ekin_bins,g_Ekin_counter(1:end),'b-.','LineWidth',4)
% plot(Ekin_bins,g_Ekin_co(1:end),'r','LineWidth',4)
% set(gca,'XTick',Ekin_v_bins); % Change x-axis ticks
% set(gca,'XTickLabel',Ekin_bins);

% xlim([200 2200])
% ylim([0 0.45])
% xlabel('Ekin (keV)');
ylabel('\omega_{v_D}/\omega_{\psi}');
legend('counter-passing','co-passing')
xlabel('E_{kin} (keV)');

annotation('textbox',...
    [0.52 0.41 0.0483 0.02],...
    'EdgeColor','none',...
    'String',{'E_{cc}'},'fontsize',16,...
    'FitBoxToText','off','color',[0 0.4 0.1]);

annotation('textbox',...
    [0.52 0.85 0.0483 0.02],...
    'EdgeColor','none',...
    'String',{'E_{cc}'},'fontsize',16,...
    'FitBoxToText','off','color',[0 0.4 0.1]);

annotation('textbox',...
    [0.24 0.66 0.0483 0.02],...
    'EdgeColor','none',...
    'String',{'(1)'},'fontsize',18,...
    'FitBoxToText','off','color',[0 0 0]);

annotation('textbox',...
    [0.24 0.84 0.0483 0.02],...
    'EdgeColor','none',...
    'String',{'(2)'},'fontsize',18,...
    'FitBoxToText','off','color',[0 0 0]);

% figure(2)
% subplot(3,1,3)
% set(gca,'FontSize',26);
% hold on
% grid on
% plot(Ekin_bins,omega_vD_Ekin_counter(1:end),'b-.','LineWidth',4)
% plot(Ekin_bins,omega_vD_Ekin_co(1:end),'r','LineWidth',4)
% % plot([ 200 400 800 2200 3200],omega_vD_Ekin_co(2:end)*0+omega_crash,'g','LineWidth',4)
% % set(gca,'XTick',Ekin_v_bins); % Change x-axis ticks
% % set(gca,'XTickLabel',Ekin_bins);
% xlim([200 2200])
% xlabel('E_{kin} (keV)');
% ylabel('|\omega_{v_D}|');
% ylim([0 2.3e5])
% 

%%
figure(3)
set(gca,'FontSize',26);
hold on
grid on
plot(Ekin_bins,omega_vD_Ekin_counter(1:end),'b-.','LineWidth',4)
plot(Ekin_bins,omega_vD_Ekin_co(1:end),'r','LineWidth',4)
xlabel('Ekin (keV)');
% ylabel('\delta_{r}');
ylabel('|\omega_{v_D}|');
