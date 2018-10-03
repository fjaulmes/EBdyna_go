reset_data_analysis_environment
load('../data_tokamak/B_fields.mat', 'Bphi0')
crash_duration_ini=tau_cr
crash_duration=crash_duration_ini
XI0_INI=0.01
mbulk=mD
MAX_DELTA_RADIUS=0.1*ceil(10*max(rx_evol_lin))

TE_profile_radial_ini=interp1((0:Nradial-1)/(Nradial-1),TE_profile_interp_ini,rho_tor_scale);
TI_profile_radial_ini=interp1((0:Nradial-1)/(Nradial-1),TI_profile_interp_ini,rho_tor_scale);
NE_profile_radial_ini=interp1((0:Nradial-1)/(Nradial-1),NE_profile_interp_ini,rho_tor_scale);
NI_profile_radial_ini=interp1((0:Nradial-1)/(Nradial-1),NI_profile_interp_ini,rho_tor_scale);
PTOT_profile_radial_ini=interp1((0:Nradial-1)/(Nradial-1),PTOT_profile_interp_ini,rho_tor_scale);


der_q=gradient(q_initial_profile,radial_r_value_flux);

ksi_dot=gradient(ksi0_evol_lin,time_scale_lin*crash_duration);
rx_dot=gradient(rx_evol_lin,time_scale_lin*crash_duration);
u_sep=ksi_dot-rx_dot;
r_core=rx_evol_lin-ksi0_evol_lin;
% f1=u_sep./ksi_dot;
% f2=rx_dot./ksi_dot;
% f1=f1.*min(r_value_q1_mean./(2*pi*r_core),1);
% f2=f2.*r_value_q1_mean./(2*pi*rx_evol_lin);

TRANSITION_FRAME=length(Bstar1_evol_lin);
TRANSITION_FRAME=TRANSITION_FRAME+1
rmix=ksi0_evol_lin(TRANSITION_FRAME)
MAX_DELTA_ALLOWED=0.08*rmix

Bstar1_evol_lin(end+1)=Bstar1_evol_lin(end);
Bstar2_evol_lin(end+1)=Bstar2_evol_lin(end);
Bstar3_evol_lin(end+1)=Bstar3_evol_lin(end);

volume1_evol_lin(end+1)=0.5*volume1_evol_lin(end);
volume2_evol_lin(end+1)=1.5*volume2_evol_lin(end)-0.5*volume2_evol_lin(end-1);
volume3_evol_lin(end+1)=volume1_evol_lin(end)+volume2_evol_lin(end);
% this is the transition that we are really interested in
% but the issue here is that Bstar is not mapped for this part 
% so we cannot really use it....
[max_ksi0 TRANSITION_FRAME_XI]=max(ksi0_evol_lin);


ksi0_evol_lin_ini=ksi0_evol_lin;
rx_evol_lin_ini=rx_evol_lin;

ksi_dot_recalc=ksi_dot;
% ksi_dot_recalc=ksi0_evol_lin*0;
rx_recalc_evol=ksi0_evol_lin*0;
rx_dot_recalc_evol=ksi0_evol_lin*0;

Ne_profile=NE_profile_radial_ini;
Ni_profile=NI_profile_radial_ini;
Te_profile=TE_profile_radial_ini*eV;
Ti_profile=TI_profile_radial_ini*eV;
%
C0=sqrt(1/(epsilon0*mu0));
Ne1=interp1(1:Nradial,Ne_profile,psi_rank_q1)
Ni1=interp1(1:Nradial,Ni_profile,psi_rank_q1)
Te1=interp1(1:Nradial,Te_profile,psi_rank_q1)
Ti1=interp1(1:Nradial,Ti_profile,psi_rank_q1)
vthe1=sqrt(2*Te1/me)
omegape=sqrt((Ne1*eV^2)/(epsilon0*me))

omegaci=eV*Bavg/mbulk;
omegace=eV*Bavg/me;
rholi=sqrt(Ti1/mbulk)/omegaci
rhole=sqrt(Te1/me)/omegace
vA1=Bavg/sqrt(mu0*Ni1*mbulk)

dvolume1_evol_lin=gradient(volume1_evol_lin,1:length(volume1_evol_lin));
dvolume2_evol_lin=gradient(volume2_evol_lin,1:length(volume1_evol_lin));
dvolume3_evol_lin=gradient(volume3_evol_lin,1:length(volume1_evol_lin));
dvolume1_evol_lin=abs(dvolume1_evol_lin);

% E1=0.5*(dvolume1_evol_lin./dvolume3_evol_lin).*Bstar1_evol_lin.^2/mu0;
% E2=0.5*(dvolume2_evol_lin./dvolume3_evol_lin).*Bstar2_evol_lin.^2/mu0;
E1=0.5*Bstar1_evol_lin.^2/mu0;
E2=0.5*Bstar2_evol_lin.^2/mu0;
Efinal=0.5*Bstar3_evol_lin.^2/mu0;
E1ini=E1;
E2ini=E2;
Efinalini=Efinal;


% deltaE=(r_core(1:length(E1))/r_value_q1_mean).*(E1);
% deltaE=(f1(1:length(E1)).*E1+f2(1:length(E1)).*E2)-Efinal;
deltaE=0.5*(E1+E2)-Efinal;
% deltaE=(0.8*E1+0.2*E2)-Efinal;
deltaE_vA=deltaE;
% deltaE=(E1+E2)-Efinal;
vA_out=sqrt(2*deltaE_vA)/sqrt(Ni1*mbulk);
vA_out_ini=vA_out;
vA1=max(vA_out_ini)

% l_out=(Bphi0./sqrt(2*mu0*deltaE))*r_value_q1_mean;
l_out=R0/(1-q_initial_profile(1));

tau_star=vA_out*0+mean(r_value_q1_mean./vA_out);

delta_evol=((C0/omegape)*(1./tau_star+(vthe1./l_out)).^(1/2)).*(tau_star.^(1/2));

delta_avg=mean(delta_evol(8:end-2));
delta_e=(C0/omegape);

%initial guess for Delta
Delta_SP=0*delta_evol+r_value_q1_mean;

%resistivity of the layer
lambda_D=sqrt(epsilon0*Te_profile./(Ne_profile*eV^2));
Lambda_ei=(12*pi*Ne_profile.*lambda_D.^3);
log_lambda=log(Lambda_ei);
% tau_ei=(3*(2*pi)^1.5)*(epsilon0^2*sqrt(me)*(Te_profile).^1.5)./(Ne_profile.*log_lambda*eV^4);
% nu_ei_1=interp1(1:Nradial,1./tau_ei,psi_rank_q1);
Zeff=1.8;
% eta1=0.51*Zeff*me*nu_ei_1/(Ne_profile(psi_rank_q1)*eV^2)
eta_profile=0.51*(sqrt(me)*eV*log_lambda)./(3*epsilon0*(2*pi*Te_profile).^1.5);
eta_profile=(eV/epsilon0)*eta_profile;
eta1=Zeff*interp1(1:Nradial,eta_profile,psi_rank_q1)
 % Bstar_diff_evol=abs(Bstar1_evol_lin-Bstar3_evol_lin);
% jpll_evol=2*Bstar_diff_evol./(mu0*Delta_SP);
% joule_heating=eta1*(jpll_evol.^2.*(2*pi*R0*delta_evol.*Delta_SP));
% deltaE_joules=max(deltaE.*dvolume3_evol_lin-joule_heating,0);
% deltaE=deltaE_joules./dvolume3_evol_lin;

Bstar_diff_evol=abs(Bstar1_evol_lin-Bstar3_evol_lin);
jpll_evol=Bstar_diff_evol./(2*mu0*delta_evol);
joule_heating=eta1*(jpll_evol.^2.*(2*pi*R0*delta_e*delta_e));
deltaE_joules=max(deltaE.*dvolume3_evol_lin-joule_heating,0);
deltaE_vA=deltaE_joules./dvolume3_evol_lin;

rho_s=sqrt(2*Te_profile(psi_rank_q1)/mbulk)/omegaci
reconnection_time_zakharov=(r_value_q1_mean/vA1)/(rho_s*der_q(psi_rank_q1))


tau_eta=mu0*r_value_q1_mean^2/eta1;
tau_A=r_value_q1_mean/vA1;
betae=Ne1*Te1/(Bavg^2/(2*mu0))
reconnection_time_wesson=tau_A/sqrt(tau_A/tau_eta+(C0/(r_value_q1_mean*omegape))^2*(1+(betae/(me/mbulk))^(1/2)))

g_xi=ksi_dot(1:TRANSITION_FRAME);
g_xi_error=g_xi*0;
g_xi_error_integ=g_xi_error*0;
delta_error=delta_evol*0;
delta_error_integ=delta_error*0;
vA_out_error=vA_out*0;
vA_out_error_integ=vA_out_error*0;


%%

evaluate_Delta_from_xi;
recalculate_xi_evol;
recalculate_tau_star;

for (NB_ITERATIONS=1:11)
    evaluate_Delta_from_xi;
    recalculate_xi_evol;
    recalculate_tau_star;
end




close all
%%
deltaPE_evol=(PE2_evol_lin-PE3_evol_lin)./(0.5*Delta_SP(1:length(P1_evol_lin)));
deltaPI_evol=(PI2_evol_lin-PI3_evol_lin)./(0.5*Delta_SP(1:length(P1_evol_lin)));
deltaPtot_evol=(PI2_evol_lin+PE2_evol_lin-PE3_evol_lin-PI3_evol_lin)./(0.5*Delta_SP(1:length(P1_evol_lin)));

vE_evol=deltaPE_evol./(-eV*Ne1*Bavg);
vI_evol=deltaPI_evol./(eV*Ni1*Bavg);
% vdia_evol=2*deltaPtot_evol./(eV*(Ne1+Ni1)*Bavg);
vdia_evol=abs(vI_evol-vE_evol);

figure(10)
grid on
hold on
% plot(time_scale_lin(1:length(vE_evol))*crash_duration,vE_evol)
plot(time_scale_lin(1:length(vI_evol))*crash_duration,vdia_evol,'r')
plot(time_scale_lin(1:length(vA_out))*crash_duration,vA_out,'g')

%%
Efield_eta_evol=eta1*jpll_evol;
Efield_eta_evol(TRANSITION_FRAME)=Efield_eta_evol(TRANSITION_FRAME-2);
Efield_eta_evol(TRANSITION_FRAME-1)=Efield_eta_evol(TRANSITION_FRAME-2);
% djpll_dt_evol=gradient(jpll_evol,time_scale_lin(1:length(vA_out)));
Efield_me_evol=(me/(Ne1*eV^2))*(ksi_dot_recalc(1:length(vA_out))./delta_evol(1:length(vA_out))+(1-q_initial_profile(1)).*jpll_evol*vthe1/R0);
Efield_me_evol(TRANSITION_FRAME)=Efield_me_evol(TRANSITION_FRAME-2);
Efield_me_evol(TRANSITION_FRAME-1)=Efield_me_evol(TRANSITION_FRAME-2);

% Efield_me_evol=(me/(Ne1*eV^2))*(djpll_dt_evol+(1-q_initial_profile(1)).*jpll_evol*vthe1/R0);
Efield_pll_evol=Efield_eta_evol+Efield_me_evol;

Ecritical=Ne1*eV^3*log_lambda(psi_rank_q1)/(4*pi*epsilon0^2*me*C0^2)
EDreicer=(2*pi*eV^2*Ne1*log_lambda(psi_rank_q1))/((4*pi*epsilon0)^2*me*vthe1^2/eV)


figure(9)
grid on
hold on
Efield_me_evol(length(vA_out):length(time_scale_lin))=0;
plot(time_scale_lin*crash_duration,Efield_me_evol)

% plot(time_scale_lin(1:length(vA_out))*crash_duration,Efield_eta_evol)

%%

% load('core_displacement.mat')
% load('xi_evolution.mat')
load('tomas_SXR_data.mat')


figure(1)
set(gca,'fontsize',22)
subplot(2,1,1)
grid on
hold on
% plot(time_scale_lin(1:length(vA_out))*4.0*1e-4,deltaE(1:length(vA_out)))
plot(time_scale_lin(1:length(vA_out))*crash_duration,vA_out(1:length(vA_out)))
% xlim([0.2 1.6]*1e-4)
ylabel('v_A')

subplot(2,1,2)
grid on
hold on
plot(time_scale_lin(1:length(vA_out))*crash_duration,delta_evol(1:length(vA_out)))
% xlim([0.2 1.6]*1e-4)
ylabel('\delta')

figure(2)
subplot(2,1,1)
grid on
hold on
% plot(time_scale_lin(1:length(vA_out))*4.0*1e-4,deltaE(1:length(vA_out)))
plot(time_scale_lin(1:length(vA_out))*crash_duration,ksi_dot_recalc(1:length(vA_out)))
% xlim([0.2 1.6]*1e-4)
ylabel('$\dot{\xi}$','interpreter','latex')

subplot(2,1,2)
grid on
hold on
plot(time_scale_lin(1:length(vA_out))*crash_duration,Delta_SP(1:length(vA_out)))
% xlim([0.2 1.6]*1e-4)
ylabel('\Delta')


%%
figure(3)
set(gca,'fontsize',22)

grid on
hold on;
plot(time_core_evol-time_core_evol(680),1.0*r_value_evol,'r','linewidth',2)
ksi_recalc_evol(51:length(time_scale_lin))=max(ksi_recalc_evol);
plot(time_scale_lin*crash_duration,ksi_recalc_evol,'b--','linewidth',2)
%  plot(core_xi_data(:,1)-core_xi_data(85,1),0.25/0.4*((core_xi_data(:,4)-core_xi_data(85,4))),'r')
xlim([-0.1 1.2]*crash_duration)
ylim([0 1.2]*rmix)

xlabel('time (s)')
ylabel('\xi (m)')

% ylim([0 0.3])
legend('soft Xray','numerical calculation')

%%
save displacement_dynamics_sweet_parker.mat crash_duration time_scale_lin ksi_recalc_evol Delta_SP vA_out ksi_dot_recalc delta_evol dvolume3_recalc Efield_eta_evol Efield_me_evol vdia_evol

%%

figure(4)

%error of the algorithm
plot(2*delta_evol.*vA_out-Delta_SP.*ksi_dot_recalc(1:TRANSITION_FRAME))


reconnection_time=(TRANSITION_FRAME/101)*crash_duration
reconnection_time_zakharov
reconnection_time_wesson


%%
figure(8)
set(gca,'fontsize',22)

grid on
hold on;
plot(time_core_evol(20:80)-time_core_evol(30),log(r_value_evol(20:80)),'r','linewidth',2)
plot(time_core_evol(70:200)-time_core_evol(30),log(r_value_evol(70:200)),'r--','linewidth',2)
% plot(time_scale_lin*crash_duration,log(ksi_recalc_evol),'b--','linewidth',2)
xlim([-0.1 1.2]*crash_duration)
xlabel('time (s)')
