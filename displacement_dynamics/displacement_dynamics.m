reset_data_analysis_environment
load('../data_tokamak/B_fields.mat', 'Bphi0')
crash_duration=tau_cr
XI0_INI=0.01
mbulk=mD
MAX_DELTA_RADIUS=0.1*ceil(10*max(rx_evol_lin))

der_q=gradient(q_initial_profile,radial_r_value_flux);

ksi_dot=gradient(ksi0_evol_lin,time_scale_lin*crash_duration);
rx_dot=gradient(rx_evol_lin,time_scale_lin*crash_duration);
u_sep=ksi_dot-rx_dot;
r_core=rx_evol_lin-ksi0_evol_lin;
% f1=u_sep./ksi_dot;
% f2=rx_dot./ksi_dot;
% f1=f1.*min(r_value_q1_mean./(2*pi*r_core),1);
% f2=f2.*r_value_q1_mean./(2*pi*rx_evol_lin);

TRANSITION_FRAME=length(Bstar1_evol_lin)
rmix=ksi0_evol_lin(TRANSITION_FRAME)

% this is the transition that we are really interested in
% but the issue here is that Bstar is not mapped for this part 
% so we cannot really use it....
[minval TRANSITION_FRAME_XI]=max(ksi0_evol_lin);


ksi0_evol_lin_ini=ksi0_evol_lin;
rx_evol_lin_ini=rx_evol_lin;

ksi_recalc_evol=ksi0_evol_lin*0;
ksi_dot_recalc=ksi0_evol_lin*0;
rx_recalc_evol=ksi0_evol_lin*0;
rx_dot_recalc_evol=ksi0_evol_lin*0;


%
C0=sqrt(1/(epsilon0*mu0));
Ne1=interp1(1:Nradial,Ne_profile,psi_rank_q1)
Te1=interp1(1:Nradial,Te_profile,psi_rank_q1)
vthe1=sqrt(2*Te1/me)
omegape=sqrt((Ne1*eV^2)/(epsilon0*me))
Ni1=Ne1

omegaci=eV*Bavg/mbulk;
rholi=sqrt(Te1/mbulk)/omegaci
vA1=Bavg/sqrt(mu0*Ne1*mbulk)

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
deltaEini=deltaE;
% deltaE=(E1+E2)-Efinal;
vA_out=sqrt(2*mu0*deltaE)/sqrt(mu0*Ne1*mbulk);
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
Lambda_ei=9*(4*pi/3*Ne_profile.*lambda_D.^3);
log_lambda=log(Lambda_ei);
% tau_ei=(3*(2*pi)^1.5)*(epsilon0^2*sqrt(me)*(Te_profile).^1.5)./(Ne_profile.*log_lambda*eV^4);
% nu_ei_1=interp1(1:Nradial,1./tau_ei,psi_rank_q1);
Zeff=1.2;
% eta1=0.51*Zeff*me*nu_ei_1/(Ne_profile(psi_rank_q1)*eV^2)
eta_profile=0.51*(sqrt(me)*log_lambda*eV^2)./(3*epsilon0^2*(2*pi*Te_profile).^1.5);
eta1=Zeff*interp1(1:Nradial,eta_profile,psi_rank_q1)
 % Bstar_diff_evol=abs(Bstar1_evol_lin-Bstar3_evol_lin);
% jpll_evol=2*Bstar_diff_evol./(mu0*Delta_SP);
% joule_heating=eta1*(jpll_evol.^2.*(2*pi*R0*delta_evol.*Delta_SP));
% deltaE_joules=max(deltaE.*dvolume3_evol_lin-joule_heating,0);
% deltaE=deltaE_joules./dvolume3_evol_lin;

Bstar_diff_evol=abs(Bstar1_evol_lin-Bstar3_evol_lin);
jpll_evol=Bstar_diff_evol./(mu0*delta_evol);
joule_heating=eta1*(jpll_evol.^2.*(2*pi*R0*pi*delta_e^2/4));
deltaE_joules=max(deltaE.*dvolume3_evol_lin-joule_heating,0);
deltaE=deltaE_joules./dvolume3_evol_lin;

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

for (NB_ITERATIONS=1:12)
    evaluate_Delta_from_xi;
    recalculate_xi_evol;
    recalculate_tau_star;
end




close all
%%

% load('core_displacement.mat')
load('xi_evolution.mat')


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
plot(time_scale_xi-time_scale_xi(20),1.35*xi_r_values,'r','linewidth',2)
plot(time_scale_lin*crash_duration,ksi_recalc_evol,'b--','linewidth',2)
% plot(core_xi_data(:,1)-core_xi_data(85,1),0.25/0.4*((core_xi_data(:,4)-core_xi_data(85,4))),'r')
xlim([-0.1 1.2]*crash_duration)

xlabel('time (s)')
ylabel('\xi (m)')

% ylim([0 0.3])
legend('soft Xray','numerical calculation')

save displacement_dynamics.mat crash_duration time_scale_lin ksi_recalc_evol Delta_SP vA_out ksi_dot_recalc delta_evol

%%

figure(4)

%error of the algorithm
plot(2*delta_evol.*vA_out-Delta_SP.*ksi_dot_recalc(1:TRANSITION_FRAME))


reconnection_time=(TRANSITION_FRAME/101)*crash_duration
reconnection_time_zakharov
reconnection_time_wesson
