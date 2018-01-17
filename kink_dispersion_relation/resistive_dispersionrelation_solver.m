
% load('deltaW_NBI_data.mat', 'deltaW_fast_hat')

% Bphi0=Bavg;
Te_profile=TE_profile_interp_ini*eV;
Ne_profile=NE_profile_interp_ini;
Ne0=NE_profile_interp_ini(1)
Te0=TE_profile_interp_ini(1)*eV

%
Btot_PR_map=sqrt(BX_PR_map.^2+BZ_PR_map.^2+Btor_PR_map.^2);
derq=gradient(q_initial_profile,radial_r_value_flux);
derP=gradient(P_initial_profile,radial_r_value_flux);
derPi=gradient((NI_profile_interp_ini.*TI_profile_interp_ini*eV),radial_r_value_flux);
derPe=gradient((NE_profile_interp_ini.*TE_profile_interp_ini*eV),radial_r_value_flux);

derTe=gradient(Te_profile,radial_r_value_flux);
derne=gradient(Ne_profile,radial_r_value_flux);
derni=gradient(NI_profile_interp_ini,radial_r_value_flux);

B1=mean(Btot_PR_map(1:end-1,psi_rank_q1))
Ti1=TI_profile_interp_ini(psi_rank_q1)*eV;
Te1=Te_profile(psi_rank_q1);
ne1=Ne_profile(psi_rank_q1)
ni1=NI_profile_interp_ini(psi_rank_q1)
vA1=B1/sqrt(mu0*mD*ni1);
derq1=derq(psi_rank_q1);
omegaA1_tilde=vA1/(sqrt(3)*r_value_q1_mean*R0*derq1);
omegaA1=vA1/R0

omegaA1_porcelli=vA1/(sqrt(3)*R0);

omegaci=eV*B1/mD;
%be careful, Pi is half of ( total pressure - NBI contribution)


Pe0=Ne0*Te0
frac_Pi=Pe0/P0
derPi1=derPi(psi_rank_q1);
derPe1=derPe(psi_rank_q1);
derTe1=derTe(psi_rank_q1);
wdia=(1/(mD*ni1*r_value_q1_mean*omegaci))*derPi1

C0=sqrt(1/(epsilon0*mu0));
omegape=sqrt((ne1*eV^2)/(epsilon0*me))
delta_e=C0/omegape


% large parallel heat conductivity LOWERS the electron diamagnetic
% frequency

we=-((1/(mD*ne1*r_value_q1_mean*omegaci))*derPe1+0.71*(1/(mD*r_value_q1_mean*omegaci))*derTe1)
% we=-((T1/(mD*r_value_q1_mean*omegaci*n1))*dern1)

%electron ion collision frequency
lambda_D=sqrt(epsilon0*Te_profile./(Ne_profile*eV^2));
Lambda_ei=9*(4*pi/3*Ne_profile.*lambda_D.^3);
log_lambda=log(Lambda_ei);
tau_ei=(3*(2*pi)^1.5*(sqrt(me)*(Te_profile).^1.5)./(Ne_profile.*log_lambda*eV^2))*epsilon0^2/eV^2;
% tau_ei=(4*pi)*(epsilon0^2*sqrt(me)*(Te_profile).^1.5)./(Ne_profile.*log_lambda*4*eV^4);
nu_ei_1=interp1(1:Nradial,1./tau_ei,psi_rank_q1)
Zeff=1.5
eta1=0.51*Zeff*me*nu_ei_1/(ne1*eV^2)
tau_eta=mu0*r_value_q1_mean^2/eta1
% derq1=5*derq1
shear1=r_value_q1_mean*derq1
omegaA1_tilde=vA1/(sqrt(3)*R0*r_value_q1_mean*derq1);
omegaA1_hat=(vA1*r_value_q1_mean*derq1)/(sqrt(3)*R0);
omegaA1_wesson=(vA1*r_value_q1_mean*derq1)/(R0);
% omegaA1_hat=(vA1*r_value_q1_mean*derq1)/(R0);

% it has to be this expression for consistency of lambdaH
% omegaA1_hat=omegaA1_tilde
% tau_H=2*pi/omegaA1_tilde
% epsilon_eta=tau_H/tau_eta
tau_A=1/omegaA1_porcelli
tau_H=1/omegaA1_hat
tau_H_wesson=1/omegaA1_wesson
epsilon_eta=(1/(tau_eta*omegaA1_hat))

wdia_abs=abs(wdia)

Te1=Te_profile(psi_rank_q1);
rhoLi=sqrt(Te1/mD)/omegaci

Sporcelli=tau_eta/tau_A

%gamma_eta_simple=(epsilon_eta*omegaA1_tilde^3)^(1/3)
gamma_eta_wesson=(tau_H_wesson)^(-2/3)*(tau_eta)^(-1/3)
gamma_eta_consistent=(epsilon_eta*omegaA1_hat^3)/abs(wdia*we)

% Porcelli expression identical to wesson chapter 6.11
gamma_eta_simple_porcelli=(Sporcelli)^(-1/3)*(r_value_q1_mean*derq1)^(2/3)/tau_A
gamma_eta_kinetic=1.1*(rhoLi/r_value_q1_mean)^(4/7)*(Sporcelli)^(-1/7)*(r_value_q1_mean*derq1)^(6/7)/tau_A

delta_layer_porcelli=(Sporcelli*r_value_q1_mean*derq1)^(-1/3)*r_value_q1_mean



%%
b_values=(0:0.000025:1)*2*wdia_abs;

figure(5)
set(gca,'fontsize',20)
% grid 
hold on

wsum=wdia+we

for a_rank=1:400
    a_value=(a_rank-201)*0.005*wdia_abs;
    a_values(a_rank)=a_value;
%     imag_disp=a_value.*b_values.*(2*a_value-we-wdia)+b_values.*(a_value-we).*(a_value-wdia)-b_values.^3;
    imag_disp=b_values.*(3*a_value.^2-2*a_value*wsum+wdia*we)-b_values.^3;
    [imag_value gamma_rank]=min(abs(imag_disp+epsilon_eta*omegaA1_hat^3));
    gamma_values(a_rank)=b_values(gamma_rank);
    plot(b_values/wdia_abs,imag_disp)
end


%%
% now solve the real side
% real_disp=a_values.*(a_values-we).*(a_values-wdia)-a_values.*gamma_values.^2-(gamma_values.^2).*(2*a_values-we-wdia);
real_disp=a_values.*(a_values-we).*(a_values-wdia)-3*a_values.*gamma_values.^2+(gamma_values.^2).*(wsum);
[real_value wr_rank]=min(abs(real_disp));


%%
figure(1)
set(gca,'fontsize',20)
grid on
hold on
plot(a_values/wdia_abs,gamma_values/wdia_abs,'b','linewidth',2)

% plot(a_values/wdia,gamma_values*0+gamma_eta_simple/we,'r--','linewidth',2)
plot(a_values/wdia_abs,gamma_values*0+gamma_eta_simple_porcelli/wdia_abs,'r--','linewidth',2)

plot(a_values/wdia_abs,gamma_values*0+gamma_eta_kinetic/wdia_abs,'g--','linewidth',2)

plot([a_values(wr_rank)/wdia_abs a_values(wr_rank)/wdia_abs],[0 2],'b-.','linewidth',3)


legend('dispersion relation','\gamma _{\eta0}','\gamma _{\eta k}','\rm{Re}(\omega)')


xlabel('$\omega_r/|\omega_{*i}|$','interpreter','latex')

% ylabel('$\gamma_\eta/\hat{\omega_{*e}}$','interpreter','latex')
ylabel('$\gamma_\eta/|\omega_{*i}|$','interpreter','latex')
ylim([0.0 1.1])
xlim([-0.8 0.6])


disp('value of real frequency=')
disp(a_values(wr_rank))
disp('=')
disp(a_values(wr_rank)/wdia_abs)

disp('value of growth rate =')
disp(gamma_values(wr_rank))

plot(a_values/wdia_abs,gamma_values*0+gamma_eta_consistent/wdia_abs,'k--','linewidth',2)

