
Bphi0=Bavg;

%
Btot_PR_map=sqrt(BX_PR_map.^2+BZ_PR_map.^2+Btor_PR_map.^2);
derq=gradient(q_initial_profile,radial_r_value_flux);
derP=gradient(P_initial_profile,radial_r_value_flux);
derT=gradient(Te_profile,radial_r_value_flux);

B1=mean(Btot_PR_map(1:end-1,psi_rank_q1))
P1=Ne_profile(psi_rank_q1)
n1=Ne_profile(psi_rank_q1)
vA1=B1/sqrt(mu0*mD*n1);
derq1=derq(psi_rank_q1);
omegaA1=vA1/(sqrt(3)*r_value_q1_mean*R0*derq1);
omegaA1=vA1/R0

omegaci=eV*B1/mD;
%be careful, Pi is half of the total pressure - NBI contribution
Pe0=Ne0*Te0
frac_Pi=Pe0/P0
derP1=frac_Pi*derP(psi_rank_q1);
derT1=derT(psi_rank_q1);
wdia=(1/(mD*n1*r_value_q1_mean*omegaci))*derP1

we=-((1/(mD*n1*r_value_q1_mean*omegaci))*derP1+0.71*(1/(mD*r_value_q1_mean*omegaci))*derT1)

%electron ion collision frequency
lambda_D=sqrt(epsilon0*Te_profile./(Ne_profile*eV^2));
Lambda_ei=(12*pi*Ne_profile.*lambda_D.^3);
log_lambda=log(Lambda_ei);
tau_ei=(3*(2*pi)^1.5)*(epsilon0^2*sqrt(me)*(Te_profile).^1.5)./(Ne_profile.*log_lambda*4*eV^4);
nu_ei_1=interp1(1:Nradial,1./tau_ei,psi_rank_q1)
Zeff=1.2
eta1=0.51*Zeff*me*nu_ei_1/(n1*eV^2)
tau_eta=mu0*r_value_q1_mean^2/eta1
omegaA1_tilde=vA1/(sqrt(3)*R0*r_value_q1_mean*derq1);
omegaA1_hat=r_value_q1_mean*derq1*vA1/(sqrt(3)*R0);
omegaA1_porcelli=vA1/(sqrt(3)*R0);
tau_A=1/omegaA1_porcelli
tau_H=1/omegaA1_hat
epsilon_eta=(tau_H/tau_eta)

wdia_abs=abs(wdia)

Te1=Te_profile(psi_rank_q1);
rhoLi=sqrt(Te1/mD)/omegaci

Sporcelli=tau_eta/tau_A

gamma_eta_simple=(epsilon_eta*omegaA1_tilde^3)^(1/3)
gamma_eta_consistent=((epsilon_eta*omegaA1_tilde^3))/abs(wdia*we)
gamma_eta_simple_porcelli=(Sporcelli)^(-1/3)*(r_value_q1_mean*derq1)^(2/3)/tau_A
gamma_eta_kinetic=1.1*(rhoLi/r_value_q1_mean)^(4/7)*(Sporcelli)^(-1/7)*(r_value_q1_mean*derq1)^(6/7)/tau_A



%%
b_values=(-1:0.001:1)*8*wdia_abs;

figure(5)
set(gca,'fontsize',20)
grid on
hold on

% plot(b_values/wdia,3*b_values.^5-12*we^2*b_values.^3+b_values.^2*3*cste+we^4*b_values,'linewidth',2)
% a_value=-13063.5
% a_value=-10000
% plot(b_values/wdia,3*a_value^2*b_values-4*a_value*b_values*we+we^2*b_values-b_values.^3,'r','linewidth',2)
% a_value=-3000
% plot(b_values/wdia,3*a_value^2*b_values-4*a_value*b_values*we+we^2*b_values-b_values.^3,'r--','linewidth',2)
a_value=0
% image_disp=a_value.*b_values.*(2*a_value-we-wdia)+b_values.*(a_value-we).*(a_value-wdia)-b_values.^3;
    image_disp=b_values.*(3*a_value.^2-2*a_value*wsum+wdia*we)-b_values.^3;

plot(b_values/wdia_abs,image_disp,'g','linewidth',2)
a_value=(we+wdia)/3
% image_disp=a_value.*b_values.*(2*a_value-we-wdia)+b_values.*(a_value-we).*(a_value-wdia)-b_values.^3;
    image_disp=b_values.*(3*a_value.^2-2*a_value*wsum+wdia*we)-b_values.^3;

plot(b_values/wdia_abs,image_disp,'g--','linewidth',2)
%
a_value=0.6*wdia_abs
% image_disp=a_value.*b_values.*(2*a_value-we-wdia)+b_values.*(a_value-we).*(a_value-wdia)-b_values.^3;
    image_disp=b_values.*(3*a_value.^2-2*a_value*wsum+wdia*we)-b_values.^3;
plot(b_values/wdia_abs,image_disp,'b--','linewidth',2)

a_value=0.8*wdia_abs
% image_disp=a_value.*b_values.*(2*a_value-we-wdia)+b_values.*(a_value-we).*(a_value-wdia)-b_values.^3;
    image_disp=b_values.*(3*a_value.^2-2*a_value*wsum+wdia*we)-b_values.^3;

plot(b_values/wdia_abs,image_disp,'b','linewidth',2)
plot(b_values/wdia_abs,b_values*0-epsilon_eta*omegaA1_hat^3,'k--')
% hl=legend('$\omega_r=-0.2 \omega_{*i}$','$\omega_r=-0.06 \omega_{*i}$','$\omega_r=0$','$\omega_r=0.06 \omega_{*i}$','$\omega_r=0.2 \omega_{*i}$','-$\epsilon_\eta\tilde{\omega_A}^2$')
hl=legend('$\omega_r=0$','$\omega_r=(\omega_{*i}+\hat{\omega_{*e}})/3$','$\omega_r=0.6 |\omega_{*i}|$','$\omega_r=0.8 |\omega_{*i}|$','-$\epsilon_\eta\hat{\omega_A}^3$')
set(hl,'fontsize',20)
set(hl,'interpreter','latex')

ylabel('$\rm{Im}[\omega (\omega-\omega_{*i})(\omega-\hat{\omega_{*e}})]$','interpreter','latex')

xlim([-0.02 0.02])
ylim([-0.7 0.7]*1e12)

% xlabel('$\gamma_\eta/\hat{\omega_{*e}}$','interpreter','latex')
xlabel('$\gamma_\eta/|\omega_{*i}|$','interpreter','latex')