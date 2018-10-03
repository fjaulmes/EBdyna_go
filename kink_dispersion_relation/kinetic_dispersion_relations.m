
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
omegaA1_porcelli=vA1/(sqrt(3)*R0);
omegaA1=vA1/R0

omegaci=eV*B1/mD;
%be careful, Pi is half of ( total pressure - NBI contribution)
Pe0=Ne0*Te0
frac_Pi=Pe0/P0
derP1=frac_Pi*derP(psi_rank_q1);
derT1=derT(psi_rank_q1);
wdia=(1/(mD*n1*r_value_q1_mean*omegaci))*derP1

C0=sqrt(1/(epsilon0*mu0));
omegape=sqrt((n1*eV^2)/(epsilon0*me))
delta_e=C0/omegape

we=-((1/(mD*n1*r_value_q1_mean*omegaci))*derP1+0.71*(1/(mD*r_value_q1_mean*omegaci))*derT1)

%electron ion collision frequency
lambda_D=sqrt(epsilon0*Te_profile./(Ne_profile*eV^2));
Lambda_ei=9*(4*pi/3*Ne_profile.*lambda_D.^3);
log_lambda=log(Lambda_ei);
tau_ei=(3*(2*pi)^1.5)*(epsilon0^2*sqrt(me)*(Te_profile).^1.5)./(Ne_profile.*log_lambda*eV^4);
% tau_ei=(4*pi)*(epsilon0^2*sqrt(me)*(Te_profile).^1.5)./(Ne_profile.*log_lambda*4*eV^4);
nu_ei_1=interp1(1:Nradial,1./tau_ei,psi_rank_q1)
Zeff=1.5
eta1=0.51*Zeff*me*nu_ei_1/(n1*eV^2)
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




%  gamma Bussac 
Bpol1sq=mean(Bpol_PR_map(1:end-1,psi_rank_q1).^2);
dr_avg=mean(dr_PR_map(1:end-1,:));
dr_avg(1:end-1)=radial_r_value_flux(2:end)-radial_r_value_flux(1:end-1);
integ1=0;
for r=2:psi_rank_q1
    integ1=integ1-radial_r_value_flux(r)^2*dr_avg(r)*(derP(r));
end
beta_pol1=(2*mu0/Bpol1sq)*integ1/r_value_q1_mean^2

deltaq=1-q_initial_profile(1)
r2=interp1(q_initial_profile,radial_r_value_flux,2)

beta_coef=(1-(r_value_q1_mean/r2))
deltaW_ideal_hat=(1-q_initial_profile(1))*(3*beta_coef*beta_pol1^2-13/112)

gammaI_nominal=omegaA1*sqrt(3)*pi*(r_value_q1_mean/R0^2)*(1/derq1)*deltaW_ideal_hat

xvalues=radial_r_value_flux/radial_r_value_flux(psi_rank_q1);
dx_avg=xvalues*0;
dx_avg(2:psi_rank_q1)=xvalues(2:psi_rank_q1)-xvalues(1:psi_rank_q1-1);
cp_integ=0
Pnorm=P_initial_profile/P0;
for r=2:psi_rank_q1
    cp_integ=cp_integ+((xvalues(r)+xvalues(r-1))/2)^(3/2)*dx_avg(r)*(Pnorm(r));
end

betai0=2*mu0*frac_Pi*P0/mean(Btor_PR_map(:,1))^2

deltaW_ko_hat=-betai0*(1.5/(6*pi))*cp_integ*R0*(r_value_q1_mean/R0)^(-3/2)


gammaI_nominal=omegaA1*sqrt(3)*pi*(r_value_q1_mean/R0^2)*(1/derq1)*(deltaW_ideal_hat+deltaW_ko_hat)



%%
figure(2)
hold on

gamma_values=(0.0001:0.0001:1)*abs(wdia);

gamma_I_k=sqrt(omegaA1_hat*(2^(3/2)/pi*rhoLi/r_value_q1_mean*gamma_values));
gamma_eta_k=omegaA1_hat*(rhoLi/r_value_q1_mean)^(4/7)*(4/pi)^(2/7)*epsilon_eta^(1/7);

plot(gamma_values/abs(wdia),gamma_I_k/abs(wdia))
plot(gamma_values/abs(wdia),gamma_I_k*0+gamma_eta_k/abs(wdia),'b--')


    b2=gamma_values.^2-(wdia^2)/4;
    omega_ideal_i1=-sqrt(b2);
    omega_ideal_i2=real(sqrt(b2));

plot(gamma_values/abs(wdia),omega_ideal_i2/abs(wdia),'k--')



%%
% Now loop to solve numerically equation (3.36)
gamma_values=(0.0001:0.0001:2)*abs(wdia);

omegar=gamma_values;
omegai=gamma_values;

om=zeros(length(omegar),length(omegai));
imag_value=zeros(length(omegar),1);
om_solutions=zeros(length(omegar),1);


for n=1:length(omegar)
    om(n,:)=omegar(n)+i*omegai;
end

RHS=-(2*sqrt(2)/pi)*gammaI_nominal*(rhoLi/r_value_q1_mean)*omegaA1_hat
for n=1:length(omegar)
    LHS=om(n,:).*(om(n,:)-wdia).*sqrt(1+we./om(n,:));
    [epsi_solution pos_solution]=min(abs(real(LHS)-RHS));
    imag_value(n)=omegai(pos_solution);
end


% rewrite the values with their known solutions
for n=1:length(omegar)
    om_solutions(n)=omegar(n)+i*imag_value(n);
end