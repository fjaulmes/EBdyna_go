reset_data_analysis_environment
% P_initial_profile=PTOT_profile_interp_ini;
% load('NBI_Phot_data.mat')
load('NBIco_Phot_data.mat', 'deltaW_fast_hat')
deltaW_fast_hat
% rescaling to use only the passing contribution
% deltaW_fast_hat=(1359049/1992454)*deltaW_fast_hat
% save('NBI_Phot_data.mat','-append','deltaW_fast_hat');
% save('deltaW_NBI_data.mat','-append','deltaW_fast_hat');

r_value_q1_mean=interp1(q_initial_profile,radial_r_value_flux,1)


TE_profile_radial_ini=interp1((0:Nradial-1)/(Nradial-1),TE_profile_interp_ini,rho_tor_scale);
TI_profile_radial_ini=interp1((0:Nradial-1)/(Nradial-1),TI_profile_interp_ini,rho_tor_scale);
NE_profile_radial_ini=interp1((0:Nradial-1)/(Nradial-1),NE_profile_interp_ini,rho_tor_scale);
NI_profile_radial_ini=interp1((0:Nradial-1)/(Nradial-1),NI_profile_interp_ini,rho_tor_scale);
PTOT_profile_radial_ini=interp1((0:Nradial-1)/(Nradial-1),PTOT_profile_interp_ini,rho_tor_scale);


%
Te_profile=TE_profile_radial_ini*eV;
Ne_profile=NE_profile_radial_ini;
Ne0=NE_profile_radial_ini(1)
Te0=TE_profile_radial_ini(1)*eV

%
Pbulk_profile=(NI_profile_radial_ini.*TI_profile_radial_ini*eV+NE_profile_radial_ini.*TE_profile_radial_ini*eV);
Ptot_profile=PTOT_profile_radial_ini;
% Ptot_profile=P_initial_profile;

Btot_PR_map=sqrt(BX_PR_map.^2+BZ_PR_map.^2+Btor_PR_map.^2);
derq=gradient(q_initial_profile,radial_r_value_flux);
derP=gradient(Ptot_profile,radial_r_value_flux);
derPi=gradient((NI_profile_radial_ini.*TI_profile_radial_ini*eV),radial_r_value_flux);
derPe=gradient((NE_profile_radial_ini.*TE_profile_radial_ini*eV),radial_r_value_flux);

derTe=gradient(Te_profile,radial_r_value_flux);
derne=gradient(Ne_profile,radial_r_value_flux);
derni=gradient(NI_profile_radial_ini,radial_r_value_flux);

B1=mean(Btot_PR_map(1:end-1,psi_rank_q1))
Ti1=TI_profile_radial_ini(psi_rank_q1)*eV;
Te1=Te_profile(psi_rank_q1);
ne1=Ne_profile(psi_rank_q1)
ni1=NI_profile_radial_ini(psi_rank_q1)
derni1=derni(psi_rank_q1);
derne1=derne(psi_rank_q1);
omegaci=eV*B1/mD;

omegastar_i= (Ti1/(mD*ni1*r_value_q1_mean*omegaci))*derni1
omegastar_e=-(Te1/(mD*ne1*r_value_q1_mean*omegaci))*derne1

vA1=B1/sqrt(mu0*mD*ni1);
derq1=derq(psi_rank_q1);
omegaA1_tilde=vA1/(sqrt(3)*r_value_q1_mean*R0*derq1);
omegaA1=vA1/R0

%be careful, Pi is half of the total pressure - NBI contribution
Pe0=Ne0*Te0
derPi1=derPi(psi_rank_q1);
derPe1=derPe(psi_rank_q1);
derTe1=derTe(psi_rank_q1);
wdia=(1/(mD*ni1*r_value_q1_mean*omegaci))*derPi1
% wdia=(T1/(mD*n1*r_value_q1_mean*omegaci))*dern1


C0=sqrt(1/(epsilon0*mu0));
omegape=sqrt((ne1*eV^2)/(epsilon0*me))
delta_e=C0/omegape

we=-((1/(mD*ne1*r_value_q1_mean*omegaci))*derPe1+0.71*(1/(mD*r_value_q1_mean*omegaci))*derTe1)

% we=-((T1/(mD*r_value_q1_mean*omegaci))*dern1/n1+0.71*(1/(mD*r_value_q1_mean*omegaci))*derT1)
% we=-((0.71*(1/(mD*r_value_q1_mean*omegaci))*derT1)
% we=-((T1/(mD*r_value_q1_mean*omegaci*n1))*dern1)

Bpol1sq=mean(Bpol_PR_map(1:end-1,psi_rank_q1))^2;

dr_avg=mean(dr_PR_map(1:end-1,:));
dr_avg(1:end-1)=radial_r_value_flux(2:end)-radial_r_value_flux(1:end-1);
integ1=0;
for r=2:psi_rank_q1
    integ1=integ1-radial_r_value_flux(r)^2*dr_avg(r)*(derP(r));
end
beta_pol1=(2*mu0/Bpol1sq)*integ1/r_value_q1_mean^2
% beta_pol1=(R0^2*2*mu0/(Bphi0^2*r_value_q1_mean^2))*integ1/r_value_q1_mean^2

%electron ion collision frequency
lambda_D=sqrt(epsilon0*Te_profile./(Ne_profile*eV^2));
Lambda_ei=9*(4*pi/3*Ne_profile.*lambda_D.^3);
log_lambda=log(Lambda_ei);
tau_ei=(3*(2*pi)^1.5)*(sqrt(me)*(Te_profile).^1.5)./(Ne_profile.*log_lambda*eV^2);
tau_ei=(epsilon0/eV)^2*tau_ei;
nu_ei_1=interp1(1:Nradial,1./tau_ei,psi_rank_q1);
Zeff=1.5;
eta1=0.51*Zeff*me*nu_ei_1/(ne1*eV^2)
tau_eta=mu0*r_value_q1_mean^2/eta1;
shear1=r_value_q1_mean*derq1
omegaA1_tilde=vA1/(sqrt(3)*R0*r_value_q1_mean*derq1);
omegaA1_hat=(vA1*r_value_q1_mean*derq1)/(sqrt(3)*R0);
omegaA1_wesson=(vA1*r_value_q1_mean*derq1)/(R0);
tau_H=1/omegaA1_hat;
epsilon_eta=tau_H/tau_eta;

%  gamma Bussac 
deltaq=1-q_initial_profile(1)
r2=interp1(q_initial_profile,radial_r_value_flux,2)
% r2=r_value_q1_mean*(1/deltaq+1)^(1/3)
% beta_coef=(1-(r_value_q1_mean/r2)^2)
% 
% gammaI_nominal=omegaA1*sqrt(3)*pi*(r_value_q1_mean^2/R0^3)*(1/derq1)*deltaq*(3*beta_coef*beta_pol1^2-13/112)
beta_coef=(1-(r_value_q1_mean/r2))
deltaW_ideal_hat=-(1-q_initial_profile(1))*(3*beta_coef*beta_pol1^2-13/112)

% gammaI_nominal=-omegaA1*sqrt(3)*pi*(r_value_q1_mean/R0^2)*(1/derq1)*deltaW_ideal_hat


xvalues=radial_r_value_flux/radial_r_value_flux(psi_rank_q1);
dx_avg=xvalues*0;
dx_avg(2:psi_rank_q1)=xvalues(2:psi_rank_q1)-xvalues(1:psi_rank_q1-1);

% PI_profile=(NI_profile_radial_ini.*TI_profile_radial_ini*eV);
PI_profile=(P_initial_profile-NE_profile_radial_ini.*TE_profile_radial_ini*eV);

%only half of NBI pressure considered contributing to trapped stabilization
PI_profile=0.5*(PI_profile+NI_profile_radial_ini.*TI_profile_radial_ini*eV);

%we now neglect the fast fraction in the trapped contribution
PI_profile=NI_profile_radial_ini.*TI_profile_radial_ini*eV;

cp_integ=0;
Pnorm=PI_profile/PI_profile(1);
for r=2:psi_rank_q1
    cp_integ=cp_integ+((xvalues(r)+xvalues(r-1))/2)^(3/2)*dx_avg(r)*(Pnorm(r));
end

betai0=2*mu0*PI_profile(1)/mean(Btor_PR_map(:,1))^2

deltaW_ko_hat=betai0*(1.5/(6*pi))*cp_integ*(r_value_q1_mean/R0)^(-3/2)

evaluate_li;

deltaW_el_hat=-3*(li_profile(psi_rank_q1)-0.5)^3*(0.5*(kappa1-1))^2

deltaW_hat_tot_el=deltaW_ideal_hat+deltaW_ko_hat+deltaW_fast_hat+deltaW_el_hat
deltaW_hat_tot=deltaW_ideal_hat+deltaW_ko_hat+deltaW_fast_hat

gammaI_nominal=-omegaA1*sqrt(3)*pi*(r_value_q1_mean/R0^2)*(1/derq1)*(deltaW_hat_tot)
deltaW_ideal_ksi0=(deltaW_hat_tot)*6*pi^2*r_value_q1_mean^4*Bphi0^2/(mu0*R0^3)


gammaI_nominal_wdia=gammaI_nominal/abs(wdia)



rhoLi=sqrt(Ti1/mD)/omegaci
rho_s=sqrt(Te1/mD)/omegaci


% gamma_no_coll_porcelli=((2/pi)*(1+Te1/Ti1))^(1/3)*(delta_e/r_value_q1_mean)*(rhoLi/delta_e)^(2/3)
gamma_no_coll_porcelli=((2/pi)*(1+Te1/Ti1))^(1/3)*(delta_e/r_value_q1_mean)*(rho_s/delta_e)^(2/3)*omegaA1_hat

gamma_no_coll_porcelli_dia=((2/pi)*(1+Te1/Ti1))^(1/3)*(delta_e/r_value_q1_mean)*(rhoLi/delta_e)^(2/3)*omegaA1_hat-sqrt(abs(omegastar_e*omegastar_i))

gamma_no_coll_pegoraro=omegaA1_hat*(rhoLi/r_value_q1_mean)^(4/7)*epsilon_eta^(1/7)
gamma_no_coll_pegoraro_dia=omegaA1_hat*(omegaA1_hat^(6)*epsilon_eta*(rhoLi/r_value_q1_mean)^(4)*(1/(abs(omegastar_i)^(1)*abs(omegastar_e)^3*abs(wdia)^(2))))


NORMALIZATION_NOCOLL_GROWTH=omegaA1_hat^(3)*(2*(1+Te1/Ti1))*(rhoLi^2)*delta_e/(pi*r_value_q1_mean^3)

%%
disp('colisionless growth rate')

polynom_matrix=diag([i*omegastar_e i*omegastar_e i*omegastar_e i*omegastar_i i*omegastar_i i*wdia]/gamma_no_coll_porcelli);
polynom_coefs=poly(polynom_matrix);

polynom_coefs(end)=polynom_coefs(end)-1

gamma_solutions=roots(polynom_coefs)

residual_values=polyval(polynom_coefs,gamma_solutions);


sol=gamma_solutions;

residual_sol=(sol-i*omegastar_e/gamma_no_coll_porcelli).*(sol-i*omegastar_i/gamma_no_coll_porcelli).^(2/3).*(sol-i*wdia/gamma_no_coll_porcelli).^(1/3)-1
residual_sol2=(sol-i*omegastar_e/gamma_no_coll_porcelli).^3.*(sol-i*omegastar_i/gamma_no_coll_porcelli).^2.*(sol-i*wdia/gamma_no_coll_porcelli)-1;




%%
disp('semi-collisional growth rate')

RHSeq18=(epsilon_eta*(rhoLi/r_value_q1_mean)^(4))*omegaA1_hat*(omegaA1_hat^(6));
RHSeq18_norm=(epsilon_eta*(rhoLi/r_value_q1_mean)^(4))^(1/7)*omegaA1_hat;
RHSeq18_norm=(((2/pi)*(1+Te1/Ti1))^2*epsilon_eta*(rho_s/r_value_q1_mean)^(4))^(1/7)*omegaA1_hat

polynom_matrix=diag([0 i*omegastar_e i*omegastar_e i*omegastar_e i*omegastar_i i*omegastar_i i*wdia]/RHSeq18_norm);
polynom_coefs=poly(polynom_matrix);

% polynom_coefs(end)=polynom_coefs(end)-(epsilon_eta*(rhoLi/r_value_q1_mean)^(4))
polynom_coefs(end)=polynom_coefs(end)-1

gamma_solutions=roots(polynom_coefs)

% polyval(polynom_coefs,gamma_solutions)


sol=gamma_solutions;

residual_sol=sol.^(1/3).*(sol-i*(omegastar_e/RHSeq18_norm));
residual_sol=residual_sol.*(sol-i*(omegastar_i/RHSeq18_norm)).^(2/3);
residual_sol=residual_sol.*(sol-i*(wdia/RHSeq18_norm)).^(1/3)-1

residual_sol2=sol.*(sol-i*(omegastar_e/RHSeq18_norm)).^3;
residual_sol2=residual_sol2.*(sol-i*(omegastar_i/RHSeq18_norm)).^2;
residual_sol2=residual_sol2.*(sol-i*(wdia/RHSeq18_norm))-1;



%%

disp('-------------------------------------------')
deltaW_ideal_hat
deltaW_ko_hat
deltaW_fast_hat

deltaW_hat_nf_tot=deltaW_ideal_hat+deltaW_ko_hat;

gammaI_nf=-omegaA1*sqrt(3)*pi*(r_value_q1_mean/R0^2)*(1/derq1)*(deltaW_hat_nf_tot)
gammaI_nominal
gammaI_nominal_wdia

disp('-------------------------------------------')
gamma_etak0=RHSeq18_norm
gamma_de=gamma_no_coll_porcelli


%%
disp('-------------------------------------------')
disp('Resistive growth Small Larmor orbit limit')

gamma_eta_small_FLR=epsilon_eta*omegaA1_hat^3/abs(wdia*we)

% disp('-------------------------------------------')
% disp('Ideal growth Small Larmor orbit limit')
% RHSeq333=-((2*sqrt(2))/(pi))*(rhoLi/r_value_q1_mean)*gammaI_nominal*omegaA1_hat;
% RHSeq333_poly=RHSeq333.^2;
% polynom_matrix=diag([0 wdia wdia -we ]/RHSeq333_poly);
% polynom_coefs=poly(polynom_matrix);
% % substract the RHS
% polynom_coefs(end)=polynom_coefs(end)-1
% disp('gamma = imag(roots)')
% gamma_solutions=roots(polynom_coefs)
% sol=gamma_solutions*RHSeq333;
% residual_sol=(sol.*(sol-wdia).*sqrt(1+we./sol))/RHSeq333-1
% 
% residual_sol2=sol.*((sol-wdia).^2).*(sol+we)-RHSeq333_poly

% oiv=i*(-0.2:0.001:2)*abs(wdia);
% eq333=oiv.*(oiv-wdia).*sqrt(1+we./oiv)-RHSeq333;
% 
% [oiv_res oiv_pos]=min(abs(real(eq333)))
% 
% gamma_ideal_small_FLR=imag(oiv(oiv_pos))
% 
% sol=i*gamma_ideal_small_FLR;
% residual_sol=(sol.*(sol-wdia).*sqrt(1+we./sol))/RHSeq333-1

%%
figure(2)
disp('-------------------------------------------')
disp('Small Larmor orbit limit')

delta_layer=delta_e*(0.1:0.1:15);

gamma_small_FLR=sqrt(wdia*we+ (delta_layer/r_value_q1_mean).^2*omegaA1_hat^2);

plot(delta_layer/delta_e,gamma_small_FLR);

xlabel('\delta / d_e')
ylabel('\gamma')

