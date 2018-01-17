
%%
Btot_PR_map=sqrt(BX_PR_map.^2+BZ_PR_map.^2+Btor_PR_map.^2);
derq=gradient(q_initial_profile,radial_r_value_flux);
derP=gradient(P_initial_profile,radial_r_value_flux);

rcore=mean(radial_r_value_flux(1:psi_rank_q1))
r1=r_value_q1_mean
B1=mean(Btot_PR_map(1:end-1,psi_rank_q1));
n1=Ne_profile(psi_rank_q1);
vA1=B1/sqrt(mu0*mDT*n1);
derq1=derq(psi_rank_q1);
omegaA1_tilde=vA1/(sqrt(3)*r_value_q1_mean*R0*derq1);
omegaA1=vA1/R0

omegaci=eV*B1/mDT;
%be careful, Pi is half of the total pressure - NBI contribution
Pe0=Ne0*Te0
frac_Pi=Pe0/P0
derP1=frac_Pi*derP(psi_rank_q1);
wdia=-(1/(mDT*n1*r1*omegaci))*derP1
% 
% Bpol1sq=mean(Bpol_PR_map(1:end-1,psi_rank_q1).^2);
% dr_avg=mean(dr_PR_map(1:end-1,:));
% integ1=0;
% for r=2:psi_rank_q1
%     integ1=integ1-dr_avg(r)*derP(r)*(radial_r_value_flux(r)/r_value_q1_mean)^2;
% end
% beta_pol1=(2*mu0/Bpol1sq)*integ1;
Bpol1sq=mean(Bpol_PR_map(1:end-1,psi_rank_q1).^2);
% Bpol1sq=mean(Bpol_PR_map(1:end-1,psi_rank_q1))^2;
dr_avg=mean(dr_PR_map(1:end-1,:));
dr_avg(1:end-1)=radial_r_value_flux(2:end)-radial_r_value_flux(1:end-1);

integ1=0;
for r=2:psi_rank_q1
    integ1=integ1-radial_r_value_flux(r)^2*dr_avg(r)*(derP(r));
end
beta_pol1=(2*mu0/Bpol1sq)*integ1/r_value_q1_mean^2

%  gamma Bussac 
% gammaI=(2*pi*omegaA1*r_value_q1_mean^2)/(R0^2)*(1-q_initial_profile(1))*(beta_pol1^2-13/144)

% gammaI=omegaA1*sqrt(3)*pi*((r_value_q1_mean/R0)^2)*(1/(derq1*r_value_q1_mean))*(1-q_initial_profile(1))*(beta_pol1^2-13/144)
deltaq=1-q_initial_profile(1);

% r2=r_value_q1_mean*(1/deltaq+1)^(1/3)
r2=interp1(q_initial_profile,radial_r_value_flux,2)

beta_coef=(1-(r_value_q1_mean/r2))

gammaI_nominal=omegaA1*sqrt(3)*pi*(r_value_q1_mean/R0^2)*(1/derq1)*(1-q_initial_profile(1))*(3*beta_coef*beta_pol1^2-13/112)
gammaI_nominal=gammaI_nominal/2

% beta_coef=(1-(r_value_q1_mean/r2)^2);
% gammaI_nominal=omegaA1*sqrt(3)*pi*(r_value_q1_mean/R0^2)*(1/derq1)*(1-q_initial_profile(1))*(beta_pol1^2-13/144)
% gammaI=omegaA1*sqrt(3)*pi*((r_value_q1_mean/R0)^2)*(1/(derq1*r_value_q1_mean))*(1-q_initial_profile(1))*(3*0.8646*beta_pol1^2-13/112)

% back to White  notation
omegaA1_tilde=vA1/(sqrt(3)*r_value_q1_mean*R0*derq1)


figure(1)
set(gca,'fontsize',22)
hold on
grid on

%%
% using approximate gamma Bussac can yield a scan in r1 or in pressure
omega_dm=70*gammaI_nominal

NUMBER_POINTS=20000;
NB_RATIO=4000;
gamma_ratio_values=zeros(NB_RATIO,1);
omega1_values=zeros(NB_RATIO,1);
omega2_values=zeros(NB_RATIO,1);
betah1_values=zeros(NB_RATIO,1);
betah2_values=zeros(NB_RATIO,1);

for grank=1:NB_RATIO
    gamma_ratio_values(grank)=10/(grank);
    gammaI=omega_dm*gamma_ratio_values(grank);
%     omega_dm=gammaI/gamma_ratio_values(grank);
    
    % solving the real side for gammaI<wdia/2
    
    clear omega_value eqton
    
    for n=1:NUMBER_POINTS
        omega_value(n)=wdia+(0.5*omega_dm-wdia)*(n-1)/(NUMBER_POINTS-1);
        eqton(n)=(1/pi)*sqrt(omega_value(n)*(omega_value(n)-wdia))*log(omega_dm/omega_value(n)-1)-gammaI;
    end
    [max_value max_pos ]=max(eqton);
    if max(eqton)>=0
        omega1=interp1(eqton(1:max_pos),omega_value(1:max_pos),0);
        omega2=interp1(eqton(max_pos:end),omega_value(max_pos:NUMBER_POINTS),0);
        omega1_values(grank)=omega1;
        omega2_values(grank)=omega2;
        betah1=(r1/R0)*omega_dm*sqrt(1-wdia/omega1)/(pi*omegaA1_tilde);
        betah2=(r1/R0)*omega_dm*sqrt(1-wdia/omega2)/(pi*omegaA1_tilde);
        betah1_values(grank)=betah1;
        betah2_values(grank)=betah2;
%         plot(omega1_values(grank),gamma_ratio_values(grank),'b.');
%         plot(omega2_values(grank),gamma_ratio_values(grank),'b.');
        plot(betah1_values(grank),gamma_ratio_values(grank),'b.');
        plot(betah2_values(grank),gamma_ratio_values(grank),'b.');
    end
    
    
end
% gammaI=omegaA1*sqrt(3)*pi*((r_value_q1_mean/R0)^2)*(1/(derq1*r_value_q1_mean))*(1-q_initial_profile(1))*(beta_pol1^2-13/144);
% gammaI=omegaA1*sqrt(3)*pi*((r_value_q1_mean/R0)^2)*(1/(derq1*r_value_q1_mean))*(1-q_initial_profile(1))*(3*0.8646*beta_pol1^2-13/112)
% gammaI_nominal=omegaA1*sqrt(3)*pi*(r_value_q1_mean^4/R0^3)*(1/derq1)*(1-q_initial_profile(1))*(3*beta_coef*beta_pol1^2-13/112)

ratio_ITER=gammaI_nominal/omega_dm

xlabel('\beta_h')
ylabel('\gamma_I/\omega_{dfast}')

plot([0 1],[gammaI_nominal gammaI_nominal]/omega_dm,'r--')
betah_nominal=0.55*(P0-2*Pe0)/(0.5*Bavg^2/mu0)
plot([betah_nominal betah_nominal],[0 1],'r--')

xlim([0.002 0.007])
ylim([0 0.1])
% %%
% figure(1)
% set(gca,'fontsize',20)
% hold on
% grid on
%         plot(betah1_values,gamma_ratio_values,'b.');
%         plot(betah2_values,gamma_ratio_values,'b.');
