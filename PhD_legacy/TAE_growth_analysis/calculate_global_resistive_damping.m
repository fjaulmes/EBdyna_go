load('../data_tokamak/psi_profiles_kadomtsev.mat', 'xPsih_zero')
load('ITER_hbeta_m_drive_amplitude.mat')
Btot_PR_map=sqrt(Bpol_PR_map.^2+Btor_PR_map.^2);
Bavg_radial=mean(Btot_PR_map(1:end-1,:),1);


MAX_Q_VALUE=2.5

dqdr=gradient(q_initial_profile,radial_r_value_flux);
shear_ini=(radial_r_value_flux./q_initial_profile).*dqdr;
dqdr=gradient(q_final_profile_diff,radial_r_value_flux);
shear_end=(radial_r_value_flux./q_final_profile_diff).*dqdr;

vA_ini=Bavg_radial./sqrt(mu0*Ne_profile*mDT);
vA_end=Bavg_radial./sqrt(mu0*Ne_final_kprofile*mDT);
% vA_end(isnan(vA_end))=(vA_end(POS_MAX_VA));
% vA_end(isinf(vA_end))=(vA_end(POS_MAX_VA));
% vA_end_pol=polyfit(radial_r_value_flux(1:POS_MAX_VA),vA_end(1:POS_MAX_VA),11);
% vA_end=polyval(vA_end_pol,radial_r_value_flux);

% LDebye=sqrt(epsilon0*Te_profile./(Ne_profile*eV^2));
% lambda=LDebye.*Te_profile/(eV^2);
% coulomb_log=log(lambda);
coulomb_log_ini=14.9-0.5*log(Ne_profile*1e-20)+log((1e-3)*Te_profile/eV);
tau_ee_ini=(3*(2*pi)^1.5)*(epsilon0^2*sqrt(me)*(Te_profile).^1.5)./(coulomb_log_ini.*Ne_profile*eV^4);
eta_s_ini=0.51*me./(eV^2*Ne_profile.*tau_ee_ini);
coulomb_log_end=14.9-0.5*log(Ne_final_kprofile*1e-20)+log((1e-3)*Te_final_kprofile/eV);
tau_ee_end=(3*(2*pi)^1.5)*(epsilon0^2*sqrt(me)*(Te_final_kprofile).^1.5)./(coulomb_log_end.*Ne_final_kprofile*eV^4);
eta_s_end=0.51*me./(eV^2*Ne_profile.*tau_ee_end);

NMAX=56
MMAX=12

resistive_damping_TAE_ini=zeros(NMAX,MMAX);
resistive_damping_TAE_end=zeros(NMAX,MMAX);


r_mix=interp1(1:Nradial,radial_r_value_flux,xPsih_zero)

figure(5)

subplot(2,1,1)

rho_scale=radial_r_value_flux/max(radial_r_value_flux);
rho_q1=interp1(psi_scale,rho_scale,psi_q1)
rho_mix=interp1(radial_r_value_flux,rho_scale,r_mix)

set(gca,'fontsize',16)
grid on;
hold on
plot([rho_q1 rho_q1],[-1 1],'g--','linewidth',3)
plot([rho_mix rho_mix],[-1 1],'g--','linewidth',3)

MAXNPLOT=56
MAXMPLOT=12
FLR_damping=zeros(1,MAXMPLOT);

for n=2:MAXNPLOT
    mval=(1:MAXMPLOT);
    qTAE=(n+mval-0.5)/n;
    qTAE=min(qTAE,MAX_Q_VALUE);
    m0TAE=n*qTAE+0.5;
    rho_TAE=interp1(q_initial_profile,rho_scale,qTAE);
    s_TAE=interp1(q_initial_profile,shear_ini,qTAE);
    r_TAE=interp1(q_initial_profile,radial_r_value_flux,qTAE);
    psi_pos_TAE_m=interp1(q_initial_profile,1:Nradial,qTAE);
    vA_TAE=interp1(q_initial_profile,vA_ini,qTAE);
    wTAE=(1./(2*R0*qTAE)).*vA_TAE;
    nu_TAE=interp1(q_initial_profile,eta_s_ini,qTAE);
    nu_TAE=nu_TAE.*((n*qTAE./r_TAE).^2)./wTAE;
    
   
    resistive_damping_TAE_ini(n,:)=-0.5*(4*n^2*(s_TAE.^2).*nu_TAE).^(1/3);

    plot(rho_TAE,resistive_damping_TAE_ini(n,:),'color',[0.2 .2 1.0]*n/MAXNPLOT);
    %     plot(rTAE,drive_TAE_ini(n,:)./w_TAE_ini(n,:),'color',[0.1 .1 1.0]*n/40,'linewidth',2);
end
title('before sawtooth')



xlim([0.2 0.53])
ylim([-0.3 0])
xlabel('\rho')
ylabel('\gamma / \omega_{TAE}')



%%
subplot(2,1,2)
set(gca,'fontsize',16)
grid on;
hold on
plot([rho_q1 rho_q1],[-1 1],'g--','linewidth',3)
plot([rho_mix rho_mix],[-1 1],'g--','linewidth',3)





FLR_damping=zeros(1,MAXMPLOT);

for n=2:MAXNPLOT
    mval=(1:MAXMPLOT);
    qTAE=(n+mval-0.5)/n;
    qTAE=min(qTAE,MAX_Q_VALUE);
    m0TAE=n*qTAE+0.5;
    rho_TAE=interp1(q_final_profile_diff,rho_scale,qTAE);
    s_TAE=interp1(q_final_profile_diff,shear_end,qTAE);
    r_TAE=interp1(q_final_profile_diff,radial_r_value_flux,qTAE);
    psi_pos_TAE_m=interp1(q_final_profile_diff,1:Nradial,qTAE);
    vA_TAE=interp1(q_final_profile_diff,vA_end,qTAE);
    wTAE=(1./(2*R0*qTAE)).*vA_TAE;
    nu_TAE=interp1(q_final_profile_diff,eta_s_ini,qTAE);
    nu_TAE=nu_TAE.*((n*qTAE./r_TAE).^2)./wTAE;
    
    resistive_damping_TAE_end(n,:)=-0.5*(4*n^2*(s_TAE.^2).*nu_TAE).^(1/3);

    plot(rho_TAE,resistive_damping_TAE_end(n,:),'color',[0.2 .2 1.0]*n/MAXNPLOT);

end


title('after sawtooth')

xlim([0.2 0.53])
ylim([-0.7 0.4])
xlabel('\rho')
ylabel('\gamma / \omega_{TAE}')


save ITER_hbeta_resistive_damping_amplitude.mat resistive_damping_TAE_ini resistive_damping_TAE_end

