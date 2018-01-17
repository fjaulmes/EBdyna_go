load('../data_tokamak/psi_profiles_kadomtsev.mat', 'xPsih_zero')
load('ITER_hbeta_m_drive_amplitude.mat')
Btot_PR_map=sqrt(Bpol_PR_map.^2+Btor_PR_map.^2);
Bavg_radial=mean(Btot_PR_map(1:end-1,:),1);


POS_MAX_VA=192

dqdr=gradient(q_initial_profile,radial_r_value_flux);
shear_ini=(radial_r_value_flux./q_initial_profile).*dqdr;
dqdr=gradient(q_final_profile_diff,radial_r_value_flux);
shear_end=(radial_r_value_flux./q_final_profile_diff).*dqdr;

vA_ini=Bavg_radial./sqrt(mu0*Ne_profile*mDT);
vA_end=Bavg_radial./sqrt(mu0*Ne_final_kprofile*mDT);
vA_end(isnan(vA_end))=(vA_end(POS_MAX_VA));
vA_end(isinf(vA_end))=(vA_end(POS_MAX_VA));
vA_end_pol=polyfit(radial_r_value_flux(1:POS_MAX_VA),vA_end(1:POS_MAX_VA),11);
vA_end=polyval(vA_end_pol,radial_r_value_flux);

dqvAdq_ini=gradient(q_initial_profile./vA_ini,q_initial_profile);
dqvAdq_end=gradient(q_final_profile_diff./vA_end,q_final_profile_diff);
epscorr_ini=1./(vA_ini.*dqvAdq_ini);
epscorr_end=1./(vA_end.*dqvAdq_end);

fun_s=[0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.2 1.4 1.6 1.8 2.0 2.2 2.4 2.6 2.8 3.0];
G_fun_s=[0.1366 0.1810 0.2337 0.2968 0.3412 0.4585 0.5594 0.6747 0.9520 1.2970 1.7161 2.2154 2.8010 3.4787 4.2544 5.1336 6.1220 7.2250];
Hp_fun_s=[4.515 2.467 1.270 0.6250 0.3072 0.1560 0.08320 0.04667 0.01686 0.00709 0.00335 0.00173 0.000958 0.000563 0.000347 0.000222 0.000147 0.0001];
Hm_fun_s=[0.09512 0.07606 0.05844 0.04472 0.03455 0.02705 0.02150 0.01733 0.01169 0.00824 0.00602 0.00453 0.00349 0.00275 0.0022 0.00179 0.00148 0.00123];;


NMAX=56
MMAX=12

cont_damping_TAE_ini=zeros(NMAX,MMAX);
cont_damping_TAE_end=zeros(NMAX,MMAX);


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
    m0TAE=n*qTAE+0.5;
    rho_TAE=interp1(q_initial_profile,rho_scale,qTAE);
    s_TAE=interp1(q_initial_profile,shear_ini,qTAE);
    r_TAE=interp1(q_initial_profile,radial_r_value_flux,qTAE);
    psi_pos_TAE_m=interp1(q_initial_profile,1:Nradial,qTAE);
    epsTAE=2*r_TAE/Raxis;
    eps_corr=interp1(q_initial_profile,epscorr_ini,qTAE);
    epsTAE_corr=epsTAE.*eps_corr;
    if s_TAE<=0.3
        Gs_TAE=sqrt(s_TAE/2)*pi^(-1.5);
        Hps_TAE=2./s_TAE;
        Hms_TAE=(pi^2)*s_TAE/8;
    else
        Gs_TAE=interp1(fun_s,G_fun_s,s_TAE);
        Hps_TAE=interp1(fun_s,Hp_fun_s,s_TAE);
        Hms_TAE=interp1(fun_s,Hm_fun_s,s_TAE);
    end
   
    cont_damping_TAE_ini(n,:)=-(epsTAE/2).*(Gs_TAE./(m0TAE.*abs(epsTAE_corr)).^1.5).*...
        (exp(-m0TAE.*abs(epsTAE_corr).*Hps_TAE)+exp(-m0TAE.*abs(epsTAE_corr).*Hms_TAE));

    plot(rho_TAE,cont_damping_TAE_ini(n,:),'color',[0.2 .2 1.0]*n/MAXNPLOT);
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
    m0TAE=n*qTAE+0.5;
    rho_TAE=interp1(q_final_profile_diff,rho_scale,qTAE);
    s_TAE=interp1(q_final_profile_diff,shear_end,qTAE);
    r_TAE=interp1(q_final_profile_diff,radial_r_value_flux,qTAE);
    psi_pos_TAE_m=interp1(q_final_profile_diff,1:Nradial,qTAE);
    epsTAE=2*r_TAE/Raxis;
    eps_corr=interp1(q_final_profile_diff,epscorr_end,qTAE);
    epsTAE_corr=epsTAE.*eps_corr;
    if s_TAE<=0.3
        Gs_TAE=sqrt(s_TAE/2)*pi^(-1.5);
        Hps_TAE=2./s_TAE;
        Hms_TAE=(pi^2)*s_TAE/8;
    else
        Gs_TAE=interp1(fun_s,G_fun_s,s_TAE);
        Hps_TAE=interp1(fun_s,Hp_fun_s,s_TAE);
        Hms_TAE=interp1(fun_s,Hm_fun_s,s_TAE);
    end
   
    cont_damping_TAE_end(n,:)=-(epsTAE/2).*(Gs_TAE./(m0TAE.*abs(epsTAE_corr)).^1.5).*...
        (exp(-m0TAE.*abs(epsTAE_corr).*Hps_TAE)+exp(-m0TAE.*abs(epsTAE_corr).*Hms_TAE));

    plot(rho_TAE,cont_damping_TAE_end(n,:),'color',[0.2 .2 1.0]*n/MAXNPLOT);

end


title('after sawtooth')

xlim([0.2 0.53])
ylim([-0.7 0.4])
xlabel('\rho')
ylabel('\gamma / \omega_{TAE}')


save ITER_hbeta_cont_damping_amplitude.mat cont_damping_TAE_ini cont_damping_TAE_end

