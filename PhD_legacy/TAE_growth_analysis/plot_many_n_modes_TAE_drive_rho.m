load('../data_tokamak/psi_profiles_kadomtsev.mat', 'xPsih_zero')


load('ITER_hbeta_m_drive_amplitude.mat')
PSI_BIN_SIZE=10

alphas_vperp=sqrt(2*alphas_Ekin*eV/mHe-alphas_vpll.^2);
alphas_Eperp=0.5*(mHe/eV)*alphas_vperp.^2;
alphas_Bfiel=alphas_Eperp./alphas_mm;
alphas_Bfield=alphas_Eperp./alphas_mm;
alphas_rhoL=alphas_vperp*mHe./(ZHe*eV*alphas_Bfield);


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

MAXNPLOT=52
MAXMPLOT=12
FLR_damping=zeros(1,MAXMPLOT);

for n=2:2:MAXNPLOT
    mval=(1:MAXMPLOT);
    qTAE=(n+mval-0.5)/n;
    rho_TAE=interp1(q_initial_profile,rho_scale,qTAE);
    r_TAE=interp1(q_initial_profile,radial_r_value_flux,qTAE);
    psi_pos_TAE_m=interp1(q_initial_profile,1:Nradial,qTAE);
    
    %     plot(rTAE,gamma_radiative_ini(n,:));
    %     plot(rTAE,drive_TAE_ini(n,:)./w_TAE_ini(n,:),'r');
    %     plot(rTAE,gamma_radiative_ini(n,:));
    
    for j=1:MAXMPLOT
        Fpsi_pop=(alphas_psi>=psi_pos_TAE_m(j)-0.5*PSI_BIN_SIZE).*(alphas_psi<=psi_pos_TAE_m(j)+0.5*PSI_BIN_SIZE);
        RES_POP_SLICE=find(Fpsi_pop);
        FLR_damping(j)=1/(1+(n-0.5+mval(j)/r_TAE(j))*mean(alphas_rhoL(RES_POP_SLICE)));
    end
    FLR_damping(isnan(FLR_damping))=1;
    plot(rho_TAE,drive_TAE_ini_Ekin(n,:)./w_TAE_ini(n,:),'color',[0.2 .2 1.0]*n/MAXNPLOT,'linewidth',2);
    %     plot(rTAE,drive_TAE_ini(n,:)./w_TAE_ini(n,:),'color',[0.1 .1 1.0]*n/40,'linewidth',2);
end
title('dF/dE')



xlim([0.2 0.53])
ylim([-0.2 0.2])
xlabel('\rho')
ylabel('\gamma / \omega_{TAE}')

subplot(2,1,2)
set(gca,'fontsize',16)
grid on;
hold on
plot([rho_q1 rho_q1],[-1 1],'g--','linewidth',3)
plot([rho_mix rho_mix],[-1 1],'g--','linewidth',3)




FLR_damping=zeros(1,MAXMPLOT);

for n=2:2:MAXNPLOT
    mval=(1:MAXMPLOT);
    qTAE=(n+mval-0.5)/n;
    rho_TAE=interp1(q_initial_profile,rho_scale,qTAE);
    r_TAE=interp1(q_initial_profile,radial_r_value_flux,qTAE);    
    psi_pos_TAE_m=interp1(q_initial_profile,1:Nradial,qTAE);
    
    %     plot(rTAE,gamma_radiative_ini(n,:));
    %     plot(rTAE,drive_TAE_ini(n,:)./w_TAE_ini(n,:),'r');
    %     plot(rTAE,gamma_radiative_ini(n,:));
    
    for j=1:MAXMPLOT
        Fpsi_pop=(alphas_psi>=psi_pos_TAE_m(j)-0.5*PSI_BIN_SIZE).*(alphas_psi<=(psi_pos_TAE_m(j)+0.5*PSI_BIN_SIZE));
        RES_POP_SLICE=find(Fpsi_pop);
        FLR_damping(j)=1/(1+(n-0.5+mval(j)/r_TAE(j))*mean(alphas_rhoL(RES_POP_SLICE)));
    end
    FLR_damping(isnan(FLR_damping))=1;
    
    plot(rho_TAE,drive_TAE_ini_gpsi(n,:)./w_TAE_ini(n,:),'color',[0.2 .2 1.0]*n/MAXNPLOT,'linewidth',2);
%     plot(rTAE,drive_TAE_end(n,:)./w_TAE_end(n,:),'color',[1.0 .1 0.1]*n/40,'linewidth',2);
end


title('dF/d\psi')

xlim([0.2 0.7])
ylim([-0.3 0.3])
xlabel('\rho')
ylabel('\gamma / \omega_{TAE}')







load('alphas_vA_all_collapse_Glisa_fc0p8h2_G160414.mat')
PSI_BIN_SIZE=10

alphas_vperp=sqrt(2*alphas_Ekin*eV/mHe-alphas_vpll.^2);
alphas_Eperp=0.5*(mHe/eV)*alphas_vperp.^2;
alphas_Bfiel=alphas_Eperp./alphas_mm;
alphas_Bfield=alphas_Eperp./alphas_mm;
alphas_rhoL=alphas_vperp*mHe./(ZHe*eV*alphas_Bfield);


r_mix=interp1(1:Nradial,radial_r_value_flux,xPsih_zero)

figure(6)

subplot(2,1,1)

rho_scale=radial_r_value_flux/max(radial_r_value_flux);
rho_q1=interp1(psi_scale,rho_scale,psi_q1)
rho_mix=interp1(radial_r_value_flux,rho_scale,r_mix)

set(gca,'fontsize',16)
grid on;
hold on
plot([rho_q1 rho_q1],[-1 1],'g--','linewidth',3)
plot([rho_mix rho_mix],[-1 1],'g--','linewidth',3)

MAXNPLOT=52
MAXMPLOT=12
FLR_damping=zeros(1,MAXMPLOT);

for n=2:2:MAXNPLOT
    mval=(1:MAXMPLOT);
    qTAE=(n+mval-0.5)/n;
    rho_TAE=interp1(q_final_profile_diff,rho_scale,qTAE);
    r_TAE=interp1(q_final_profile_diff,radial_r_value_flux,qTAE);
    psi_pos_TAE_m=interp1(q_final_profile_diff,1:Nradial,qTAE);
    
    %     plot(rTAE,gamma_radiative_ini(n,:));
    %     plot(rTAE,drive_TAE_ini(n,:)./w_TAE_ini(n,:),'r');
    %     plot(rTAE,gamma_radiative_ini(n,:));
    
    for j=1:MAXMPLOT
        Fpsi_pop=(alphas_psi>=psi_pos_TAE_m(j)-0.5*PSI_BIN_SIZE).*(alphas_psi<=psi_pos_TAE_m(j)+0.5*PSI_BIN_SIZE);
        RES_POP_SLICE=find(Fpsi_pop);
        FLR_damping(j)=1/(1+(n-0.5+mval(j)/r_TAE(j))*mean(alphas_rhoL(RES_POP_SLICE)));
    end
    FLR_damping(isnan(FLR_damping))=1;
    plot(rho_TAE,drive_TAE_end_Ekin(n,:)./w_TAE_end(n,:),'color',[0.2 .2 1.0]*n/MAXNPLOT,'linewidth',2);
    %     plot(rTAE,drive_TAE_ini(n,:)./w_TAE_ini(n,:),'color',[0.1 .1 1.0]*n/40,'linewidth',2);
end
title('dF/dE')



xlim([0.2 0.7])
ylim([-0.2 0.2])
xlabel('\rho')
ylabel('\gamma / \omega_{TAE}')

subplot(2,1,2)
set(gca,'fontsize',16)
grid on;
hold on
plot([rho_q1 rho_q1],[-1 1],'g--','linewidth',3)
plot([rho_mix rho_mix],[-1 1],'g--','linewidth',3)




FLR_damping=zeros(1,MAXMPLOT);

for n=2:2:MAXNPLOT
    mval=(1:MAXMPLOT);
    qTAE=(n+mval-0.5)/n;
    rho_TAE=interp1(q_final_profile_diff,rho_scale,qTAE);
    r_TAE=interp1(q_final_profile_diff,radial_r_value_flux,qTAE);    
    psi_pos_TAE_m=interp1(q_final_profile_diff,1:Nradial,qTAE);
    
    %     plot(rTAE,gamma_radiative_ini(n,:));
    %     plot(rTAE,drive_TAE_ini(n,:)./w_TAE_ini(n,:),'r');
    %     plot(rTAE,gamma_radiative_ini(n,:));
    
    for j=1:MAXMPLOT
        Fpsi_pop=(alphas_psi>=psi_pos_TAE_m(j)-0.5*PSI_BIN_SIZE).*(alphas_psi<=(psi_pos_TAE_m(j)+0.5*PSI_BIN_SIZE));
        RES_POP_SLICE=find(Fpsi_pop);
        FLR_damping(j)=1/(1+(n-0.5+mval(j)/r_TAE(j))*mean(alphas_rhoL(RES_POP_SLICE)));
    end
    FLR_damping(isnan(FLR_damping))=1;
    
    plot(rho_TAE,drive_TAE_end_gpsi(n,:)./w_TAE_end(n,:),'color',[0.2 .2 1.0]*n/MAXNPLOT,'linewidth',2);
%     plot(rTAE,drive_TAE_end(n,:)./w_TAE_end(n,:),'color',[1.0 .1 0.1]*n/40,'linewidth',2);
end


title('dF/d\psi')

xlim([0.2 0.7])
ylim([-0.3 0.3])
xlabel('\rho')
ylabel('\gamma / \omega_{TAE}')
