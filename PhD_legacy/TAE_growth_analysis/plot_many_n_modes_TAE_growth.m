load('../data_tokamak/psi_profiles_kadomtsev.mat', 'xPsih_zero')
load('ITER_hbeta_m_e_damping_amplitude.mat')

r_mix=interp1(1:Nradial,radial_r_value_flux,xPsih_zero)

figure(5)

subplot(2,1,1)

set(gca,'fontsize',16)
grid on;
hold on
plot([r_value_q1_mean r_value_q1_mean],[-1 1],'g--','linewidth',3)
plot([r_mix r_mix],[-1 1],'g--','linewidth',3)

MAXNPLOT=40

for n=2:2:MAXNPLOT
    mval=(1:8);
    qTAE=(n+mval-0.5)/n;
    rTAE=interp1(q_initial_profile,radial_r_value_flux,qTAE);
    
%     plot(rTAE,gamma_radiative_ini(n,:));
%     plot(rTAE,drive_TAE_ini(n,:)./w_TAE_ini(n,:),'r');
%     plot(rTAE,gamma_radiative_ini(n,:));
    plot(rTAE,ion_damping_TAE_ini(n,:)./w_TAE_ini(n,:)+e_damping_TAE_ini(n,:)./w_TAE_ini(n,:)+gamma_radiative_ini_v2(n,:)+drive_TAE_ini(n,:)./w_TAE_ini(n,:),'color',[0.2 .2 1.0]*n/40,'linewidth',2);
%     plot(rTAE,drive_TAE_ini(n,:)./w_TAE_ini(n,:),'color',[0.1 .1 1.0]*n/40,'linewidth',2);
end
title('before sawtooth')


xlim([0.7 1.2])
ylim([-0.3 0])
xlabel('r (m)')
ylabel('\gamma / \omega_{TAE}')

subplot(2,1,2)
set(gca,'fontsize',16)
grid on;
hold on
plot([r_value_q1_mean r_value_q1_mean],[-1 1],'g--','linewidth',3)
plot([r_mix r_mix],[-1 1],'g--','linewidth',3)

for n=2:2:MAXNPLOT
    mval=(1:8);
    qTAE=(n+mval-0.5)/n;
    rTAE=interp1(q_final_profile_diff,radial_r_value_flux,qTAE);
    
%     plot(rTAE,gamma_radiative_ini(n,:));
%     plot(rTAE,drive_TAE_ini(n,:)./w_TAE_ini(n,:),'r');
%     plot(rTAE,gamma_radiative_ini(n,:));
    plot(rTAE,ion_damping_TAE_ini(n,:)./w_TAE_ini(n,:)+e_damping_TAE_end(n,:)./w_TAE_end(n,:)+gamma_radiative_end_v2(n,:)+drive_TAE_end(n,:)./w_TAE_end(n,:),'color',[1.0 .1 0.1]*n/40,'linewidth',2);
%     plot(rTAE,drive_TAE_end(n,:)./w_TAE_end(n,:),'color',[1.0 .1 0.1]*n/40,'linewidth',2);
end

title('after sawtooth')

xlim([0.7 1.2])
ylim([-0.7 0.4])
xlabel('r (m)')
ylabel('\gamma / \omega_{TAE}')
