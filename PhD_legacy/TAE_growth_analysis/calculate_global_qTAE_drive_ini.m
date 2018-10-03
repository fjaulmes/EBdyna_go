% reset_data_analysis_environment;

load('initialG_alphas_vA_all_pre_collapse.mat')
load('ITER_hbm_TAE_radiative_damping.mat')
% load('alphas_vA_all_collapse_Glisa_fc1h2_G110414.mat')
% load('alphas_vA_all_collapse_Glisa_fc2h2_G100214.mat')


frac_fusion=0.007
Nalphas_simulated=length(find(~alphas_ejected))
volume_vessel=sum(sum(volume_tor_diff));
avg_density_alphas=frac_fusion*Ne0
P_profile=n_alphas_initial_profile.*T_alphas_initial_profile;
n_profile=n_alphas_initial_profile;
n_bulk_profile=Ne_profile;




q_values=interp1(psi_scale,q_initial_profile,psi_bins);
% Ne_values=interp1(psi_scale,Ne_profile,psi_bins);
% B_values=interp1(psi_scale,Ne_profile,psi_bins);

% alphas_vperp=sqrt(2*Bavg*eV*alphas_mm/mDT);
alphas_vperp=sqrt(2*alphas_Ekin*eV/mHe-alphas_vpll.^2);
alphas_speed=sqrt(alphas_vpll.^2+alphas_vperp.^2);


drive_TAE_ini=zeros(NMAX,MMAX);
w_TAE_ini=zeros(NMAX,MMAX);
drive_TAE_ini_Ekin=zeros(NMAX,MMAX);
drive_TAE_ini_gpsi=zeros(NMAX,MMAX);

for n=2:NMAX
    n
    for mindex=1:MMAX
        m=n+mindex-1;
        calculate_qTAE_drive;
        drive_TAE_ini(n,mindex)=gamma_TAE;
        drive_TAE_ini_Ekin(n,mindex)=gamma_Ekin;
        drive_TAE_ini_gpsi(n,mindex)=gamma_TAE-gamma_Ekin;
        w_TAE_ini(n,mindex)=wTAE;
    end
end

save ITER_hbeta_m_drive_amplitude.mat drive_TAE_ini drive_TAE_end w_TAE_ini w_TAE_end drive_TAE_ini_Ekin drive_TAE_ini_gpsi drive_TAE_end_Ekin drive_TAE_end_gpsi

%%

POS_DRIVE=(drive_TAE_end-gamma_radiative_end>0);
figure(3)
set(gca,'fontsize',18)
imagesc((1:MMAX)-0.5,(1:NMAX),POS_DRIVE.*drive_TAE_end./drive_TAE_ini,[0 12])
set(gca,'XTick',[0.5:1:7.5])
colorbar

% imagesc(POS_DRIVE.*(drive_TAE_end-gamma_radiative_end)./(drive_TAE_ini-gamma_radiative_ini))

% set(gca,'XTickLabel',['0.5';'1.5';'2.5';'3.5';'4.5';'5.5';'6.5';'7.5'])
ylabel('n')
xlabel('m=n+val')

%%
figure(4)
set(gca,'fontsize',18)
imagesc((1:MMAX)-0.5,(1:NMAX),drive_TAE_end-drive_TAE_ini)
set(gca,'XTick',[0.5:1:11.5])
colorbar


% set(gca,'XTickLabel',['0.5';'1.5';'2.5';'3.5';'4.5';'5.5';'6.5';'7.5'])
ylabel('n')
xlabel('m=n+val')
