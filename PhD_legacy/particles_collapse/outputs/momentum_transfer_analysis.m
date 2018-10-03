close all
% reset_data_analysis_environment

load('../../data_tokamak/q_profile.mat', 'psi_rank_q1')
psi_mix=size_r-4
delta_psi_q1=20
psi_core=round(psi_rank_q1-delta_psi_q1)
psi_outer=round(1.1*psi_mix)
r_mix=interp1(1:257,radial_r_value_flux,psi_mix);
r_q1=interp1(1:257,radial_r_value_flux,psi_rank_q1);
r_core=interp1(1:257,radial_r_value_flux,psi_core);
r_outer=interp1(1:257,radial_r_value_flux,psi_outer);



rmax=max(radial_r_value_flux);

load('initial_MB_D_m_pre_collapse_all.mat')
load('initial_MB_D_m_precession_stats_all.mat')

psipos_avg=interp1(radial_r_value_flux,1:Nradial,r_avg);
radial_bin_size=30;
radial_bin_half_size=0.5*(radial_bin_size);

radial_bins=[radial_bin_half_size+1:radial_bin_size:simulation_size_r+1*radial_bin_half_size+1]';
N_radial_bins=size(radial_bins,1)
radial_bins_lims=[1:radial_bin_size:radial_bin_size*N_radial_bins+1]';
radial_bins_rho_values=radial_r_value_flux(radial_bins)/rmax;

psipos_profile_ini=zeros(N_radial_bins,1);
phidot_profile_ini=zeros(N_radial_bins,1);
trapped_fraction_profile_ini=zeros(N_radial_bins,1);

phidot_global_ini=mean(phidot_avg(ALL_PASSING_POP|STAGNATION_POP))

for n=1:N_radial_bins
    BINPOP=find((psipos_avg>=radial_bins_lims(n)).*(psipos_avg<radial_bins_lims(n+1)).*(ALL_PASSING_POP|STAGNATION_POP));
    psipos_profile_ini(n)=mean(psipos_avg(BINPOP));
    phidot_profile_ini(n)=mean(phidot_avg(BINPOP));
    BINPOPTOT=find((psipos_avg>=radial_bins_lims(n)).*(psipos_avg<radial_bins_lims(n+1)));
    BINPOPTR=find((psipos_avg>=radial_bins_lims(n)).*(psipos_avg<radial_bins_lims(n+1)).*ALL_TRAPPED_POP);
    trapped_fraction_profile_ini(n)=length(BINPOPTR)/length(BINPOPTOT);
end


%radial bins for toroidal momentum
load('final_MB_D_m_post_collapse_all.mat')
load('final_MB_D_m_precession_stats_all.mat')

psipos_avg=interp1(radial_r_value_flux,1:Nradial,r_avg);

psipos_profile_end=zeros(N_radial_bins,1);
phidot_profile_end=zeros(N_radial_bins,1);
trapped_fraction_profile_end=zeros(N_radial_bins,1);

phidot_global_end=mean(phidot_avg(ALL_PASSING_POP|STAGNATION_POP))

for n=1:N_radial_bins
    BINPOP=find((psipos_avg>=radial_bins_lims(n)).*(psipos_avg<radial_bins_lims(n+1)).*(ALL_PASSING_POP|STAGNATION_POP));
    psipos_profile_end_end(n)=mean(psipos_avg(BINPOP));
    phidot_profile_end(n)=mean(phidot_avg(BINPOP));
    BINPOPTOT=find((psipos_avg>=radial_bins_lims(n)).*(psipos_avg<radial_bins_lims(n+1)));
    BINPOPTR=find((psipos_avg>=radial_bins_lims(n)).*(psipos_avg<radial_bins_lims(n+1)).*ALL_TRAPPED_POP);
    trapped_fraction_profile_end(n)=length(BINPOPTR)/length(BINPOPTOT);
end


%%
figure(1)
set(gca,'fontsize',22)
hold on
grid on
plot(radial_bins_rho_values(1:end-1),-phidot_profile_ini(1:end-1),'b--','linewidth',4)
plot(radial_bins_rho_values(1:end-1),-phidot_profile_end(1:end-1),'r','linewidth',4)
plot([r_q1 r_q1]/rmax,[0 3]*1e4,'k--','LineWidth',4)
% plot([r_mix r_mix]/rmax,[0 3]*1e4,'g--','LineWidth',4)

plot(radial_bins_rho_values(1:end-1),-phidot_profile_ini(1:end-1),'b--','linewidth',4)
plot(radial_bins_rho_values(1:end-1),-phidot_profile_end(1:end-1),'r','linewidth',4)
xlabel('\rho')
ylabel('\Omega (rad/s)')
legend('before sawtooth','after sawtooth')

ylim([1 2.2]*1e4)