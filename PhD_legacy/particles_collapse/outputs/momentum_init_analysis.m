close all
rmax=max(radial_r_value_flux);

load('initial_MB_D_pre_collapse_all.mat')
load('initial_MB_D_precession_stats_all.mat')

psipos_avg=interp1(radial_r_value_flux,1:Nradial,r_avg);
radial_bin_size=32;
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
load('initial_MB_D_m_pre_collapse_all.mat')
load('initial_MB_D_m_precession_stats_all.mat')

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
plot(radial_bins_rho_values(1:end-1),-phidot_profile_ini(1:end-1),'g--','linewidth',4)
plot(radial_bins_rho_values(1:end-1),-phidot_profile_end(1:end-1),'b','linewidth',4)
xlabel('\rho')
ylabel('\Omega (rad/s)')
legend('natural rotation','shifted vpll')