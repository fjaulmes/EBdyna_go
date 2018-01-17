
%radial bins for toroidal momentum
radial_bin_size=12;
radial_bin_half_size=0.5*(radial_bin_size);

radial_bins=[radial_bin_half_size+1:radial_bin_size:simulation_size_r+3*radial_bin_half_size+1]';
N_radial_bins=size(radial_bins,1)
radial_bins_lims=[1:radial_bin_size:radial_bin_size*N_radial_bins+1]';

psipos_profile=zeros(N_radial_bins,1);
volume_profile=zeros(N_radial_bins,1);
Npart_ini_profile=zeros(N_radial_bins,1);
Ekin_ini_profile=zeros(N_radial_bins,1);


load('initial_NBI60keV_precession_stats_all.mat')
load('initial_NBI60keV_pre_collapse_all.mat')
psipos_avg_ini=interp1(radial_r_value_flux,1:Nradial,r_avg);

for n=1:N_radial_bins
    BINPOP=find((psipos_avg>=radial_bins_lims(n)).*(psipos_avg<radial_bins_lims(n+1)));
    psipos_profile_ini(n)=mean(psipos_avg_ini(BINPOP));
    volume_profile(n)=sum(sum(volume_tor_diff(radial_bins_lims(n):radial_bins_lims(n+1)-1,:)));
    Npart_ini_profile(n)=length(BINPOP);
    Ekin_ini_profile(n)=mean(alphas_Ekin(BINPOP));
end
close all

% figure(1)
% % plot(psipos_profile,Npart_ini_profile./volume_profile,'b')
% plot(psipos_profile,(2/3)*Ekin_ini_profile,'b')
% hold on




%
psipos_avg_ini=interp1(radial_r_value_flux,1:Nradial,r_avg);

psipos_profile_mom=zeros(N_radial_bins,1);
phidot_profile_mom=zeros(N_radial_bins,1);
trapped_fraction_profile=zeros(N_radial_bins,1);
% current_ini_profile=zeros(N_radial_bins,1);
% momentum_ini_profile=zeros(N_radial_bins,1);
phidot_profile_mom_ini=zeros(N_radial_bins,1);

for n=1:N_radial_bins
    BINPOP=find((psipos_avg_ini>=radial_bins_lims(n)).*(psipos_avg_ini<radial_bins_lims(n+1)));
    psipos_profile_mom(n)=mean(psipos_avg_ini(BINPOP));
    phidot_profile_mom_ini(n)=mean(v_phi(BINPOP)./(alphas_pos_x(BINPOP)+R0));
    psi_ini_profile(n)=mean(alphas_psi_star_ini_value(BINPOP));
%     current_ini_profile(n)=mean(alphas_current(BINPOP));
%     momentum_ini_profile(n)=mean(alphas_momentum(BINPOP));
    BINPOPTOT=find((psipos_avg>=radial_bins_lims(n)).*(psipos_avg<radial_bins_lims(n+1)));
    BINPOPTR=find((psipos_avg>=radial_bins_lims(n)).*(psipos_avg<radial_bins_lims(n+1)).*ALL_TRAPPED_POP);
    trapped_fraction_profile(n)=length(BINPOPTR)/length(BINPOPTOT);
end



density_profile_ini=Npart_ini_profile./volume_profile;


load('final_NBI60keV_precession_stats_all.mat')
load('final_NBI60keV_post_collapse_all.mat')

psipos_avg_end=interp1(radial_r_value_flux,1:Nradial,r_avg);


phidot_profile=zeros(N_radial_bins,1);
psi_final_profile=zeros(N_radial_bins,1);
Npart_end_profile=zeros(N_radial_bins,1);
Ekin_end_profile=zeros(N_radial_bins,1);
% current_final_profile=zeros(N_radial_bins,1);
% momentum_final_profile=zeros(N_radial_bins,1);

phidot_global=mean(phidot_avg(ALL_PASSING_POP|STAGNATION_POP));

for n=1:N_radial_bins
    BINPOP=find((psipos_avg_end>=radial_bins_lims(n)).*(psipos_avg_end<radial_bins_lims(n+1)));
    psipos_profile_end(n)=mean(psipos_avg_end(BINPOP));
    phidot_profile(n)=mean(phidot_avg(BINPOP));
    psi_final_profile(n)=mean(alphas_psi_star_ini_value(BINPOP));
%     current_final_profile(n)=mean(alphas_current(BINPOP));
%     momentum_final_profile(n)=mean(alphas_momentum(BINPOP));
    Npart_end_profile(n)=length(BINPOP);
    Ekin_end_profile(n)=mean(alphas_Ekin(BINPOP));
end
figure(1)
% plot(psipos_profile,Npart_end_profile./volume_profile,'r')
% plot(psipos_profile,(2/3)*Ekin_end_profile,'r')

[psipos_dist_profile hbins]=histc(psipos_avg,radial_bins_lims);

psipos_profile_mom=zeros(N_radial_bins,1);
phidot_profile_mom_end=zeros(N_radial_bins,1);

for n=1:N_radial_bins
    BINPOP=find((psipos_avg_end>=radial_bins_lims(n)).*(psipos_avg_end<radial_bins_lims(n+1)));
    psipos_profile_mom(n)=mean(psipos_avg(BINPOP));
    phidot_profile_mom_end(n)=mean(v_phi(BINPOP)./(alphas_pos_x(BINPOP)+R0));
end

density_profile_end=Npart_end_profile./volume_profile;
