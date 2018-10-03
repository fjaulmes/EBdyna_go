reset_data_analysis_environment;

load('initialG_DT_MB_all_pre_collapse.mat')
load('initialG_DT_MB_all_precession_stats.mat')
load('DT_MBall_collapse_Glisa_fc0p8h2_G160414.mat')
% load('Pkinetic_final_profiles.mat')
P_final_kprofile(isnan(P_final_kprofile))=0

% assuming Te == Ti
alphas_vpll=sqrt(mDT/me)*alphas_vpll;

% load('alphas_vA_all_collapse_Glisa_fc2h2_G100214.mat')

frac_fusion=0.007
Nalphas_simulated=length(find(~alphas_ejected));
volume_vessel=sum(sum(volume_tor_diff));
avg_density_alphas=frac_fusion*Ne0
P_profile=P_final_kprofile;
n_profile=Ne_final_kprofile;
n_bulk_profile=n_profile;
Btot_PR_map=sqrt(Bpol_PR_map.^2+Btor_PR_map.^2);
Btot_radial_profile=mean(Btot_PR_map(1:NP-1,:),1);

PSI_BIN_SIZE=10
PSI_BIN_SIZE_HALF=round(0.5*PSI_BIN_SIZE);

psi_pos_bins=PSI_BIN_SIZE_HALF+1:PSI_BIN_SIZE:Nradial-PSI_BIN_SIZE_HALF;
psi_bins=psi_scale(PSI_BIN_SIZE_HALF+1:PSI_BIN_SIZE:end-PSI_BIN_SIZE_HALF);
NB_PSI_BINS=length(psi_bins)

psi_bins_inf=psi_scale(1:PSI_BIN_SIZE:end-PSI_BIN_SIZE+1);
psi_bins_sup=psi_scale(PSI_BIN_SIZE+1:PSI_BIN_SIZE:end);
psi_bins_lim=(1:PSI_BIN_SIZE:257);

q_values=interp1(psi_scale,q_final_profile_diff,psi_bins);
% Ne_values=interp1(psi_scale,Ne_profile,psi_bins);
% B_values=interp1(psi_scale,Ne_profile,psi_bins);

% alphas_vperp=sqrt(2*Bavg*eV*alphas_mm/mDT);
alphas_vperp=sqrt(2*alphas_Ekin*eV/me-alphas_vpll.^2);
alphas_speed=sqrt(alphas_vpll.^2+alphas_vperp.^2);


NMAX=66
MMAX=22
e_damping_TAE_end=zeros(NMAX,MMAX);
w_TAE_end=zeros(NMAX,MMAX);


VBIN_SIZE=40*1e6;
NB_VBINS=20
V_HALFBIN_SIZE=0.5*VBIN_SIZE;
vperp_bins_lim=(0:VBIN_SIZE:NB_VBINS*VBIN_SIZE);
vperp_bins=(0:VBIN_SIZE:(NB_VBINS-1)*VBIN_SIZE)+V_HALFBIN_SIZE;
VPERP0_BIN_POS=1;
vpll_bins_lim=(-NB_VBINS*VBIN_SIZE:VBIN_SIZE:(NB_VBINS)*VBIN_SIZE)+V_HALFBIN_SIZE;
vpll_bins=(-NB_VBINS*VBIN_SIZE:VBIN_SIZE:(NB_VBINS-1)*VBIN_SIZE)+V_HALFBIN_SIZE;
VPLL0_BIN_POS=round(interp1(vpll_bins,1:length(vpll_bins),0));

DV=3*1e6
SPEED_BIN_SIZE=8*1e6;
NB_SBINS=120
speed_bins_lim=(0:SPEED_BIN_SIZE:NB_SBINS*SPEED_BIN_SIZE);
speed_bins=(0:SPEED_BIN_SIZE:(NB_SBINS-1)*SPEED_BIN_SIZE)+0.5*SPEED_BIN_SIZE;

% vpll_speed_bins_lim=(-NB_SBINS*SPEED_BIN_SIZE:SPEED_BIN_SIZE:NB_SBINS*SPEED_BIN_SIZE);
% vpll_speed_bins=(-NB_SBINS*SPEED_BIN_SIZE:SPEED_BIN_SIZE:(NB_SBINS-1)*SPEED_BIN_SIZE)+0.5*SPEED_BIN_SIZE;

VPERP0_BIN_POS=ceil(interp1(speed_bins,1:length(speed_bins),min(vperp_bins)));
VPERP0_MAX_BIN_POS=floor(interp1(speed_bins,1:length(speed_bins),max(vperp_bins)));


for n=2:NMAX
    for mindex=1:MMAX
        m=n+mindex-1
        calculate_qTAE_e_damping;
        e_damping_TAE_end(n,mindex)=gamma_TAE;
        w_TAE_end(n,mindex)=wTAE;
    end
end

