reset_data_analysis_environment;

load('initialG_alphas_vA_all_pre_collapse.mat')
load('alphas_vA_all_collapse_Glisa_fc1h2_G110414.mat')
% load('alphas_vA_all_collapse_Glisa_fc2h2_G100214.mat')

frac_fusion=0.007
Nalphas_simulated=length(find(~alphas_ejected))
volume_vessel=sum(sum(volume_tor_diff));
avg_density_alphas=frac_fusion*Ne0
Frescale=(2/3)*frac_fusion/Nalphas_simulated

Btot_PR_map=sqrt(Bpol_PR_map.^2+Btor_PR_map.^2);
Btot_radial_profile=mean(Btot_PR_map(1:NP-1,:),1);

PSI_BIN_SIZE=10
PSI_BIN_SIZE_HALF=round(0.5*PSI_BIN_SIZE);

psi_pos_bins=PSI_BIN_SIZE_HALF+1:PSI_BIN_SIZE:Nradial-PSI_BIN_SIZE_HALF
psi_bins=psi_scale(PSI_BIN_SIZE_HALF+1:PSI_BIN_SIZE:end-PSI_BIN_SIZE_HALF);
NB_PSI_BINS=length(psi_bins)

psi_bins_inf=psi_scale(1:PSI_BIN_SIZE:end-PSI_BIN_SIZE+1);
psi_bins_sup=psi_scale(PSI_BIN_SIZE+1:PSI_BIN_SIZE:end);
psi_bins_lim=(1:PSI_BIN_SIZE:257);

q_values=interp1(psi_scale,q_final_profile_diff,psi_bins);
% Ne_values=interp1(psi_scale,Ne_profile,psi_bins);
% B_values=interp1(psi_scale,Ne_profile,psi_bins);

% alphas_vperp=sqrt(2*Bavg*eV*alphas_mm/mDT);
alphas_vperp=sqrt(2*alphas_Ekin*eV/mHe-alphas_vpll.^2);


NMAX=42
MMAX=8
drive_TAE_end=zeros(NMAX,MMAX)
w_TAE_end=zeros(NMAX,MMAX)

for n=2:NMAX
    for mindex=1:MMAX
        m=n+mindex-1
        calculate_qTAE_drive;
        drive_TAE_end(n,mindex)=gamma_TAE;
        w_TAE_end(n,mindex)=wTAE;
    end
end

