% reset_data_analysis_environment;

load('initialG_DT_MB_all_pre_collapse.mat')

% assuming Te == Ti

% load('alphas_vA_all_collapse_Glisa_fc1h2_G110414.mat')
% load('alphas_vA_all_collapse_Glisa_fc2h2_G100214.mat')

Nalphas_simulated=length(find(~alphas_ejected))
volume_vessel=sum(sum(volume_tor_diff));
P_profile=P_initial_profile;
n_profile=Ne_profile;




q_values=interp1(psi_scale,q_initial_profile,psi_bins);
% Ne_values=interp1(psi_scale,Ne_profile,psi_bins);
% B_values=interp1(psi_scale,Ne_profile,psi_bins);

% alphas_vperp=sqrt(2*Bavg*eV*alphas_mm/mDT);
alphas_vperp=sqrt(2*alphas_Ekin*eV/mDT-alphas_vpll.^2);
alphas_speed=sqrt(alphas_vpll.^2+alphas_vperp.^2);


ion_damping_TAE_ini=zeros(NMAX,MMAX);
w_TAE_ini=zeros(NMAX,MMAX);

for n=2:NMAX
    for mindex=1:MMAX
        m=n+mindex-1
        calculate_qTAE_ion_damping;
        ion_damping_TAE_ini(n,mindex)=gamma_TAE;
        w_TAE_ini(n,mindex)=wTAE;
    end
end

save ITER_hbeta_m_ion_damping_amplitude.mat ion_damping_TAE_ini ion_damping_TAE_end w_TAE_ini w_TAE_end

%%

figure(3)
set(gca,'fontsize',18)
imagesc((1:MMAX)-0.5,(1:NMAX),ion_damping_TAE_end./w_TAE_end-ion_damping_TAE_ini./w_TAE_ini)
set(gca,'XTick',[0.5:1:7.5])
colorbar


% set(gca,'XTickLabel',['0.5';'1.5';'2.5';'3.5';'4.5';'5.5';'6.5';'7.5'])
ylabel('n')
xlabel('m=n+val')




figure(4)
set(gca,'fontsize',18)
imagesc((1:MMAX)-0.5,(1:NMAX),(ion_damping_TAE_ini./w_TAE_ini))
set(gca,'XTick',[0.5:1:7.5])
colorbar


% set(gca,'XTickLabel',['0.5';'1.5';'2.5';'3.5';'4.5';'5.5';'6.5';'7.5'])
ylabel('n')
xlabel('m=n+val')