initialize_folder_names;
filename=strcat(DATA_FOLDER,'XZsmall_fields_tokamak_pre_collapse.mat');
load(filename);

alphas_Bfield=interp2(scale_X,scale_Z,Btot_XZ_map',alphas_pos_x,alphas_pos_z,'*linear');
alphas_Eperp=alphas_Bfield.*alphas_mm;
alphas_Epll=0.5*(mHe/eV)*alphas_vpll.^2;
alphas_Eperp=abs(alphas_Ekin-alphas_Epll);
alphas_vperp=sqrt(2*alphas_Eperp*eV/mHe);

VBIN_SIZE=0.7*1e5;
V_HALFBIN_SIZE=0.5*VBIN_SIZE;
vperp_bins_lim=(0:VBIN_SIZE:16*VBIN_SIZE);
vperp_bins=(0:VBIN_SIZE:15*VBIN_SIZE)+V_HALFBIN_SIZE;
VPERP0_BIN_POS=1;
vpll_bins=(-15*VBIN_SIZE:VBIN_SIZE:15*VBIN_SIZE)+V_HALFBIN_SIZE;
VPLL0_BIN_POS=round(interp1(vpll_bins,1:length(vpll_bins),0));

close all
figure(3)
grid on
set(gca,'fontsize',20)
hold on
mHist = hist2d ([alphas_vpll, alphas_vperp], vpll_bins, vperp_bins);
contour(vpll_bins(1:end-1)+V_HALFBIN_SIZE,vperp_bins(1:end-1)+V_HALFBIN_SIZE,mHist')

xlim([-vperp_bins(end-4) vperp_bins(end-4)])
ylim([vpll_bins(17) vpll_bins(end-4)])
xlabel('v_{||} (m/s)')
ylabel('v_{perp} (m/s)')

