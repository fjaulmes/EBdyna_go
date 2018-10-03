ZHe=1
mHe=mD

RFILD=2.16;
ZFILD=0.32;

XF=RFILD-1.65;
ZF=0.32-0.07;

DF=0.15

RHOLIBIN=0.002
PANGLEBIN=5

rhoLi_bin_lims=(0:RHOLIBIN:0.05);
pangle_bin_lims=(0:PANGLEBIN:90);
rhoLi_bins=rhoLi_bin_lims(1:end-1)+0.5*RHOLIBIN;
pangle_bins=pangle_bin_lims(1:end-1)+0.5*PANGLEBIN;

load('NBI60keV_losses2_fc1p6h1p6_all.mat','alphas_ejected','alphas_eject_posX','alphas_eject_posZ','alphas_eject_vpll');
load('initial_NBI60keV_pre_collapse_all.mat');
alphas_Bfield=interp2(scale_X,scale_Z,Btot_XZ_map',alphas_pos_x,alphas_pos_z,'*linear');


% only considering losses coming from the sawtooth region
INNER_PART_BOUNDARY=200


FILDPOP=((alphas_psi<INNER_PART_BOUNDARY).*(alphas_eject_posX>XF-DF).*(alphas_eject_posX<XF+DF).*(alphas_eject_posZ>ZF-DF).*(alphas_eject_posZ<ZF+DF));
FILDPART=find(alphas_ejected.*FILDPOP);
total_number_losses_FILD=length(FILDPART)

pitch_angle=acos(-alphas_eject_vpll./sqrt(2*alphas_Ekin*eV/mHe));
pitch_angle=pitch_angle*180/pi;

alphas_Eperp=alphas_Ekin-0.5*(mHe/eV)*alphas_vpll.^2;

alphas_omegaci=ZHe*eV*alphas_Bfield/mHe;
alphas_rhoLi=sqrt(2*alphas_Eperp*eV/mHe)./alphas_omegaci;

hist_classification=zeros(length(rhoLi_bins),length(pangle_bins));

for rb=1:length(rhoLi_bins)
    for pb=1:length(pangle_bins)
        PARTBIN=find(alphas_ejected.*FILDPOP.*(alphas_rhoLi>=rhoLi_bin_lims(rb)).*(alphas_rhoLi<rhoLi_bin_lims(rb+1)).*...
            (pitch_angle>=pangle_bin_lims(pb)).*(pitch_angle<pangle_bin_lims(pb+1)));
        hist_classification(rb,pb)=length(PARTBIN);
    end
end



figure(1);
set(gca,'fontsize',20)
contourf(pangle_bins,rhoLi_bins,hist_classification)
xlabel('pitch angle');
ylabel('Larmor radius (m)');
hold on
grid on

%%
figure(2)
set(gca,'fontsize',20)
hold on
grid on
contour(scale_X,scale_Z,psi_XZsmall_map',10,'k')
plot(alphas_eject_posX(find(alphas_ejected.*(alphas_psi<INNER_PART_BOUNDARY))),alphas_eject_posZ(find(alphas_ejected.*(alphas_psi<INNER_PART_BOUNDARY))),'r.')
plot(alphas_pos_x(find(alphas_ejected.*(alphas_psi<INNER_PART_BOUNDARY))),alphas_pos_z(find(alphas_ejected.*(alphas_psi<INNER_PART_BOUNDARY))),'g.')


xlabel('X (m)');
ylabel('Z (m)');

NBI_eject_posX=alphas_eject_posX(FILDPART);
NBI_eject_posZ=alphas_eject_posZ(FILDPART);
NBI_Ekin=alphas_Ekin(FILDPART);
NBI_vpll=alphas_eject_vpll(FILDPART);
NBI_ini_posX=alphas_pos_x(FILDPART);
NBI_ini_posZ=alphas_pos_z(FILDPART);

save('FILD_lost_NBI_ions_AUG30382_2p5.mat','NBI_eject_posX','NBI_eject_posZ','NBI_Ekin','NBI_vpll','NBI_ini_posX','NBI_ini_posZ')


%%
figure(3)
set(gca,'fontsize',20)
hold on
grid on
contour(scale_X,scale_Z,psi_XZsmall_map',10,'k')
plot(alphas_pos_x(FILDPART),alphas_pos_z(FILDPART),'r.')
plot(alphas_eject_posX(FILDPART),alphas_eject_posZ(FILDPART),'b.')


xlabel('X (m)');
ylabel('Z (m)');

