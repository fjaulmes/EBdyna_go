filename=strcat(DATA_FOLDER,'q_profile.mat');
load(filename);
close all
load_particles_init_stats;
NB_PART_RESCALE=3.96e12


PART_SELECT=(~alphas_ejected).*CO_PASSING_POP.*(vpll_avg<vA_TAE/2);
PART_POP=find(PART_SELECT);

psiTAE=interp1(q_initial_profile,psi_scale,qTAE)
pphiTAE=(mHe/eV)*R0*vA_TAE-ZHe*psiTAE

tmax=time_stamp

Dpphi=pphi_output(tmax,:)-pphi_output(1,:);
% Dpphibis=pphi_evol(tmax,:)-pphi_evol(3,:);
DEkin=Ekin_output(tmax,:)-Ekin_output(2,:);
% DEkingc=Ekin_gc_output(tmax,:)-Ekin_gc_output(3,:);
DEpot=Epot_output(tmax,:)-Epot_output(1,:);
DE=Etot_output(tmax,:)-Etot_output(1,:);
DEbis=DEkin+ZHe*DEpot;
% DEth=alphas_Etot_th-Etot_output(2,:)';

RSIZE=0.05

radial_bins_lim=0.1:RSIZE:0.8;
radial_bins_pos=radial_bins_lim(1:end-1)+0.5*RSIZE;
DEkin_r=radial_bins_pos*0;
Dpphi_r=radial_bins_pos*0;
volume_radial_bins=pi*R0*pi*radial_bins_lim.^2;
volume_radial=volume_radial_bins(2:end)-volume_radial_bins(1:end-1);

dn_r_ini=zeros(1,length(radial_bins_pos));
dn_r_end=zeros(1,length(radial_bins_pos));
dn_r_2=zeros(1,length(radial_bins_pos));
dn_r_3=zeros(1,length(radial_bins_pos));
dn_r_4=zeros(1,length(radial_bins_pos));

psi_ini=squeeze(psi_value_output(1,:)');
r_avg_ini=interp1(psi_scale,radial_r_value_flux,psi_ini);
psi_end=squeeze(psi_value_output(end,:)');
r_avg_end=interp1(psi_scale,radial_r_value_flux,psi_end);
psi_2=squeeze(psi_value_output(2,:)');
r_avg_2=interp1(psi_scale,radial_r_value_flux,psi_2);
psi_3=squeeze(psi_value_output(3,:)');
r_avg_3=interp1(psi_scale,radial_r_value_flux,psi_3);
psi_4=squeeze(psi_value_output(4,:)');
r_avg_4=interp1(psi_scale,radial_r_value_flux,psi_4);

for r=1:length(radial_bins_pos)
    RPOP=find(PART_SELECT.*(r_avg_ini>=radial_bins_lim(r)).*(r_avg_ini<radial_bins_lim(r+1)));
    dn_r_ini(r)=sum(alphas_weight(RPOP));
    RPOP=find(PART_SELECT.*(r_avg_end>=radial_bins_lim(r)).*(r_avg_end<radial_bins_lim(r+1)));
    dn_r_end(r)=sum(alphas_weight(RPOP));
    RPOP=find(PART_SELECT.*(r_avg_2>=radial_bins_lim(r)).*(r_avg_2<radial_bins_lim(r+1)));
    dn_r_2(r)=sum(alphas_weight(RPOP));
    RPOP=find(PART_SELECT.*(r_avg_3>=radial_bins_lim(r)).*(r_avg_3<radial_bins_lim(r+1)));
    dn_r_3(r)=sum(alphas_weight(RPOP));
    RPOP=find(PART_SELECT.*(r_avg_4>=radial_bins_lim(r)).*(r_avg_4<radial_bins_lim(r+1)));
    dn_r_4(r)=sum(alphas_weight(RPOP));
end

dn_r_ini=NB_PART_RESCALE*dn_r_ini./volume_radial;
dn_r_end=NB_PART_RESCALE*dn_r_end./volume_radial;
dn_r_2=NB_PART_RESCALE*dn_r_2./volume_radial;
dn_r_3=NB_PART_RESCALE*dn_r_3./volume_radial;
dn_r_4=NB_PART_RESCALE*dn_r_4./volume_radial;

dn_r_ini_dr=-gradient(dn_r_ini,radial_r_value_flux);
dn_r_end_dr=-gradient(dn_r_end,radial_r_value_flux);

%%
figure(2)
set(gca,'fontsize',22)
hold on
grid on
plot(radial_bins_pos,dn_r_ini_dr,'b--','linewidth',3)
% plot(radial_bins_pos,dn_r_2,'b--','linewidth',3)
% plot(radial_bins_pos,dn_r_3,'k--','linewidth',3)
% plot(radial_bins_pos,dn_r_4,'r--','linewidth',3)
plot(radial_bins_pos,dn_r_end_dr,'r','linewidth',3)
% plot([r1 r1],[-1.5 1.5],'g--','linewidth',3)
% plot([r2 r2],[-1.5 1.5],'g--','linewidth',3)
plot([rTAE rTAE],[0 5]*1e18,'g--','linewidth',3)
% ylim([0.2 3.5]*1e17)
xlim([0.25 0.72])
% ylim([0.2 3.5]*1e18)

xlabel('r (m)')
ylabel('-dn_f /dr (m^{-4})')

% legend('initial H 400 keV density','H 400 keV at t=2','H 400 keV at t=3','H 400 keV at t=4','H 400 keV at t=5')
legend('initial D 300 keV','after 12 oscillations')

