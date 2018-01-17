% reset_data_analysis_environment
rho_tor_sale=sqrt(tor_flux_profile/max(tor_flux_profile));
close all

weight_transp_NBI=1.4e+13

PRE_COLLAPSE_STATS_FILENAME='initial_NBI60keV_R_precession_stats_all.mat'

load(PRE_COLLAPSE_STATS_FILENAME);
% NB_PART_SIM_TRANSP_RATIO=202016/length(r_avg)
% NB_PART_SIM_TRANSP_RATIO=0.5*(0.25+499988/length(r_avg))
NB_PART_SIM_TRANSP_RATIO=(999975/length(r_avg))

NB_PART_SIM_TRANSP_RATIO=1.05*NB_PART_SIM_TRANSP_RATIO



% density_correction=(1e-6)*NB_PART_SIM_TRANSP_RATIO*weight_transp_NBI

%%
PRE_COLLAPSE_FILENAME='initial_NBI60keV_Rlab_pre_collapse_all.mat'
POST_COLLAPSE_FILENAME='NBI60keV_Rlab_fc1h1p6_all.mat'


load(PRE_COLLAPSE_FILENAME);

load('weight_correction_profileP.mat')
alphas_weight_radial_profile_old=alphas_weight_radial_profile;
alphas_weight_corr=interp1(radial_r_value_flux,alphas_weight_radial_profile,r_avg);
alphas_weight=alphas_weight_corr;

%%
% 
RMIN_TOMO=1.5-R0
RMAX_TOMO=1.75-R0
ZMIN_TOMO=-0.4
ZMAX_TOMO=0.4
% RMIN_TOMO=1.65-R0
% RMAX_TOMO=1.9-R0
% ZMIN_TOMO=-0.5
% ZMAX_TOMO=0.5


% EKIN_BIN_SIZE=6*1e3;
% EKIN_BINS=(0:EKIN_BIN_SIZE:84*1e3);
% Ekin_values=EKIN_BINS(1:end-1)+0.5*EKIN_BIN_SIZE;

% EKIN_BIN_SIZE=16*1e3;
EKIN_BINS=[0 20 26 33 34 44]*1e3
Ekin_values=[10 23 29.5 33.5 39 ]*1e3

PICH_BIN_SIZE=0.25;
PITCH_BINS=(-0.8:PICH_BIN_SIZE:1.1);
pitch_values=PITCH_BINS(1:end-1)+0.5*PICH_BIN_SIZE;

R_BIN_SIZE=0.6;
R_BINS=(-0.6:R_BIN_SIZE:0.6)+R0+X_axis;
R_values=R_BINS(1:end-1)+0.5*R_BIN_SIZE;

Z_BIN_SIZE=0.8;
Z_BINS=(-0.8:Z_BIN_SIZE:0.8);
Z_values=Z_BINS(1:end-1)+0.5*Z_BIN_SIZE;

RADIAL_BIN_SIZE=0.1;
% RADIAL_BINS=(0:RADIAL_BIN_SIZE:0.6);
% RADIAL_values=RADIAL_BINS(1:end-1)+0.5*RADIAL_BIN_SIZE;

RADIAL_BINS=[0.00 0.01 0.17 0.25 0.35 0.45]
RADIAL_values=[0.0025 0.09 0.22 0.3 0.4]

psi_radial_values=interp1(radial_r_value_flux,1:Nradial,RADIAL_values);
rho_tor_radial_values=interp1(radial_r_value_flux,rho_tor_sale,RADIAL_values)
volume_radial_bins_values=interp1(radial_r_value_flux,volume_flux,RADIAL_BINS);
volume_radial_values=volume_radial_bins_values(2:end)-volume_radial_bins_values(1:end-1);

%%

% rhist_ini=histc(r_avg,RADIAL_BINS);


l1=length(Ekin_values); 
l2=length(pitch_values)
l3=length(R_values) 
l4=length(Z_values)
l5=length(RADIAL_values)

trapped_fraction_profile_ini=zeros(l5,1);
trapped_fraction_profile_end=zeros(l5,1);
lambda_profile_ini=zeros(l5,1);
lambda_profile_end=zeros(l5,1);


% density_correction_factor=(1e-6)*NB_PART_SIM_TRANSP_RATIO*weight_transp_NBI/(EKIN_BIN_SIZE)/PICH_BIN_SIZE

volume_RZ=zeros(length(R_values),length(Z_values));


for x=1:length(R_values)
    for z=1:length(Z_values)
        volume_RZ(x,z) = Z_BIN_SIZE*R_BIN_SIZE*2*pi*R_values(x);
    end
end

%R_BIN_SIZE=0.06;
%R_BINS=(0:R_BIN_SIZE:0.6);
%r_values=R_BINS(1:end-1)+0.5*R_BIN_SIZE;

alphas_vpll_ini=alphas_vpll;

alphas_B_ini=interp2(scale_X,scale_Z,Btot_XZ_map',alphas_pos_x,alphas_pos_z,'*linear');
alphas_Eperp=max(alphas_Ekin-0.5*(mHe/eV)*alphas_vpll.^2,0);
alphas_vperp=sqrt(2*(eV/mHe)*alphas_Eperp);
alphas_vtot=sqrt(2*(eV/mHe)*alphas_Ekin);
alphas_pitch=-alphas_vpll./alphas_vtot;
alphas_lambda_ini=atan(alphas_vperp./alphas_vpll);
alphas_lambda0_ini=Bavg*alphas_mm./alphas_Ekin;
alphas_r=interp1(1:257,radial_r_value_flux,alphas_psi);
alphas_Ekin_ini=alphas_Ekin;

% here can try alphas_r or r_avg or even pphi [if small range of Ekin]
poloidal_hist_ini=zeros(length(Ekin_values),length(pitch_values),length(R_values),length(Z_values));
dist_EpitchRZ_ini=poloidal_hist_ini;

SAVE_INI=0
SAVE_FINAL=0


% FILENAME_INI='dist_NBI_PEcorr_AUG30815_2p29_EpitchRZ_ini.dat'
% FILENAME_END='dist_NBI_PEcorr_AUG30815_2p29_EpitchRZ_125_end.dat'
% 
% string_title='# initial NBI distribution : dimensions ; scales (4 lines) ; [Ekin pitch R Z] '
% 
Raxis=R0+X_axis

%save dist_EpitchRZ_ini.dat -append  string_title -ASCII
% if SAVE_INI==1
%     save(FILENAME_INI,'-append','Raxis','Z_axis','-ASCII')
%     save(FILENAME_INI,'-append','l1','l2','l3','l4','-ASCII')
%     save(FILENAME_INI,'-append','Ekin_values','-ASCII')
%     save(FILENAME_INI,'-append','pitch_values','-ASCII')
%     save(FILENAME_INI,'-append','R_values','-ASCII')
%     save(FILENAME_INI,'-append','Z_values','-ASCII')
% end
% 
% for eb=1:length(Ekin_values)
%     for pb=1:length(pitch_values)
%         PART_POP=find((alphas_Ekin>=EKIN_BINS(eb)).*(alphas_Ekin<EKIN_BINS(eb+1)).*(alphas_pitch>=PITCH_BINS(pb)).*(alphas_pitch<PITCH_BINS(pb+1)));
%         poloidal_hist_ini(eb,pb,:,:) = hist2d_corr(pos_X_gc(PART_POP)+R0, pos_Z_gc(PART_POP), R_BINS, Z_BINS,alphas_weight(PART_POP));
%         data_array=density_correction_factor*squeeze(poloidal_hist_ini(eb,pb,:,:))./ volume_RZ;
%         dist_EpitchRZ_ini(eb,pb,:,:) =data_array;
%         data_array=data_array';
%         if SAVE_INI==1
%             save(FILENAME_INI,'-append','data_array','-ASCII')
%         end
%     end
% end
% 

clear Ekin_profile_ini rhist_ini

for rb=1:length(RADIAL_values)
    PART_POP=((r_avg>=RADIAL_BINS(rb)).*(r_avg<RADIAL_BINS(rb+1)));
    TR_BIN_POP=find(PART_POP.*ALL_TRAPPED_POP);
    PASS_BIN_POP=find(PART_POP.*ALL_PASSING_POP);
    trapped_fraction_profile_ini(rb)=length(TR_BIN_POP)/length(find(PART_POP));
%     lambda_profile_ini(rb)=mean(abs(alphas_lambda(find(PART_POP.*(alphas_lambda>=0)))));
%     lambda_profile_ini(rb)=mean(alphas_vpll(find(PART_POP.*(alphas_lambda>=0)))./alphas_vtot(find(PART_POP.*(alphas_lambda>=0))));
%     Ekin_profile_ini(rb)=mean(alphas_Ekin(find(PART_POP)));
    Ekin_profile_ini(rb)=mean(alphas_weight(find(PART_POP)).*alphas_Ekin(find(PART_POP)));
    rhist_ini(rb)=sum(alphas_weight(find(PART_POP)));
end



data_hist_radial_ini=zeros(length(RADIAL_values),length(Ekin_values),length(pitch_values));
data_n_radial_ini=zeros(length(RADIAL_values),length(Ekin_values),length(pitch_values));

TOMO_POS_PART=(alphas_pos_x>=RMIN_TOMO).*(alphas_pos_x<=RMAX_TOMO).*(alphas_pos_z>=ZMIN_TOMO).*(alphas_pos_z<=ZMAX_TOMO);


for rb=1:length(RADIAL_values)
    for eb=1:length(Ekin_values)
        for pb=1:length(pitch_values)
            PART_POP=find((r_avg>=RADIAL_BINS(rb)).*(r_avg<RADIAL_BINS(rb+1)).*(alphas_Ekin>=EKIN_BINS(eb)).*(alphas_Ekin<EKIN_BINS(eb+1))...
                .*(alphas_pitch>=PITCH_BINS(pb)).*(alphas_pitch<PITCH_BINS(pb+1))...
                .*TOMO_POS_PART);
            data_hist_radial_ini(rb,eb,pb)=sum(alphas_weight(PART_POP));
            data_n_radial_ini(rb,eb,pb)=data_hist_radial_ini(rb,eb,pb)/volume_radial_values(rb);
        end
    end
end

%%
% pcolor (r_values, pitch_values, mHist);
% shading faceted; 


POST_COLLAPSE_FILENAME
load(POST_COLLAPSE_FILENAME);
Nalphas_simulated=length(alphas_Ekin)
%load('final_NBI60keV_precession_stats_all.mat');
r_avg=interp2(scale_X,scale_Z,radial_XZsmall_map',pos_X_gc,pos_Z_gc);
r_avg=interp1(1:Nradial,radial_r_value_flux,r_avg);
% rhist_end=histc(r_avg,RADIAL_BINS);

alphas_B_end=interp2(scale_X,scale_Z,Btot_XZ_map',pos_X_gc,pos_Z_gc,'*linear');
alphas_theta_end=interp2(scale_X,scale_Z,theta_XZsmall_map',pos_X_gc,pos_Z_gc,'*linear');

alphas_Eperp=max(alphas_Ekin-0.5*(mHe/eV)*alphas_vpll.^2,0);
alphas_vperp=sqrt(2*(eV/mHe)*alphas_Eperp);
alphas_vtot=sqrt(2*(eV/mHe)*alphas_Ekin);
alphas_pitch=-alphas_vpll./alphas_vtot;
alphas_lambda_end=atan(alphas_vperp./alphas_vpll);
alphas_r=interp1(1:257,radial_r_value_flux,alphas_psi);
alphas_lambda0_end=Bavg*alphas_mm./alphas_Ekin;
alphas_Ekin_end=alphas_Ekin;

poloidal_hist_end=zeros(length(Ekin_values),length(pitch_values),length(R_values),length(Z_values));
dist_EpitchRZ_end=poloidal_hist_end;

% string_title='# final NBI distribution : dimensions ; scales (4 lines) ; [Ekin pitch R Z] '
% %save dist_EpitchRZ_end.dat -append string_title -ASCII
% if SAVE_FINAL==1
%     save(FILENAME_END,'-append','Raxis','Z_axis','-ASCII')
%     save(FILENAME_END,'-append','l1','l2','l3','l4','-ASCII')
%     save(FILENAME_END,'-append','Ekin_values','-ASCII')
%     save(FILENAME_END,'-append','pitch_values','-ASCII')
%     save(FILENAME_END,'-append','R_values','-ASCII')
%     save(FILENAME_END,'-append','Z_values','-ASCII')
% end
% 
% for eb=1:length(Ekin_values)
%     for pb=1:length(pitch_values)
%         PART_POP=find((~alphas_ejected).*(alphas_Ekin>=EKIN_BINS(eb)).*(alphas_Ekin<EKIN_BINS(eb+1)).*(alphas_pitch>=PITCH_BINS(pb)).*(alphas_pitch<PITCH_BINS(pb+1)));
% %         poloidal_hist_end(eb,pb,:,:) = hist2d ([pos_X_gc(PART_POP)+R0 pos_Z_gc(PART_POP)], R_BINS, Z_BINS);
%         poloidal_hist_end(eb,pb,:,:) = hist2d_corr(pos_X_gc(PART_POP)+R0, pos_Z_gc(PART_POP), R_BINS, Z_BINS,alphas_weight(PART_POP));
%         data_array=density_correction_factor*squeeze(poloidal_hist_end(eb,pb,:,:))./ volume_RZ;
%         dist_EpitchRZ_end(eb,pb,:,:) =data_array;  % cm3
%         data_array=data_array';
%         if SAVE_FINAL==1
%             save(FILENAME_END,'-append','data_array','-ASCII')
%         end
%     end
% end

clear Ekin_profile_end rhist_end

for rb=1:length(RADIAL_values)
    PART_POP=((~alphas_ejected).*(r_avg>=RADIAL_BINS(rb)).*(r_avg<RADIAL_BINS(rb+1)));
    TR_BIN_POP=find(PART_POP.*ALL_TRAPPED_POP(1:Nalphas_simulated));
    PASS_BIN_POP=find(PART_POP.*ALL_PASSING_POP(1:Nalphas_simulated));
    trapped_fraction_profile_end(rb)=length(TR_BIN_POP)/length(find(PART_POP));
%     lambda_profile_end(rb)=mean(abs(alphas_lambda(find(PART_POP.*(alphas_lambda>=0)))));
%     lambda_profile_end(rb)=mean(alphas_vpll(find(PART_POP.*(alphas_lambda>=0)))./alphas_vtot(find(PART_POP.*(alphas_lambda>=0))));
    Ekin_profile_end(rb)=mean(alphas_weight(find(PART_POP)).*alphas_Ekin(find(PART_POP)));
    rhist_end(rb)=sum(alphas_weight(find(PART_POP)));
end


data_hist_radial_end=zeros(length(RADIAL_values),length(Ekin_values),length(pitch_values));
data_n_radial_end=zeros(length(RADIAL_values),length(Ekin_values),length(pitch_values));

%%
TOMO_POS_PART=(alphas_pos_x>=RMIN_TOMO).*(alphas_pos_x<=RMAX_TOMO).*(alphas_pos_z>=ZMIN_TOMO).*(alphas_pos_z<=ZMAX_TOMO);

for rb=1:length(RADIAL_values)
    for eb=1:length(Ekin_values)
        for pb=1:length(pitch_values)
            PART_POP=find((~alphas_ejected).*(r_avg>=RADIAL_BINS(rb)).*(r_avg<RADIAL_BINS(rb+1)).*(alphas_Ekin>=EKIN_BINS(eb)).*(alphas_Ekin<EKIN_BINS(eb+1))...
                .*(alphas_pitch>=PITCH_BINS(pb)).*(alphas_pitch<PITCH_BINS(pb+1))...
                .*TOMO_POS_PART);
            data_hist_radial_end(rb,eb,pb)=sum(alphas_weight(PART_POP));
            data_n_radial_end(rb,eb,pb)=data_hist_radial_end(rb,eb,pb)/volume_radial_values(rb);           
        end
    end
end
% save large_crash_NBI_distributions.mat dist_EpitchRZ_ini dist_EpitchRZ_end Ekin_values pitch_values R_values Z_values trapped_fraction_profile_ini trapped_fraction_profile_end RADIAL_values

% %%
% 
% figure(1)
% set(gca,'fontsize',20)
% contourf(R_values,Z_values,(squeeze(sum(sum(poloidal_hist_end(1:end,:,:,:),2),1))'),12);axis xy;hold on
% contour(scale_X+R0,scale_Z,psi_norm_XZsmall_map','k','linewidth',2)
% 
% %%
% figure(2)
% set(gca,'fontsize',20)
% contourf(R_values,Z_values,(squeeze(sum(sum(poloidal_hist_ini(1:end,:,:,:),2),1))'),12);axis xy;hold on
% contour(scale_X+R0,scale_Z,psi_norm_XZsmall_map','k','linewidth',2)



%%
figure(5)
set(gca,'fontsize',20)

hold on
grid on

rb_inf=2
rb_sup=2

EKINBIN1=2
EKINBIN2=2
E1_value=0.5*(Ekin_values(EKINBIN1)+Ekin_values(EKINBIN2))
E1_value=E1_value*1e-3

density_NBI_given_EKin_ini=squeeze(sum(sum(data_n_radial_ini(rb_inf:rb_sup,EKINBIN1:EKINBIN2,:),1),2));
density_NBI_given_EKin_end=squeeze(sum(sum(data_n_radial_end(rb_inf:rb_sup,EKINBIN1:EKINBIN2,:),1),2));

plot(pitch_values,-(density_NBI_given_EKin_end-density_NBI_given_EKin_ini)./density_NBI_given_EKin_ini,'b--','linewidth',2)


EKINBIN1=3
EKINBIN2=3
E2_value=0.5*(Ekin_values(EKINBIN1)+Ekin_values(EKINBIN2))
E2_value=E2_value*1e-3

density_NBI_given_EKin_ini=squeeze(sum(sum(data_n_radial_ini(rb_inf:rb_sup,EKINBIN1:EKINBIN2,:),1),2));
density_NBI_given_EKin_end=squeeze(sum(sum(data_n_radial_end(rb_inf:rb_sup,EKINBIN1:EKINBIN2,:),1),2));


plot(pitch_values,-(density_NBI_given_EKin_end-density_NBI_given_EKin_ini)./density_NBI_given_EKin_ini,'r--','linewidth',2)




EKINBIN1=5
EKINBIN2=5
E3_value=0.5*(Ekin_values(EKINBIN1)+Ekin_values(EKINBIN2))
E3_value=E3_value*1e-3

density_NBI_given_EKin_ini=squeeze(sum(sum(data_n_radial_ini(rb_inf:rb_sup,EKINBIN1:EKINBIN2,:),1),2));
density_NBI_given_EKin_end=squeeze(sum(sum(data_n_radial_end(rb_inf:rb_sup,EKINBIN1:EKINBIN2,:),1),2));


plot(pitch_values,-(density_NBI_given_EKin_end-density_NBI_given_EKin_ini)./density_NBI_given_EKin_ini,'g--','linewidth',2)




xlim([-0.6 1])

% legend('24 keV','36 keV','51 keV')

xlabel('v_{||} / v_{tot}')
ylabel('n_{NBI} (m^{-3}) on [16-40] keV')
yl=ylabel('$-\Delta n / n_{ini}$')
set(yl,'Interpreter','latex')




Erank1=2;
E1=E_tomo_values(Erank1)
% plot(sign(p_values).*(abs(p_values)),(BEFORE_OFFSET*tomo_before(:,Erank1)-tomo_after(:,Erank1))./tomo_before(:,Erank1),'b','linewidth',3);
errorbar(sign(p_values).*(abs(p_values)),(BEFORE_OFFSET*tomo_before(:,Erank1)-tomo_after(:,Erank1))./tomo_before(:,Erank1),tomo_rel_error(:,Erank1),'b','linewidth',3);

% Erank2=5;
% E2=E_tomo_values(Erank2)
% plot(sign(p_values).*(abs(p_values)),(BEFORE_OFFSET*tomo_before(:,Erank2)-tomo_after(:,Erank2))./tomo_before(:,Erank2),'k--','linewidth',3);

Erank3=4;
E3=E_tomo_values(Erank3)
% plot(sign(p_values(3:end)).*(abs(p_values(3:end))),(BEFORE_OFFSET*tomo_before(3:end,Erank3)-tomo_after(3:end,Erank3))./tomo_before(3:end,Erank3),'r','linewidth',3);
errorbar(sign(p_values(3:end)).*(abs(p_values(3:end))),(BEFORE_OFFSET*tomo_before(3:end,Erank3)-tomo_after(3:end,Erank3))./tomo_before(3:end,Erank3),tomo_rel_error(3:end,Erank3),'r','linewidth',3);

Erank4=7;
E4=E_tomo_values(Erank4)
% plot(sign(p_values(3:end)).*(abs(p_values(3:end))),(BEFORE_OFFSET*tomo_before(3:end,Erank4)-tomo_after(3:end,Erank4))./tomo_before(3:end,Erank4),'g','linewidth',3);
errorbar(sign(p_values(3:end)).*(abs(p_values(3:end))),(BEFORE_OFFSET*tomo_before(3:end,Erank4)-tomo_after(3:end,Erank4))./tomo_before(3:end,Erank4),tomo_rel_error(3:end,Erank4),'g','linewidth',3);


% ylabel('-\Delta n / n_{ini}')
% hl=legend(strcat('EBdyna :',num2str(E1_value),' keV'),strcat('EBdyna :',num2str(E2_value),' keV'),strcat('EBdyna :',num2str(E3_value),' keV'))
% set(hl,'fontsize',16)

hl=legend(strcat('EBdyna :',num2str(E1_value),' keV'),strcat('EBdyna :',num2str(floor(E2_value)),' keV'),strcat('EBdyna :',num2str(E3_value),' keV'), ...
    strcat('tomo :',num2str(round(E1)),' keV'),strcat('tomo :',num2str(round(E3)),' keV'),strcat('tomo :',num2str(round(E4)),' keV'))
set(hl,'fontsize',16)


xlim([-0.6 0.9])
ylim([-0.2 0.45])


%%
figure(4)
hold on
grid on
rho_lin_scale=(0:256)/256;

plot(rho_lin_scale,PFAST_ini,'b','linewidth',2)
plot(rho_lin_scale,PFAST_end,'r','linewidth',2)
xlabel('\rho _{tor}')
ylabel('Ptot (Pa)')

hold on
rho_tor_scale=sqrt(tor_flux_profile/max(tor_flux_profile));
rho_tor_scale_values=interp1(radial_r_value_flux,rho_tor_scale,RADIAL_values);

density_NBI_ini=rhist_ini(1:end)./volume_radial_values;
density_NBI_end=rhist_end(1:end)./volume_radial_values;
PNBI_ini=(2/3)*eV*NB_PART_SIM_TRANSP_RATIO*weight_transp_NBI*density_NBI_ini.*Ekin_profile_ini;
PNBI_end=(2/3)*eV*NB_PART_SIM_TRANSP_RATIO*weight_transp_NBI*density_NBI_end.*Ekin_profile_end;


PNBI_ini_rholin=interp1(rho_tor_scale_values,PNBI_ini,(0:257-1)/(257-1));
r_rholin=interp1(rho_tor_scale,radial_r_value_flux,rho_lin_scale);

weight_old_rholin=interp1(radial_r_value_flux,alphas_weight_radial_profile_old,r_rholin);

PNBI_norm=1.0;
PNBI_ini=PNBI_ini/PNBI_norm;
PNBI_end=PNBI_end/PNBI_norm;





hold on
grid on
plot(rho_tor_scale_values,PNBI_ini,'b--')
plot(rho_tor_scale_values,PNBI_end,'r--')
xlabel('\rho_{tor}')
% ylabel('P/P_{tot}')
ylabel('P_{NBI} (Pa)')




