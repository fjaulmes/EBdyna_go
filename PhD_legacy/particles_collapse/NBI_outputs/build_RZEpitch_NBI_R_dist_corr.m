reset_data_analysis_environment
close all

weight_transp_NBI=1.4e+13
%  31557         # shotnumber
%  2.24330       # time
%  1000014       # fast-ion markers
%   1.40089e+13  # weight [fast ions per marker]

load('initial_NBI60keV_R_precession_stats_all.mat');
% NB_PART_SIM_TRANSP_RATIO=202016/length(r_avg)
% NB_PART_SIM_TRANSP_RATIO=0.5*(0.25+499988/length(r_avg))
NB_PART_SIM_TRANSP_RATIO=(1000014/length(r_avg))

NB_PART_SIM_TRANSP_RATIO=1.05*NB_PART_SIM_TRANSP_RATIO



% density_correction=(1e-6)*NB_PART_SIM_TRANSP_RATIO*weight_transp_NBI

%%
load('initial_NBI60keV_Rlab_pre_collapse_all.mat');

load('weight_correction_profileP.mat')
alphas_weight_radial_profile_old=alphas_weight_radial_profile;
alphas_weight_corr=interp1(radial_r_value_flux,alphas_weight_radial_profile,r_avg);
alphas_weight=alphas_weight_corr;

%%

%EKIN_INF=1*1e3

% PART_POP=find((alphas_Ekin>EKIN_INF).*(alphas_pos_z<0.2).*(alphas_pos_z>-0.2));
%PART_POP=find((alphas_Ekin>EKIN_INF));

EKIN_BIN_SIZE=6*1e3;
EKIN_BINS=(0:EKIN_BIN_SIZE:84*1e3);
Ekin_values=EKIN_BINS(1:end-1)+0.5*EKIN_BIN_SIZE;

PICH_BIN_SIZE=0.08;
PITCH_BINS=(-1.0:PICH_BIN_SIZE:1.1);
pitch_values=PITCH_BINS(1:end-1)+0.5*PICH_BIN_SIZE;

R_BIN_SIZE=0.05;
R_BINS=(-0.6:R_BIN_SIZE:0.6)+R0+X_axis;
R_values=R_BINS(1:end-1)+0.5*R_BIN_SIZE;

Z_BIN_SIZE=0.05;
Z_BINS=(-0.8:Z_BIN_SIZE:0.8);
Z_values=Z_BINS(1:end-1)+0.5*Z_BIN_SIZE;

RADIAL_BIN_SIZE=0.03;
RADIAL_BINS=(0:RADIAL_BIN_SIZE:0.6);
RADIAL_values=RADIAL_BINS(1:end-1)+0.5*RADIAL_BIN_SIZE;
psi_radial_values=interp1(radial_r_value_flux,1:Nradial,RADIAL_values);
volume_radial_bins_values=interp1(radial_r_value_flux,volume_flux,RADIAL_BINS);
volume_radial_values=volume_radial_bins_values(2:end)-volume_radial_bins_values(1:end-1);

%%

rhist_ini=histc(r_avg,RADIAL_BINS);


l1=length(Ekin_values) 
l2=length(pitch_values)
l3=length(R_values) 
l4=length(Z_values)
l5=length(RADIAL_values)

trapped_fraction_profile_ini=zeros(l5,1);
trapped_fraction_profile_end=zeros(l5,1);
lambda_profile_ini=zeros(l5,1);
lambda_profile_end=zeros(l5,1);


density_correction_factor=(1e-6)*NB_PART_SIM_TRANSP_RATIO*weight_transp_NBI/(EKIN_BIN_SIZE)/PICH_BIN_SIZE

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


FILENAME_INI='dist_NBI_Pcorr_AUG31557_2p25_EpitchRZ_ini.dat'
FILENAME_END='dist_NBI_Pcorr_AUG31557_2p25_EpitchRZ_125_end.dat'

string_title='# initial NBI distribution : dimensions ; scales (4 lines) ; [Ekin pitch R Z] '

Raxis=R0+X_axis

%save dist_EpitchRZ_ini.dat -append  string_title -ASCII
if SAVE_INI==1
    save(FILENAME_INI,'-append','Raxis','Z_axis','-ASCII')
    save(FILENAME_INI,'-append','l1','l2','l3','l4','-ASCII')
    save(FILENAME_INI,'-append','Ekin_values','-ASCII')
    save(FILENAME_INI,'-append','pitch_values','-ASCII')
    save(FILENAME_INI,'-append','R_values','-ASCII')
    save(FILENAME_INI,'-append','Z_values','-ASCII')
end

for eb=1:length(Ekin_values)
    for pb=1:length(pitch_values)
        PART_POP=find((alphas_Ekin>=EKIN_BINS(eb)).*(alphas_Ekin<EKIN_BINS(eb+1)).*(alphas_pitch>=PITCH_BINS(pb)).*(alphas_pitch<PITCH_BINS(pb+1)));
        poloidal_hist_ini(eb,pb,:,:) = hist2d_corr(pos_X_gc(PART_POP)+R0, pos_Z_gc(PART_POP), R_BINS, Z_BINS,alphas_weight(PART_POP));
        data_array=density_correction_factor*squeeze(poloidal_hist_ini(eb,pb,:,:))./ volume_RZ;
        dist_EpitchRZ_ini(eb,pb,:,:) =data_array;
        data_array=data_array';
        if SAVE_INI==1
            save(FILENAME_INI,'-append','data_array','-ASCII')
        end
    end
end



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

%%
% pcolor (r_values, pitch_values, mHist);
% shading faceted; 



load('NBI60keV_Rlab_fc1p6h1p6_all.mat');
Nalphas_simulated=length(alphas_Ekin)
%load('final_NBI60keV_precession_stats_all.mat');
r_avg=interp2(scale_X,scale_Z,radial_XZsmall_map',pos_X_gc,pos_Z_gc);
r_avg=interp1(1:Nradial,radial_r_value_flux,r_avg);
rhist_end=histc(r_avg,RADIAL_BINS);

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

string_title='# final NBI distribution : dimensions ; scales (4 lines) ; [Ekin pitch R Z] '
%save dist_EpitchRZ_end.dat -append string_title -ASCII
if SAVE_FINAL==1
    save(FILENAME_END,'-append','Raxis','Z_axis','-ASCII')
    save(FILENAME_END,'-append','l1','l2','l3','l4','-ASCII')
    save(FILENAME_END,'-append','Ekin_values','-ASCII')
    save(FILENAME_END,'-append','pitch_values','-ASCII')
    save(FILENAME_END,'-append','R_values','-ASCII')
    save(FILENAME_END,'-append','Z_values','-ASCII')
end

for eb=1:length(Ekin_values)
    for pb=1:length(pitch_values)
        PART_POP=find((~alphas_ejected).*(alphas_Ekin>=EKIN_BINS(eb)).*(alphas_Ekin<EKIN_BINS(eb+1)).*(alphas_pitch>=PITCH_BINS(pb)).*(alphas_pitch<PITCH_BINS(pb+1)));
%         poloidal_hist_end(eb,pb,:,:) = hist2d ([pos_X_gc(PART_POP)+R0 pos_Z_gc(PART_POP)], R_BINS, Z_BINS);
        poloidal_hist_end(eb,pb,:,:) = hist2d_corr(pos_X_gc(PART_POP)+R0, pos_Z_gc(PART_POP), R_BINS, Z_BINS,alphas_weight(PART_POP));
        data_array=density_correction_factor*squeeze(poloidal_hist_end(eb,pb,:,:))./ volume_RZ;
        dist_EpitchRZ_end(eb,pb,:,:) =data_array;  % cm3
        data_array=data_array';
        if SAVE_FINAL==1
            save(FILENAME_END,'-append','data_array','-ASCII')
        end
    end
end

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

% save large_crash_NBI_distributions.mat dist_EpitchRZ_ini dist_EpitchRZ_end Ekin_values pitch_values R_values Z_values trapped_fraction_profile_ini trapped_fraction_profile_end RADIAL_values

%%
dist_EpitchRZ_end_value=dist_EpitchRZ_end*(1e6)*EKIN_BIN_SIZE*PICH_BIN_SIZE;
figure(1)
set(gca,'fontsize',20)
contourf(R_values,Z_values,(squeeze(sum(sum(dist_EpitchRZ_end_value(1:end,:,:,:),2),1))'),12);axis xy;hold on
colorbar
contour(scale_X+R0,scale_Z,psi_norm_XZsmall_map','k','linewidth',2)

%%
dist_EpitchRZ_ini_value=dist_EpitchRZ_ini*(1e6)*EKIN_BIN_SIZE*PICH_BIN_SIZE;
figure(2)
set(gca,'fontsize',20)
contourf(R_values,Z_values,(squeeze(sum(sum(dist_EpitchRZ_ini_value(1:end,:,:,:),2),1))'),12);axis xy;hold on
colorbar
contour(scale_X+R0,scale_Z,psi_norm_XZsmall_map','k','linewidth',2)

%%
figure(3)
hold on
grid on
plot(R_values,(1e6)*EKIN_BIN_SIZE*PICH_BIN_SIZE*(squeeze(sum(sum(sum(dist_EpitchRZ_ini(1:end,:,:,floor(0.5*(length(Z_values))-2):ceil(0.5*(length(Z_values)))+2),2),1),4)/5)'),'r')
plot(R_values,(1e6)*EKIN_BIN_SIZE*PICH_BIN_SIZE*(squeeze(sum(sum(sum(dist_EpitchRZ_end(1:end,:,:,floor(0.5*(length(Z_values))-2):ceil(0.5*(length(Z_values)))+2),2),1),4)/5)'),'b')


%%
figure(5)
set(gca,'fontsize',20)

hold on
grid on
DELTA_SQUARE=2

EKINBIN1=4
EKINBIN2=5

midRinf=floor(0.5*(length(R_values)));
midZinf=floor(0.5*(length(Z_values)))+1;
midRsup=ceil(0.5*(length(R_values)));
midZsup=ceil(0.5*(length(Z_values)))+1;

DbinR=midRsup-midRinf
DbinZ=midZsup-midZinf
% DbinR=1
% DbinZ=1

Nbsurfaxis=(1+2*DELTA_SQUARE+DbinR)*(1+2*DELTA_SQUARE+DbinZ)+2*(1+DbinR+DELTA_SQUARE-1)+2*(1+DbinZ+DELTA_SQUARE-1)

density_NBI_given_EKin_ini=(1e6)*EKIN_BIN_SIZE*PICH_BIN_SIZE*(squeeze(sum(sum(sum(dist_EpitchRZ_ini(EKINBIN1:EKINBIN2,:,midRinf-DELTA_SQUARE:midRsup+DELTA_SQUARE,midZinf-DELTA_SQUARE:midZsup+DELTA_SQUARE),1),3),4))');
density_NBI_given_EKin_end=(1e6)*EKIN_BIN_SIZE*PICH_BIN_SIZE*(squeeze(sum(sum(sum(dist_EpitchRZ_end(EKINBIN1:EKINBIN2,:,midRinf-DELTA_SQUARE:midRsup+DELTA_SQUARE,midZinf-DELTA_SQUARE:midZsup+DELTA_SQUARE),1),3),4))');

density_NBI_given_EKin_ini=density_NBI_given_EKin_ini+(1e6)*EKIN_BIN_SIZE*PICH_BIN_SIZE*squeeze(...
    sum(sum(sum(dist_EpitchRZ_ini(EKINBIN1:EKINBIN2,:,midRinf-DELTA_SQUARE+1:midRsup+DELTA_SQUARE-1,midZinf-DELTA_SQUARE-1:midZsup-DELTA_SQUARE-1),1),3),4))';
density_NBI_given_EKin_ini=density_NBI_given_EKin_ini+(1e6)*EKIN_BIN_SIZE*PICH_BIN_SIZE*squeeze(...
    sum(sum(sum(dist_EpitchRZ_ini(EKINBIN1:EKINBIN2,:,midRinf-DELTA_SQUARE+1:midRsup+DELTA_SQUARE-1,midZinf+DELTA_SQUARE+1:midZsup+DELTA_SQUARE+1),1),3),4))';
density_NBI_given_EKin_ini=density_NBI_given_EKin_ini+(1e6)*EKIN_BIN_SIZE*PICH_BIN_SIZE*squeeze(...
    sum(sum(sum(dist_EpitchRZ_ini(EKINBIN1:EKINBIN2,:,midRinf-DELTA_SQUARE-1:midRsup-DELTA_SQUARE-1,midZinf-DELTA_SQUARE+1:midZsup+DELTA_SQUARE-1),1),3),4))';
density_NBI_given_EKin_ini=density_NBI_given_EKin_ini+(1e6)*EKIN_BIN_SIZE*PICH_BIN_SIZE*squeeze(...
    sum(sum(sum(dist_EpitchRZ_ini(EKINBIN1:EKINBIN2,:,midRinf+DELTA_SQUARE+1:midRsup+DELTA_SQUARE+1,midZinf-DELTA_SQUARE+1:midZsup+DELTA_SQUARE-1),1),3),4))';

density_NBI_given_EKin_ini=density_NBI_given_EKin_ini/Nbsurfaxis;

density_NBI_given_EKin_end=density_NBI_given_EKin_end+(1e6)*EKIN_BIN_SIZE*PICH_BIN_SIZE*squeeze(...
    sum(sum(sum(dist_EpitchRZ_ini(EKINBIN1:EKINBIN2,:,midRinf-DELTA_SQUARE+1:midRsup+DELTA_SQUARE-1,midZinf-DELTA_SQUARE-1:midZsup-DELTA_SQUARE-1),1),3),4))';
density_NBI_given_EKin_end=density_NBI_given_EKin_end+(1e6)*EKIN_BIN_SIZE*PICH_BIN_SIZE*squeeze(...
    sum(sum(sum(dist_EpitchRZ_ini(EKINBIN1:EKINBIN2,:,midRinf-DELTA_SQUARE+1:midRsup+DELTA_SQUARE-1,midZinf+DELTA_SQUARE+1:midZsup+DELTA_SQUARE+1),1),3),4))';
density_NBI_given_EKin_end=density_NBI_given_EKin_end+(1e6)*EKIN_BIN_SIZE*PICH_BIN_SIZE*squeeze(...
    sum(sum(sum(dist_EpitchRZ_ini(EKINBIN1:EKINBIN2,:,midRinf-DELTA_SQUARE-1:midRsup-DELTA_SQUARE-1,midZinf-DELTA_SQUARE+1:midZsup+DELTA_SQUARE-1),1),3),4))';
density_NBI_given_EKin_end=density_NBI_given_EKin_end+(1e6)*EKIN_BIN_SIZE*PICH_BIN_SIZE*squeeze(...
    sum(sum(sum(dist_EpitchRZ_ini(EKINBIN1:EKINBIN2,:,midRinf+DELTA_SQUARE+1:midRsup+DELTA_SQUARE+1,midZinf-DELTA_SQUARE+1:midZsup+DELTA_SQUARE-1),1),3),4))';

density_NBI_given_EKin_end=density_NBI_given_EKin_end/Nbsurfaxis;

% plot(pitch_values,density_NBI_given_EKin_ini,'r--','linewidth',2)
% plot(pitch_values,density_NBI_given_EKin_end,'b','linewidth',2)
plot(pitch_values,-(density_NBI_given_EKin_end-density_NBI_given_EKin_ini)./density_NBI_given_EKin_ini,'r--','linewidth',2)


EKINBIN1=6
EKINBIN2=7

density_NBI_given_EKin_ini=(1e6)*EKIN_BIN_SIZE*PICH_BIN_SIZE*(squeeze(sum(sum(sum(dist_EpitchRZ_ini(EKINBIN1:EKINBIN2,:,midRinf-DELTA_SQUARE:midRsup+DELTA_SQUARE,midZinf-DELTA_SQUARE:midZsup+DELTA_SQUARE),1),3),4))');
density_NBI_given_EKin_end=(1e6)*EKIN_BIN_SIZE*PICH_BIN_SIZE*(squeeze(sum(sum(sum(dist_EpitchRZ_end(EKINBIN1:EKINBIN2,:,midRinf-DELTA_SQUARE:midRsup+DELTA_SQUARE,midZinf-DELTA_SQUARE:midZsup+DELTA_SQUARE),1),3),4))');

density_NBI_given_EKin_ini=density_NBI_given_EKin_ini+(1e6)*EKIN_BIN_SIZE*PICH_BIN_SIZE*squeeze(...
    sum(sum(sum(dist_EpitchRZ_ini(EKINBIN1:EKINBIN2,:,midRinf-DELTA_SQUARE+1:midRsup+DELTA_SQUARE-1,midZinf-DELTA_SQUARE-1:midZsup-DELTA_SQUARE-1),1),3),4))';
density_NBI_given_EKin_ini=density_NBI_given_EKin_ini+(1e6)*EKIN_BIN_SIZE*PICH_BIN_SIZE*squeeze(...
    sum(sum(sum(dist_EpitchRZ_ini(EKINBIN1:EKINBIN2,:,midRinf-DELTA_SQUARE+1:midRsup+DELTA_SQUARE-1,midZinf+DELTA_SQUARE+1:midZsup+DELTA_SQUARE+1),1),3),4))';
density_NBI_given_EKin_ini=density_NBI_given_EKin_ini+(1e6)*EKIN_BIN_SIZE*PICH_BIN_SIZE*squeeze(...
    sum(sum(sum(dist_EpitchRZ_ini(EKINBIN1:EKINBIN2,:,midRinf-DELTA_SQUARE-1:midRsup-DELTA_SQUARE-1,midZinf-DELTA_SQUARE+1:midZsup+DELTA_SQUARE-1),1),3),4))';
density_NBI_given_EKin_ini=density_NBI_given_EKin_ini+(1e6)*EKIN_BIN_SIZE*PICH_BIN_SIZE*squeeze(...
    sum(sum(sum(dist_EpitchRZ_ini(EKINBIN1:EKINBIN2,:,midRinf+DELTA_SQUARE+1:midRsup+DELTA_SQUARE+1,midZinf-DELTA_SQUARE+1:midZsup+DELTA_SQUARE-1),1),3),4))';

density_NBI_given_EKin_ini=density_NBI_given_EKin_ini/Nbsurfaxis;

density_NBI_given_EKin_end=density_NBI_given_EKin_end+(1e6)*EKIN_BIN_SIZE*PICH_BIN_SIZE*squeeze(...
    sum(sum(sum(dist_EpitchRZ_ini(EKINBIN1:EKINBIN2,:,midRinf-DELTA_SQUARE+1:midRsup+DELTA_SQUARE-1,midZinf-DELTA_SQUARE-1:midZsup-DELTA_SQUARE-1),1),3),4))';
density_NBI_given_EKin_end=density_NBI_given_EKin_end+(1e6)*EKIN_BIN_SIZE*PICH_BIN_SIZE*squeeze(...
    sum(sum(sum(dist_EpitchRZ_ini(EKINBIN1:EKINBIN2,:,midRinf-DELTA_SQUARE+1:midRsup+DELTA_SQUARE-1,midZinf+DELTA_SQUARE+1:midZsup+DELTA_SQUARE+1),1),3),4))';
density_NBI_given_EKin_end=density_NBI_given_EKin_end+(1e6)*EKIN_BIN_SIZE*PICH_BIN_SIZE*squeeze(...
    sum(sum(sum(dist_EpitchRZ_ini(EKINBIN1:EKINBIN2,:,midRinf-DELTA_SQUARE-1:midRsup-DELTA_SQUARE-1,midZinf-DELTA_SQUARE+1:midZsup+DELTA_SQUARE-1),1),3),4))';
density_NBI_given_EKin_end=density_NBI_given_EKin_end+(1e6)*EKIN_BIN_SIZE*PICH_BIN_SIZE*squeeze(...
    sum(sum(sum(dist_EpitchRZ_ini(EKINBIN1:EKINBIN2,:,midRinf+DELTA_SQUARE+1:midRsup+DELTA_SQUARE+1,midZinf-DELTA_SQUARE+1:midZsup+DELTA_SQUARE-1),1),3),4))';

density_NBI_given_EKin_end=density_NBI_given_EKin_end/Nbsurfaxis;

plot(pitch_values,-(density_NBI_given_EKin_end-density_NBI_given_EKin_ini)./density_NBI_given_EKin_ini,'b--','linewidth',2)




EKINBIN1=8
EKINBIN2=9

density_NBI_given_EKin_ini=(1e6)*EKIN_BIN_SIZE*PICH_BIN_SIZE*(squeeze(sum(sum(sum(dist_EpitchRZ_ini(EKINBIN1:EKINBIN2,:,midRinf-DELTA_SQUARE:midRsup+DELTA_SQUARE,midZinf-DELTA_SQUARE:midZsup+DELTA_SQUARE),1),3),4))');
density_NBI_given_EKin_end=(1e6)*EKIN_BIN_SIZE*PICH_BIN_SIZE*(squeeze(sum(sum(sum(dist_EpitchRZ_end(EKINBIN1:EKINBIN2,:,midRinf-DELTA_SQUARE:midRsup+DELTA_SQUARE,midZinf-DELTA_SQUARE:midZsup+DELTA_SQUARE),1),3),4))');

density_NBI_given_EKin_ini=density_NBI_given_EKin_ini+(1e6)*EKIN_BIN_SIZE*PICH_BIN_SIZE*squeeze(...
    sum(sum(sum(dist_EpitchRZ_ini(EKINBIN1:EKINBIN2,:,midRinf-DELTA_SQUARE+1:midRsup+DELTA_SQUARE-1,midZinf-DELTA_SQUARE-1:midZsup-DELTA_SQUARE-1),1),3),4))';
density_NBI_given_EKin_ini=density_NBI_given_EKin_ini+(1e6)*EKIN_BIN_SIZE*PICH_BIN_SIZE*squeeze(...
    sum(sum(sum(dist_EpitchRZ_ini(EKINBIN1:EKINBIN2,:,midRinf-DELTA_SQUARE+1:midRsup+DELTA_SQUARE-1,midZinf+DELTA_SQUARE+1:midZsup+DELTA_SQUARE+1),1),3),4))';
density_NBI_given_EKin_ini=density_NBI_given_EKin_ini+(1e6)*EKIN_BIN_SIZE*PICH_BIN_SIZE*squeeze(...
    sum(sum(sum(dist_EpitchRZ_ini(EKINBIN1:EKINBIN2,:,midRinf-DELTA_SQUARE-1:midRsup-DELTA_SQUARE-1,midZinf-DELTA_SQUARE+1:midZsup+DELTA_SQUARE-1),1),3),4))';
density_NBI_given_EKin_ini=density_NBI_given_EKin_ini+(1e6)*EKIN_BIN_SIZE*PICH_BIN_SIZE*squeeze(...
    sum(sum(sum(dist_EpitchRZ_ini(EKINBIN1:EKINBIN2,:,midRinf+DELTA_SQUARE+1:midRsup+DELTA_SQUARE+1,midZinf-DELTA_SQUARE+1:midZsup+DELTA_SQUARE-1),1),3),4))';

density_NBI_given_EKin_ini=density_NBI_given_EKin_ini/Nbsurfaxis;

density_NBI_given_EKin_end=density_NBI_given_EKin_end+(1e6)*EKIN_BIN_SIZE*PICH_BIN_SIZE*squeeze(...
    sum(sum(sum(dist_EpitchRZ_ini(EKINBIN1:EKINBIN2,:,midRinf-DELTA_SQUARE+1:midRsup+DELTA_SQUARE-1,midZinf-DELTA_SQUARE-1:midZsup-DELTA_SQUARE-1),1),3),4))';
density_NBI_given_EKin_end=density_NBI_given_EKin_end+(1e6)*EKIN_BIN_SIZE*PICH_BIN_SIZE*squeeze(...
    sum(sum(sum(dist_EpitchRZ_ini(EKINBIN1:EKINBIN2,:,midRinf-DELTA_SQUARE+1:midRsup+DELTA_SQUARE-1,midZinf+DELTA_SQUARE+1:midZsup+DELTA_SQUARE+1),1),3),4))';
density_NBI_given_EKin_end=density_NBI_given_EKin_end+(1e6)*EKIN_BIN_SIZE*PICH_BIN_SIZE*squeeze(...
    sum(sum(sum(dist_EpitchRZ_ini(EKINBIN1:EKINBIN2,:,midRinf-DELTA_SQUARE-1:midRsup-DELTA_SQUARE-1,midZinf-DELTA_SQUARE+1:midZsup+DELTA_SQUARE-1),1),3),4))';
density_NBI_given_EKin_end=density_NBI_given_EKin_end+(1e6)*EKIN_BIN_SIZE*PICH_BIN_SIZE*squeeze(...
    sum(sum(sum(dist_EpitchRZ_ini(EKINBIN1:EKINBIN2,:,midRinf+DELTA_SQUARE+1:midRsup+DELTA_SQUARE+1,midZinf-DELTA_SQUARE+1:midZsup+DELTA_SQUARE-1),1),3),4))';

density_NBI_given_EKin_end=density_NBI_given_EKin_end/Nbsurfaxis;

plot(pitch_values,-(density_NBI_given_EKin_end-density_NBI_given_EKin_ini)./density_NBI_given_EKin_ini,'g--','linewidth',2)




xlim([-0.6 1])

legend('24 keV','36 keV','48 keV')

xlabel('v_{||} / v_{tot}')
ylabel('n_{NBI} (m^{-3}) on [16-40] keV')
ylabel('\Delta n / n_{ini}')


%%
figure(4)
set(gca,'fontsize',22)
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

density_NBI_ini=rhist_ini(1:end-1)./volume_radial_values';
density_NBI_end=rhist_end(1:end-1)./volume_radial_values';
PNBI_ini=(2/3)*eV*NB_PART_SIM_TRANSP_RATIO*weight_transp_NBI*density_NBI_ini.*Ekin_profile_ini';
PNBI_end=(2/3)*eV*NB_PART_SIM_TRANSP_RATIO*weight_transp_NBI*density_NBI_end.*Ekin_profile_end';


PNBI_ini_rholin=interp1(rho_tor_scale_values,PNBI_ini,(0:257-1)/(257-1));
r_rholin=interp1(rho_tor_scale,radial_r_value_flux,rho_lin_scale);

weight_old_rholin=interp1(radial_r_value_flux,alphas_weight_radial_profile_old,r_rholin);

PNBI_norm=1.0;
PNBI_ini=PNBI_ini/PNBI_norm;
PNBI_end=PNBI_end/PNBI_norm;


% alphas_weight_rho_profile_new=PFAST_ini./PNBI_ini_rholin;



delta_profile=PFAST_ini-PNBI_ini_rholin;
alphas_weight_rho_profile_new=weight_old_rholin+0.7*delta_profile./PFAST_ini;
alphas_weight_rho_profile_new(isnan(alphas_weight_rho_profile_new))=max(alphas_weight_rho_profile_new);
alphas_weight_rho_profile_new(1:9)=alphas_weight_rho_profile_new(9);
alphas_weight_rho_profile_new=max(alphas_weight_rho_profile_new,0.9);
alphas_weight_radial_profile=interp1(rho_lin_scale,alphas_weight_rho_profile_new,rho_tor_scale);
alphas_weight_radial_profile=min(alphas_weight_radial_profile,2);

% alphas_weight_radial_profile=0.1+0.6*alphas_weight_radial_profile_new+0.3*alphas_weight_radial_profile;
% alphas_weight_radial_profile=min(alphas_weight_radial_profile,alphas_weight_radial_profile_new);

% adjust the outer part of the pressure
% alphas_weight_radial_profile(size_r:Nradial)=(alphas_weight_radial_profile(size_r:Nradial)+0.2*(alphas_weight_radial_profile_new(size_r:Nradial)-alphas_weight_radial_profile(size_r:Nradial)).*(1:Nradial-size_r+1)/(Nradial-size_r+1))

% use the old one when visulaization ok
% alphas_weight_radial_profile=alphas_weight_radial_profile_old;
% weight_poly=polyfit(rho_tor_scale,alphas_weight_radial_profile,4);
% alphas_weight_radial_profile=polyval(weight_poly,rho_tor_scale);
% weight_poly=polyfit(rho_tor_scale,alphas_weight_radial_profile,6);
% alphas_weight_radial_profile=polyval(weight_poly,rho_tor_scale);

% save weight_correction_profileP.mat radial_r_value_flux alphas_weight_radial_profile alphas_weight_radial_profile_old


hold on
grid on
plot(rho_tor_scale_values,PNBI_ini,'b--')
plot(rho_tor_scale_values,PNBI_end,'r--')
xlabel('\rho_{tor}')
% ylabel('P/P_{tot}')
ylabel('P_{NBI} (Pa)')

legend('TRANSP ini','TRANSP end','EBdyna ini','EBdyna end')

%%
figure(8)
hold on
grid on

plot(rho_tor_scale(1:round(1.4*size_r)),alphas_weight_radial_profile_old(1:round(1.4*size_r)),'b--')
plot(rho_tor_scale(1:round(1.4*size_r)),alphas_weight_radial_profile(1:round(1.4*size_r)))

plot(rho_tor_scale,alphas_weight_radial_profile_old,'b--')
plot(rho_tor_scale,alphas_weight_radial_profile)

weight_poly=polyfit(rho_tor_scale(1:round(1.5*size_r)),alphas_weight_radial_profile(1:round(1.5*size_r)),4);
alphas_weight_radial_profile_new=polyval(weight_poly,rho_tor_scale);
plot(rho_tor_scale,alphas_weight_radial_profile_new,'r')

alphas_weight_radial_profile_adjust=0.5*(alphas_weight_radial_profile_new+alphas_weight_radial_profile);
alphas_weight_radial_profile=0.4+0.6*alphas_weight_radial_profile_new;
alphas_weight_radial_profile(1:128)=min(alphas_weight_radial_profile(1:128),1.1);
alphas_weight_radial_profile(128:end)=min(alphas_weight_radial_profile(128:end),1.5);
% alphas_weight_radial_profile=alphas_weight_radial_profile_new;

% alphas_weight_radial_profile(1:15)=alphas_weight_radial_profile_adjust(1:15);
save weight_correction_profileP.mat radial_r_value_flux alphas_weight_radial_profile alphas_weight_radial_profile_old

plot(rho_tor_scale,alphas_weight_radial_profile,'g')

