reset_data_analysis_environment
close all

weight_transp_W=8e10


load('initial_MB_W_precession_stats_all.mat');
NB_PART_SIM_TRANSP_RATIO=1


%%
load('initial_MB_W_pre_collapse_all.mat');



%EKIN_INF=1*1e3

% PART_POP=find((alphas_Ekin>EKIN_INF).*(alphas_pos_z<0.2).*(alphas_pos_z>-0.2));
%PART_POP=find((alphas_Ekin>EKIN_INF));

EKIN_BIN_SIZE=80*1e3;
EKIN_BINS=(0:EKIN_BIN_SIZE:160*1e3);
Ekin_values=EKIN_BINS(1:end-1)+0.5*EKIN_BIN_SIZE;

PICH_BIN_SIZE=2.0;
PITCH_BINS=(-1.0:PICH_BIN_SIZE:1.0);
pitch_values=PITCH_BINS(1:end-1)+0.5*PICH_BIN_SIZE;

R_BIN_SIZE=0.02;
R_BINS=(-0.5:R_BIN_SIZE:0.5)+R0;
R_values=R_BINS(1:end-1)+0.5*R_BIN_SIZE;

Z_BIN_SIZE=0.05;
Z_BINS=(-0.8:Z_BIN_SIZE:0.8);
Z_values=Z_BINS(1:end-1)+0.5*Z_BIN_SIZE;

RADIAL_BIN_SIZE=0.05;
RADIAL_BINS=(0:R_BIN_SIZE:0.5);
RADIAL_values=RADIAL_BINS(1:end-1)+0.5*RADIAL_BIN_SIZE;


l1=length(Ekin_values) 
l2=length(pitch_values)
l3=length(R_values) 
l4=length(Z_values)
l5=length(RADIAL_values)

trapped_fraction_profile_ini=zeros(l5,1);
trapped_fraction_profile_end=zeros(l5,1);
lambda_profile_ini=zeros(l5,1);
lambda_profile_end=zeros(l5,1);


density_correction_factor=(1e-6)*weight_transp_W/(EKIN_BIN_SIZE)/PICH_BIN_SIZE

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
alphas_pitch=alphas_vpll./alphas_vtot;
alphas_lambda_ini=atan(alphas_vperp./alphas_vpll);
alphas_lambda0_ini=Bavg*alphas_mm./alphas_Ekin;
alphas_r=interp1(1:257,radial_r_value_flux,alphas_psi);
alphas_Ekin_ini=alphas_Ekin;

% here can try alphas_r or r_avg or even pphi [if small range of Ekin]
poloidal_hist_ini=zeros(length(Ekin_values),length(pitch_values),length(R_values),length(Z_values));
dist_EpitchRZ_ini=poloidal_hist_ini;

FILENAME_INI='dist_W_AUG30382_2p3_EpitchRZ_ini.dat'
FILENAME_END='dist_W_AUG30382_2p3_EpitchRZ_end.dat'

string_title='# initial W distribution : dimensions ; scales (4 lines) ; [Ekin pitch R Z] '

%save dist_EpitchRZ_ini.dat -append  string_title -ASCII
save(FILENAME_INI,'-append','l1','l2','l3','l4','-ASCII')
save(FILENAME_INI,'-append','Ekin_values','-ASCII')
save(FILENAME_INI,'-append','pitch_values','-ASCII')
save(FILENAME_INI,'-append','R_values','-ASCII')
save(FILENAME_INI,'-append','Z_values','-ASCII')

for eb=1:length(Ekin_values)
    for pb=1:length(pitch_values)
        PART_POP=find((alphas_Ekin>=EKIN_BINS(eb)).*(alphas_Ekin<EKIN_BINS(eb+1)).*(alphas_pitch>=PITCH_BINS(pb)).*(alphas_pitch<PITCH_BINS(pb+1)));
        poloidal_hist_ini(eb,pb,:,:) = hist2d ([pos_X_gc(PART_POP)+R0 pos_Z_gc(PART_POP)], R_BINS, Z_BINS);
        data_array=density_correction_factor*squeeze(poloidal_hist_ini(eb,pb,:,:))./ volume_RZ;
        dist_EpitchRZ_ini(eb,pb,:,:) =data_array;
        data_array=data_array';
        save(FILENAME_INI,'-append','data_array','-ASCII')

    end
end



for rb=1:l5
    PART_POP=((r_avg>=RADIAL_BINS(rb)).*(r_avg<RADIAL_BINS(rb+1)));
    TR_BIN_POP=find(PART_POP.*ALL_TRAPPED_POP);
    PASS_BIN_POP=find(PART_POP.*ALL_PASSING_POP);
    trapped_fraction_profile_ini(rb)=length(TR_BIN_POP)/length(find(PART_POP));
    lambda_profile_ini(rb)=mean(abs(alphas_lambda(find(PART_POP.*(alphas_lambda>=0)))));
    lambda_profile_ini(rb)=mean(alphas_vpll(find(PART_POP.*(alphas_lambda>=0)))./alphas_vtot(find(PART_POP.*(alphas_lambda>=0))));
end

%%
% pcolor (r_values, pitch_values, mHist);
% shading faceted; 



load('MB_W_m_fc1p6h1_all.mat');
%load('final_NBI60keV_precession_stats_all.mat');
r_avg=interp2(scale_X,scale_Z,radial_XZsmall_map',pos_X_gc,pos_Z_gc);

alphas_B_end=interp2(scale_X,scale_Z,Btot_XZ_map',pos_X_gc,pos_Z_gc,'*linear');
alphas_theta_end=interp2(scale_X,scale_Z,theta_XZsmall_map',pos_X_gc,pos_Z_gc,'*linear');

alphas_Eperp=max(alphas_Ekin-0.5*(mHe/eV)*alphas_vpll.^2,0);
alphas_vperp=sqrt(2*(eV/mHe)*alphas_Eperp);
alphas_vtot=sqrt(2*(eV/mHe)*alphas_Ekin);
alphas_pitch=alphas_vpll./alphas_vtot;
alphas_lambda_end=atan(alphas_vperp./alphas_vpll);
alphas_r=interp1(1:257,radial_r_value_flux,alphas_psi);
alphas_lambda0_end=Bavg*alphas_mm./alphas_Ekin;
alphas_Ekin_end=alphas_Ekin;

poloidal_hist_end=zeros(length(Ekin_values),length(pitch_values),length(R_values),length(Z_values));
dist_EpitchRZ_end=poloidal_hist_end;

string_title='# final W distribution : dimensions ; scales (4 lines) ; [Ekin pitch R Z] '
%save dist_EpitchRZ_end.dat -append string_title -ASCII
save(FILENAME_END,'-append','l1','l2','l3','l4','-ASCII')
save(FILENAME_END,'-append','Ekin_values','-ASCII')
save(FILENAME_END,'-append','pitch_values','-ASCII')
save(FILENAME_END,'-append','R_values','-ASCII')
save(FILENAME_END,'-append','Z_values','-ASCII')


for eb=1:length(Ekin_values)
    for pb=1:length(pitch_values)
        PART_POP=find((alphas_Ekin>=EKIN_BINS(eb)).*(alphas_Ekin<EKIN_BINS(eb+1)).*(alphas_pitch>=PITCH_BINS(pb)).*(alphas_pitch<PITCH_BINS(pb+1)));
        poloidal_hist_end(eb,pb,:,:) = hist2d ([pos_X_gc(PART_POP)+R0 pos_Z_gc(PART_POP)], R_BINS, Z_BINS);
        data_array=density_correction_factor*squeeze(poloidal_hist_end(eb,pb,:,:))./ volume_RZ;
        dist_EpitchRZ_end(eb,pb,:,:) =data_array;  % cm3
        data_array=data_array';
        save(FILENAME_END,'-append','data_array','-ASCII')
    end
end

for rb=1:l5
    PART_POP=((r_avg>=RADIAL_BINS(rb)).*(r_avg<RADIAL_BINS(rb+1)));
    TR_BIN_POP=find(PART_POP.*ALL_TRAPPED_POP);
    PASS_BIN_POP=find(PART_POP.*ALL_PASSING_POP);
    trapped_fraction_profile_end(rb)=length(TR_BIN_POP)/length(find(PART_POP));
    lambda_profile_end(rb)=mean(abs(alphas_lambda(find(PART_POP.*(alphas_lambda>=0)))));
    lambda_profile_end(rb)=mean(alphas_vpll(find(PART_POP.*(alphas_lambda>=0)))./alphas_vtot(find(PART_POP.*(alphas_lambda>=0))));
end

save large_crash_NBI_distributions.mat dist_EpitchRZ_ini dist_EpitchRZ_end Ekin_values pitch_values R_values Z_values trapped_fraction_profile_ini trapped_fraction_profile_end RADIAL_values


%%

figure(1)
contourf(R_values,Z_values,(squeeze(sum(poloidal_hist_end(1,:,:,:),2))'),8);axis xy;hold on
contour(scale_X+R0,scale_Z,psi_norm_XZsmall_map','k')

figure(2)
contourf(R_values,Z_values,(squeeze(sum(poloidal_hist_ini(1,:,:,:),2))'),8);axis xy;hold on
contour(scale_X+R0,scale_Z,psi_norm_XZsmall_map','k')