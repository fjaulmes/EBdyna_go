% reset_data_analysis_environment
% filename=strcat('../../data_tokamak/pressure_profile.mat');
% load(filename);
% filename=strcat(DATA_FOLDER,'../../data_tokamak/flux_geometry.mat');
% load(filename);

close all

weight_particles_MBH=3.96e+12
% 
MBTEMP=500
load(strcat('initial_MBH_',num2str(MBTEMP),'_precession_stats_all.mat'));

INI_PRE_COLLAPSE_FILENAME=strcat('initial_MBH_',num2str(MBTEMP),'_pre_collapse_all.mat')

load(INI_PRE_COLLAPSE_FILENAME);

RECALCULATE_INI_WEIGHTS=0
SAVE_NEW_WEIGHTS=0
SAVE_FINAL_WEIGHTS=1

% 
Ti_value=MBTEMP*1e3;
if RECALCULATE_INI_WEIGHTS==1
    vtot_vector=sqrt(2*(eV/mHe)*alphas_Ekin);
    weight_vector_Ekin=4*pi*((mHe/(2*pi*Ti_value*eV))^(3/2)).*(vtot_vector.^2).*exp(-alphas_Ekin/Ti_value);
    weight_vector_Ekin=weight_vector_Ekin*1.5*Ti_value/mean(alphas_Ekin.*weight_vector_Ekin);
    
    alphas_weight_density=alphas_weight./weight_vector_Ekin;
    alphas_weight_Ekin=weight_vector_Ekin;
end

    disp('-------------------------------------')
    disp('statistics')
    disp('-------------------------------------')
    disp('trapped')
    tr_pop_nb=length(ALL_TRAPPED)/Nalphas_simulated;
    disp(tr_pop_nb)
    disp('co-passing')
    cop_pop_nb=length(CO_PASSING)/Nalphas_simulated;
    disp(cop_pop_nb)
    disp('counter-passing')
    counterp_pop_nb=length(COUNTER_PASSING)/Nalphas_simulated;
    disp(counterp_pop_nb)
    disp('potatoes')
    potatoes_pop_nb=length(POTATOES)/Nalphas_simulated;
    disp(potatoes_pop_nb)
    disp('stagnation')
    stagnation_pop_nb=length(STAGNATION)/Nalphas_simulated;
    disp(stagnation_pop_nb)
    disp('-------------------------------------')
    disp('total')
    disp(tr_pop_nb+cop_pop_nb+counterp_pop_nb+potatoes_pop_nb+stagnation_pop_nb)
    



alphas_weight=alphas_weight_Ekin.*alphas_weight_density;


NB_PART_SIM_TRANSP_RATIO=1.0


TE_profile_radial_ini=Te_profile;
NE_profile_radial_ini=Ne_profile;
TI_profile_radial_ini=Te_profile;
NI_profile_radial_ini=Ne_profile;
PTOT_profile_radial_ini=P_initial_profile;

% TE_profile_radial_ini=interp1((0:Nradial-1)/(Nradial-1),TE_profile_interp_ini,rho_tor_scale);
% TI_profile_radial_ini=interp1((0:Nradial-1)/(Nradial-1),TI_profile_interp_ini,rho_tor_scale);
% NE_profile_radial_ini=interp1((0:Nradial-1)/(Nradial-1),NE_profile_interp_ini,rho_tor_scale);
% NI_profile_radial_ini=interp1((0:Nradial-1)/(Nradial-1),NI_profile_interp_ini,rho_tor_scale);
% PTOT_profile_radial_ini=interp1((0:Nradial-1)/(Nradial-1),PTOT_profile_interp_ini,rho_tor_scale);
% 
Pbulk_ini= NI_profile_radial_ini.*TI_profile_radial_ini;

rho_tor_scale=sqrt(tor_flux_profile/max(tor_flux_profile));


nH0=1.44131e17
fast_density_profile=nH0*0.521298*exp(-(0.198739/0.298228)*tanh((rho_tor_scale-0.49123)/0.198739));
PI_MBH_ini=fast_density_profile*MBTEMP*1e3*eV;

% density_correction=(1e-6)*NB_PART_SIM_TRANSP_RATIO*weight_particles_MBH

%%



%EKIN_INF=1*1e3

% PART_POP=find((alphas_Ekin>EKIN_INF).*(alphas_pos_z<0.2).*(alphas_pos_z>-0.2));
%PART_POP=find((alphas_Ekin>EKIN_INF));

EKIN_BIN_SIZE=2000*1e3;
EKIN_BINS=(0:EKIN_BIN_SIZE:4000*1e3);
Ekin_values=EKIN_BINS(1:end-1)+0.5*EKIN_BIN_SIZE;

PICH_BIN_SIZE=1.0;
PITCH_BINS=(-1.0:PICH_BIN_SIZE:1.1);
pitch_values=PITCH_BINS(1:end-1)+0.5*PICH_BIN_SIZE;

R_BIN_SIZE=0.1;
R_BINS=(-1.0:R_BIN_SIZE:1.0)+R0+X_axis;
R_values=R_BINS(1:end-1)+0.5*R_BIN_SIZE;

Z_BIN_SIZE=0.1;
Z_BINS=(-1.0:Z_BIN_SIZE:1.0);
Z_values=Z_BINS(1:end-1)+0.5*Z_BIN_SIZE;

RADIAL_BIN_SIZE=0.03;
RADIAL_BINS=(0:RADIAL_BIN_SIZE:1.0);
RADIAL_values=RADIAL_BINS(1:end-1)+0.5*RADIAL_BIN_SIZE;
psi_radial_values=interp1(radial_r_value_flux,1:Nradial,RADIAL_values);
volume_radial_bins_values=interp1(radial_r_value_flux,volume_flux,RADIAL_BINS);
volume_radial_values=volume_radial_bins_values(2:end)-volume_radial_bins_values(1:end-1);
rho_tor_scale_values=interp1(radial_r_value_flux,rho_tor_scale,RADIAL_values);
rho_tor_scale_bins=interp1(radial_r_value_flux,rho_tor_scale,RADIAL_BINS);

%%

r_avg_gc=interp2(scale_X,scale_Z,psi_norm_XZsmall_map',pos_X_gc,pos_Z_gc,'*cubic');
rho_tor_gc=interp1(1:Nradial,rho_tor_scale,r_avg_gc);
% r_avg_gc=interp1(1:Nradial,radial_r_value_flux,r_avg);

rhist_ini=histc(r_avg,RADIAL_BINS);


l5=length(RADIAL_BINS)

trapped_fraction_profile_ini=zeros(l5,1);
trapped_fraction_profile_end=zeros(l5,1);
lambda_profile_ini=zeros(l5,1);
lambda_profile_end=zeros(l5,1);


density_correction_factor=(1e-6)*NB_PART_SIM_TRANSP_RATIO*weight_particles_MBH/(EKIN_BIN_SIZE)/PICH_BIN_SIZE

volume_RZ=zeros(length(R_values),length(Z_values));


for x=1:length(R_values)
    for z=1:length(Z_values)
        volume_RZ(x,z) = Z_BIN_SIZE*R_BIN_SIZE*2*pi*R_values(x);
    end
end


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

FILENAME_INI='dist_MBH_EpitchRZ_ini.dat'
% FILENAME_END='dist_MBH_EpitchRZ_200_end.dat'

string_title='# initial MBD distribution : dimensions ; scales (4 lines) ; [Ekin pitch R Z] '

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
% %         poloidal_hist_ini(eb,pb,:,:) = hist2d ([pos_X_gc(PART_POP)+R0 pos_Z_gc(PART_POP)], R_BINS, Z_BINS);
% %         data_array=density_correction_factor*squeeze(poloidal_hist_ini(eb,pb,:,:))./ volume_RZ;
% %         dist_EpitchRZ_ini(eb,pb,:,:) =data_array;
%         data_array=data_array';
%         if SAVE_INI==1
%             save(FILENAME_INI,'-append','data_array','-ASCII')
%         end
%     end
% end

%%
clear Ekin_profile_ini rhist_ini

for rb=1:length(RADIAL_values)-1
%     PART_POP=((r_avg>=RADIAL_BINS(rb)).*(r_avg<RADIAL_BINS(rb+1)));
    PART_POP=((rho_tor_gc>=rho_tor_scale_bins(rb)).*(rho_tor_gc<rho_tor_scale_bins(rb+1)));
    TR_BIN_POP=find(PART_POP.*ALL_TRAPPED_POP);
    PASS_BIN_POP=find(PART_POP.*ALL_PASSING_POP);
    trapped_fraction_profile_ini(rb)=length(TR_BIN_POP)/length(find(PART_POP));
    lambda_profile_ini(rb)=mean(abs(alphas_lambda(find(PART_POP.*(alphas_lambda>=0)))));
    lambda_profile_ini(rb)=mean(alphas_vpll(find(PART_POP.*(alphas_lambda>=0)))./alphas_vtot(find(PART_POP.*(alphas_lambda>=0))));
    Ekin_profile_ini(rb)=mean(alphas_weight_Ekin(find(PART_POP)).*alphas_Ekin(find(PART_POP)));
    rhist_ini(rb)=sum(alphas_weight_density(find(PART_POP)));
    P_ini(rb)=sum(alphas_weight(find(PART_POP)).*alphas_Ekin(find(PART_POP)));
end

Ekin_profile_ini(rb+1)=Ekin_profile_ini(rb);
rhist_ini(rb+1)=rhist_ini(rb);
P_ini(rb+1)=P_ini(rb);

% %%
% figure(2)
% set(gca,'fontsize',20)
% contourf(R_values,Z_values,(squeeze(sum(sum(dist_EpitchRZ_ini(1:end,:,:,:),2),1))'),12);axis xy;hold on
% contour(scale_X+R0,scale_Z,psi_norm_XZsmall_map','k','linewidth',2)

% %%
% figure(3)
% hold on
% grid on
% plot(R_values,(1e6)*EKIN_BIN_SIZE*PICH_BIN_SIZE*(squeeze(sum(sum(sum(dist_EpitchRZ_ini(1:end,:,:,floor(0.5*(length(Z_values))-2):ceil(0.5*(length(Z_values)))+2),2),1),4)/5)'),'r')
% % plot(R_values,(1e6)*EKIN_BIN_SIZE*PICH_BIN_SIZE*(squeeze(sum(sum(sum(dist_EpitchRZ_end(1:end,:,:,floor(0.5*(length(Z_values))-2):ceil(0.5*(length(Z_values)))+2),2),1),4)/5)'),'b')



%%
figure(8)

hold on

density_MBH_ini=NB_PART_SIM_TRANSP_RATIO*weight_particles_MBH*rhist_ini./volume_radial_values;
% density_NBI_end=rhist_end(1:end-1)./volume_radial_values';
P_ini=(2/3)*NB_PART_SIM_TRANSP_RATIO*weight_particles_MBH*P_ini*eV;

PMBH_ini=(2/3)*eV*density_MBH_ini.*Ekin_profile_ini;
% PNBI_end=(2/3)*eV*NB_PART_SIM_TRANSP_RATIO*weight_particles_MBH*density_NBI_end.*Ekin_profile_end';


% PI_ini=TI_profile_interp_ini.*NI_profile_interp_ini*eV;
% PI_end=TI_profile_interp_end.*NI_profile_interp_end*eV;


% PMBH_ini_rholin=interp1(rho_tor_scale_values,PMBH_ini,(0:257-1)/(257-1));
PMBH_ini_rho=interp1(rho_tor_scale_values,PMBH_ini,rho_tor_scale);
PMBH_ini_rho(isnan(PMBH_ini_rho))=0;

density_ini_rho=interp1(rho_tor_scale_values,density_MBH_ini,rho_tor_scale);
density_ini_rho(isnan(density_ini_rho))=0;


%
subplot(2,1,1)

hold on
grid on
plot(rho_tor_scale,PMBH_ini_rho,'b--')
plot(rho_tor_scale,PI_MBH_ini,'r')
% plot(rho_tor_scale_values,PNBI_end,'r--')
xlabel('\rho_{tor}')
% ylabel('P/P_{tot}')
ylabel('P_{NBI} (Pa)')

legend('EBdyna ini','P_i profile')

subplot(2,1,2)

hold on
grid on
plot(rho_tor_scale,density_ini_rho,'b--')
plot(rho_tor_scale,fast_density_profile,'r')
% plot(rho_tor_scale_values,PNBI_end,'r--')
xlabel('\rho_{tor}')


%%




PMBH_ini_rho=interp1(rho_tor_scale_values,PMBH_ini,rho_tor_scale);
PMBH_ini_rho(isnan(PMBH_ini_rho))=0;
n_ini_rho=interp1(rho_tor_scale_values,density_MBH_ini,rho_tor_scale);
E_ini_rho=interp1(rho_tor_scale_values,Ekin_profile_ini,rho_tor_scale);

alphas_weight_rho_profile=fast_density_profile./n_ini_rho;
alphas_weight_rho_profile2=PI_MBH_ini./PMBH_ini_rho;
alphas_weight_rho_profile=0.5*(alphas_weight_rho_profile+alphas_weight_rho_profile2);
alphas_weight_rho_profile_Ekin=(1.5*Ti_value)./E_ini_rho;
alphas_weight_rho_profile(1:12)=min(alphas_weight_rho_profile);


alphas_weight_rho_profile(isnan(alphas_weight_rho_profile))=max(alphas_weight_rho_profile);
alphas_weight_rho_profile(1:5)=alphas_weight_rho_profile(5);
alphas_weight_rho_profile_Ekin(isnan(alphas_weight_rho_profile))=max(alphas_weight_rho_profile_Ekin);
alphas_weight_rho_profile_Ekin(1:12)=min(alphas_weight_rho_profile_Ekin);

alphas_weight_radial_profile=interp1((0:257-1)/(257-1),alphas_weight_rho_profile,rho_tor_scale);
alphas_weight_radial_profile(isnan(alphas_weight_radial_profile))=1;
alphas_weight_radial_profile(isinf(alphas_weight_radial_profile))=1;
alphas_weight_radial_profile(alphas_weight_radial_profile==0)=1;
alphas_weight_radial_Ekin_interp=interp1((0:257-1)/(257-1),alphas_weight_rho_profile_Ekin,rho_tor_scale);
alphas_weight_radial_Ekin_interp(isnan(alphas_weight_radial_Ekin_interp))=1;
alphas_weight_radial_Ekin_interp(isinf(alphas_weight_radial_Ekin_interp))=1;
alphas_weight_radial_Ekin_interp(alphas_weight_radial_Ekin_interp==0)=1;

% alphas_weight_radial_Ekin_interp=0.5*+0.5*alphas_weight_radial_Ekin_interp;

alphas_weight_radial_profile(1:127)=min(alphas_weight_radial_profile(1:127),2.2);
alphas_weight_radial_profile(128:257)=min(alphas_weight_radial_profile(128:257),2.2);

radial_profile_coefs=polyfit(rho_tor_scale,alphas_weight_radial_profile,7);
alphas_weight_radial_fit=polyval(radial_profile_coefs,rho_tor_scale);
alphas_weight_radial_profile=alphas_weight_radial_fit;
[polymax_value polymax ]=max(alphas_weight_radial_fit)
alphas_weight_radial_profile(polymax:end)=polymax_value;
[polymin_value polymin ]=min(alphas_weight_radial_fit)
alphas_weight_radial_profile(1:polymin)=polymin_value;

% alphas_weight_radial_profile=0.1+0.9*alphas_weight_radial_profile;



if SAVE_NEW_WEIGHTS==1
    
    save weight_correction_profileP.mat radial_r_value_flux alphas_weight_radial_profile
    r_avg_gc=interp1(1:length(radial_r_value_flux),radial_r_value_flux,r_avg_gc);
    r_avg_avg=0.5*(r_avg+r_avg_gc);
    
    alphas_weight_density2=interp1(radial_r_value_flux,alphas_weight_radial_profile,r_avg_avg);
    alphas_weight_Ekin2=interp1(radial_r_value_flux,alphas_weight_radial_Ekin_interp,r_avg_avg);
%     alphas_weight_density=alphas_weight_density2;
%     alphas_weight_Ekin=alphas_weight_Ekin2;
    alphas_weight_Ekin2(isnan(alphas_weight_Ekin2))=1;
    
    alphas_weight_density=alphas_weight_density.*(0.1+0.9*alphas_weight_density2);
%     alphas_weight_density=alphas_weight_density.*alphas_weight_Ekin2;
    alphas_weight_density(isnan(alphas_weight_density))=1;
    
    save(INI_PRE_COLLAPSE_FILENAME,'-append','alphas_weight_density','alphas_weight_Ekin')
    
end



if SAVE_FINAL_WEIGHTS==1
    alphas_weight=alphas_weight_density.*alphas_weight_Ekin;
    save(INI_PRE_COLLAPSE_FILENAME,'-append','alphas_weight')
end



