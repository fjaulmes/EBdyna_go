reset_data_analysis_environment
% close all
    [XXsmall ZZsmall]=meshgrid(scale_X,scale_Z);
    finesse_data_X=reshape((Rpos_PR_map(:,1:size_r)-R0),NP*size_r,1);
    finesse_data_Z=reshape(Z_PR_map(:,1:size_r),NP*size_r,1);

weight_transp_W=8e12

NBFRAMES_MOVIE=99
FRAMETSSIZE=200

load('initial_W40_precession_stats_all.mat');
NB_PART_SIM_TRANSP_RATIO=1


%%
load('initial_W40_pre_collapse_all.mat');


figure(1)

%EKIN_INF=1*1e3

% PART_POP=find((alphas_Ekin>EKIN_INF).*(alphas_pos_z<0.2).*(alphas_pos_z>-0.2));
%PART_POP=find((alphas_Ekin>EKIN_INF));

EKIN_BIN_SIZE=80*1e3;
EKIN_BINS=(0:EKIN_BIN_SIZE:80*1e3);
Ekin_values=EKIN_BINS(1:end-1)+0.5*EKIN_BIN_SIZE;

PICH_BIN_SIZE=2.0;
PITCH_BINS=(-1.0:PICH_BIN_SIZE:1.0);
pitch_values=PITCH_BINS(1:end-1)+0.5*PICH_BIN_SIZE;

R_BIN_SIZE=0.03;
R_BINS=(-0.5:R_BIN_SIZE:0.5)+R0;
R_values=R_BINS(1:end-1)+0.5*R_BIN_SIZE;

Z_BIN_SIZE=0.03;
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


density_correction_factor=weight_transp_W/(EKIN_BIN_SIZE)/PICH_BIN_SIZE

volume_RZ=zeros(length(R_values),length(Z_values));


for x=1:length(R_values)
    for z=1:length(Z_values)
        volume_RZ(x,z) = Z_BIN_SIZE*R_BIN_SIZE*2*pi*R_values(x);
    end
end

%R_BIN_SIZE=0.06;
%R_BINS=(0:R_BIN_SIZE:0.6);
%r_values=R_BINS(1:end-1)+0.5*R_BIN_SIZE;
ROTATION_RATE=pi/8
PHI_MIN_ini=(1+0.0)*2*ROTATION_RATE
PHI_MAX_ini=(2)*2*ROTATION_RATE
time_step=0

for frame_number=200:FRAMETSSIZE:19800
    %%
    f_rank=round(frame_number/FRAMETSSIZE);
    time_scale_movie(f_rank)=((frame_number-FRAMETSSIZE)/FRAMETSSIZE);
    PHI_MIN=PHI_MIN_ini+(frame_number/FRAMETSSIZE)*(ROTATION_RATE);
    PHI_MAX=PHI_MAX_ini+(frame_number/FRAMETSSIZE)*(ROTATION_RATE);
    PHI_MIN=wrap2pi(PHI_MIN);
    PHI_MAX=wrap2pi(PHI_MAX)
    if (PHI_MAX<(2*ROTATION_RATE)-1e-4) && (PHI_MIN<2*pi-1e-4) && (PHI_MIN>2*pi-(ROTATION_RATE)-1e-4)
        INTER_FRAME=1
        PHI_MIN1=PHI_MIN;
        PHI_MAX1=2*pi;
        PHI_MIN2=0;
        PHI_MAX2=PHI_MAX;    
        PHI_AVG_POS=2*pi
    else
        INTER_FRAME=0;
        PHI_AVG_POS=0.5*(PHI_MIN+PHI_MAX)
    end
    if PHI_MAX<1e-4
         INTER_FRAME=0;
         PHI_MAX=2*pi;
         PHI_AVG_POS=0.5*(PHI_MIN+PHI_MAX)
    end
    PHI_AVG_RANK=round(257*PHI_AVG_POS/(2*pi))+1
    phi_avg_movie(f_rank)=PHI_AVG_POS;
    
    % post simulation results
    INPUTNAME=strcat('./W40_fc1h1_',num2str(frame_number),'_all.mat')
    load(INPUTNAME);
    
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
    poloidal_hist_frame=zeros(length(Ekin_values),length(pitch_values),length(R_values),length(Z_values));
    dist_EpitchRZ_ini=poloidal_hist_frame;
    alphas_pos_phi_wrapped=wrap2pi(alphas_pos_phi);
    Rpos_zb=zeros(length(Z_values),1);
    Zpos_rb=zeros(length(R_values),1);
    
    for eb=1:length(Ekin_values)
        for pb=1:length(pitch_values)
            %             if PHI_MIN<PHI_MAX
            if INTER_FRAME==0
                PART_POP=find((alphas_pos_phi_wrapped>=PHI_MIN).*(alphas_pos_phi_wrapped<=PHI_MAX).*(alphas_Ekin>100).*(alphas_Ekin<EKIN_BINS(eb+1)).*(alphas_pitch>=PITCH_BINS(pb)).*(alphas_pitch<PITCH_BINS(pb+1)));
            else
                PART_POP1=((alphas_pos_phi_wrapped>=PHI_MIN1).*(alphas_pos_phi_wrapped<=PHI_MAX1).*(alphas_Ekin>100).*(alphas_Ekin<EKIN_BINS(eb+1)).*(alphas_pitch>=PITCH_BINS(pb)).*(alphas_pitch<PITCH_BINS(pb+1)));
                PART_POP2=((alphas_pos_phi_wrapped>=PHI_MIN2).*(alphas_pos_phi_wrapped<=PHI_MAX2).*(alphas_Ekin>100).*(alphas_Ekin<EKIN_BINS(eb+1)).*(alphas_pitch>=PITCH_BINS(pb)).*(alphas_pitch<PITCH_BINS(pb+1)));
                PART_POP=find((PART_POP1+PART_POP2)>0);
            end
            %             else
            %             PART_POP=find((alphas_pos_phi>PHI_MIN).*(alphas_pos_phi>PHI_MAX).*(alphas_Ekin>=EKIN_BINS(eb)).*(alphas_Ekin<EKIN_BINS(eb+1)).*(alphas_pitch>=PITCH_BINS(pb)).*(alphas_pitch<PITCH_BINS(pb+1)));
            %             end
            poloidal_hist_frame(eb,pb,:,:) = hist2d ([pos_X_gc(PART_POP)+R0 pos_Z_gc(PART_POP)], R_BINS, Z_BINS);
            data_array=squeeze(poloidal_hist_frame(eb,pb,:,:))./ volume_RZ;
            data_array_corr=density_correction_factor*data_array.^(1.5);
            dist_EpitchRZ_ini(eb,pb,:,:)=data_array_corr;
            for zb=1:length(Z_values)
                [max_value Rpos_value ]=max(data_array_corr(:,zb));
                Rpos_zb(zb)=max_value;
            end
            [max_value Zpos_zb_rank]=max(Rpos_zb);
            Zcore_movie(f_rank)=Z_values(Zpos_zb_rank);
            for rb=1:length(R_values)
                [max_value Zpos_value ]=max(data_array_corr(rb,:));
                Zpos_rb(rb)=max_value;
            end
            [max_value Rpos_rb_rank]=max(Zpos_rb);
            Rcore_movie(f_rank)=R_values(Rpos_rb_rank);
            disp('(R,Z) core ==')
            disp(R_values(Rpos_rb_rank));
            disp(Z_values(Zpos_zb_rank));
        end
    end
    
%     for rb=1:l5
%         PART_POP=((r_avg>=RADIAL_BINS(rb)).*(r_avg<RADIAL_BINS(rb+1)));
%         TR_BIN_POP=find(PART_POP.*ALL_TRAPPED_POP);
%         PASS_BIN_POP=find(PART_POP.*ALL_PASSING_POP);
%         trapped_fraction_profile_ini(rb)=length(TR_BIN_POP)/length(find(PART_POP));
%         lambda_profile_ini(rb)=mean(abs(alphas_lambda(find(PART_POP.*(alphas_lambda>=0)))));
%         lambda_profile_ini(rb)=mean(alphas_vpll(find(PART_POP.*(alphas_lambda>=0)))./alphas_vtot(find(PART_POP.*(alphas_lambda>=0))));
%     end
%     
    clf
    fig=figure(1)
    
    contourf(R_values,Z_values,(squeeze(sum(dist_EpitchRZ_ini(1,:,:,:),2))'),8);axis xy;hold on
%     colorbar
    pause(0.2)

    contour(scale_X+R0,scale_Z,psi_norm_XZsmall_map','k','linewidth',2)
    
    xlim([1.38 1.95])
    ylim([-0.38 0.38])
    
    psi_star_omega_map_half(:,:)=psi_star_2D_evol_lin(round(100*time_step/NBFRAMES_MOVIE)+1,:,:);
    
    % Using symmetry to reconstruct a poloidal turn
    psi_star_omega_map=zeros(size_r,NB_THETA);
    psi_star_omega_map(:,1:round(0.5*NB_THETA))=psi_star_omega_map_half(:,:);
    psi_star_omega_map(:,round(0.5*NB_THETA):NB_THETA)=psi_star_omega_map_half(:,round(0.5*NB_THETA):-1:1);
    
    psi_star_map_phi_rank=psi_star_omega_map*0;
    psi_star_map_phi_rank(:,:)=[psi_star_omega_map(:,PHI_AVG_RANK:NB_THETA)  psi_star_omega_map(:,1:PHI_AVG_RANK-1)];
    psi_star_map_phi_rank=psi_star_map_phi_rank(:,NB_THETA:-1:1);
    psi_star_map_phi_rank=psi_star_map_phi_rank';
    psidata=reshape(psi_star_map_phi_rank(:,1:size_r),NB_THETA*size_r,1);
    
    psistarmap=griddata(finesse_data_X,finesse_data_Z,psidata,XXsmall,ZZsmall,'cubic');
    psistarmap=psistarmap';
    psistarmap(isnan(psistarmap))=0;
    contour(scale_X+R0,scale_Z-Z_axis,psistarmap',(-9.01:1.0:0.05)*1e-3,'g','linewidth',3);
%     contour(scale_X+R0,scale_Z-Z_axis,psistarmap',7,'g','linewidth',3);
    
    time_step=time_step+1
    pause(0.2)
    frame_name='t';
    if (time_step<10)
        frame_name=strcat(frame_name,'0');
    end
    frame_name=strcat(frame_name,num2str(time_step));
    filename='cartoon/'
    filename=strcat(filename,frame_name,'.png');
   
    F=getframe(fig);
    [h, w, p] = size(F.cdata);  % use 1st frame to get dimensions

    [im,map] = frame2im(F);    %Return associated image data
    imwrite(im,filename,'png');
    saveas(gcf,filename,'png');%saving freq plot

    h
    w
    %%
    clc
    pause(0.2)
    

    
end

%%
psi_core_evol=interp2(scale_X,scale_Z,psi_norm_XZsmall_map',Rcore_movie-R0,Zcore_movie);
psi_value_evol=interp1(1:Nradial,psi_scale,psi_core_evol);
r_value_core_evol=interp1(1:Nradial,radial_r_value_flux,psi_core_evol);
rhopol_recalc_core_evol=sqrt(1-psi_value_evol/psi_global)


save Wpos_fc1_core_movie_piover8_fm_Elim1000.mat time_scale_movie Rcore_movie Zcore_movie phi_avg_movie r_value_core_evol

SAWTOOTH_DURATION=4e-4
SAWTOOTH_DURATION_ADJUST=0.95;
SAWTOOTH_DURATION=SAWTOOTH_DURATION*SAWTOOTH_DURATION_ADJUST

%
%%
figure(1)
load('tomas_SXR_data.mat')

subplot(2,1,1)
set(gca,'fontsize',16)

hold on
grid on
plot(SAWTOOTH_DURATION*time_scale_movie/100,Rcore_movie-R0,'b')
plot(time_core_evol-time_core_evol(56),R_SXR_evol-R0,'g--','linewidth',2)
ylabel('X')
xlabel('time (s)')
xlim([-0.1 1.1]*SAWTOOTH_DURATION);
legend('EBdyna','SXR')

subplot(2,1,2)
set(gca,'fontsize',16)

hold on
grid on

ylabel('Z')
plot(SAWTOOTH_DURATION*time_scale_movie/100,Zcore_movie,'r')
plot(time_core_evol-time_core_evol(56),Z_SXR_evol,'g--','linewidth',2)
xlabel('time (s)')
xlim([-0.1 1.1]*SAWTOOTH_DURATION);
legend('EBdyna','SXR')

%%
figure(2)
set(gca,'fontsize',16)
hold on
grid on

plot(time_core_evol-time_core_evol(56),r_value_evol,'r','linewidth',2)
plot(SAWTOOTH_DURATION*time_scale_movie/100,r_value_core_evol,'k')
plot(time_scale_lin*SAWTOOTH_DURATION,ksi0_evol_lin,'g--','linewidth',2)
legend('SXR','EBdyna (W40)','EBdyna (\xi)' )

xlim([-0.1 1.1]*SAWTOOTH_DURATION);
ylabel('\xi (m)')
xlabel('time (s)')