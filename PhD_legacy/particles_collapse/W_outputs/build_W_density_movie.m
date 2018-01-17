reset_data_analysis_environment
close all

weight_transp_W=8e10


load('initial_MB_W40_precession_stats_all.mat');
NB_PART_SIM_TRANSP_RATIO=1


%%
load('initial_MB_W40_pre_collapse_all.mat');


figure(1)

%EKIN_INF=1*1e3

% PART_POP=find((alphas_Ekin>EKIN_INF).*(alphas_pos_z<0.2).*(alphas_pos_z>-0.2));
%PART_POP=find((alphas_Ekin>EKIN_INF));

EKIN_BIN_SIZE=80*1e3;
EKIN_BINS=(0:EKIN_BIN_SIZE:160*1e3);
Ekin_values=EKIN_BINS(1:end-1)+0.5*EKIN_BIN_SIZE;

PICH_BIN_SIZE=2.0;
PITCH_BINS=(-1.0:PICH_BIN_SIZE:1.0);
pitch_values=PITCH_BINS(1:end-1)+0.5*PICH_BIN_SIZE;

R_BIN_SIZE=0.04;
R_BINS=(-0.5:R_BIN_SIZE:0.5)+R0;
R_values=R_BINS(1:end-1)+0.5*R_BIN_SIZE;

Z_BIN_SIZE=0.04;
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
PHI_MIN=pi/8
PHI_MAX=4*pi/8
time_step=1

for frame_number=640:640:32000
    PHI_MIN=wrap2pi(PHI_MIN+pi/4)
    PHI_MAX=wrap2pi(PHI_MAX+pi/4)
    if PHI_MAX<PHI_MIN
%         PHI_MIN=0
        PHI_MAX=2*pi
    end
    % post simulation results
    INPUTNAME=strcat('./MB_W40_fc1p6h2_',num2str(frame_number),'_all.mat')
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
    
    
    for eb=1:length(Ekin_values)
        for pb=1:length(pitch_values)
            PART_POP=find((wrap2pi(alphas_pos_phi)>PHI_MIN).*(wrap2pi(alphas_pos_phi)<PHI_MAX).*(alphas_Ekin>=EKIN_BINS(eb)).*(alphas_Ekin<EKIN_BINS(eb+1)).*(alphas_pitch>=PITCH_BINS(pb)).*(alphas_pitch<PITCH_BINS(pb+1)));
%             PART_POP=find((alphas_Ekin>=EKIN_BINS(eb)).*(alphas_Ekin<EKIN_BINS(eb+1)).*(alphas_pitch>=PITCH_BINS(pb)).*(alphas_pitch<PITCH_BINS(pb+1)));
            poloidal_hist_frame(eb,pb,:,:) = hist2d ([pos_X_gc(PART_POP)+R0 pos_Z_gc(PART_POP)], R_BINS, Z_BINS);
            data_array=density_correction_factor*squeeze(poloidal_hist_frame(eb,pb,:,:))./ volume_RZ;
            dist_EpitchRZ_ini(eb,pb,:,:)=data_array;
            data_array=data_array';
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
    
    pause(0.1)
    figure(1)
    clf
    
    contourf(R_values,Z_values,(squeeze(sum(dist_EpitchRZ_ini(1,:,:,:),2))'),8);axis xy;hold on
    contour(scale_X+R0,scale_Z,psi_norm_XZsmall_map','k','linewidth',2)
    xlim([1.35 2.05])
    ylim([-0.45 0.4])
    colorbar
    time_step=time_step+1
    pause(0.1)
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

    clc
    h
    w
end
