TITLESTRING='CU #5400 '
% clear all

try
load('../../../data_common/wallRZ_CU.mat')
catch
load('../data_common/wallRZ_CU.mat')
end

SAVE_NEW_PREC_INPUT=0;

PLOT_BASICS=1;
PLOT_SD=1;
PLOT_NEUTRALS_RELOC=0;
PLOT_CX_POWER_2D=0;
PLOT_DOTS=0;
PLOT_HIST=0;
PLOT_CXLOSS=0;
PLOT_NDD=0;
PLOT_2D_LOSSES=1;

PLOT_HIST_CX_SUMMARY=0;

if ~exist('par')
    error('Please load an EBdyna output first!')
end
if ~exist('hist2d')
    try
    run('../../../startup_EBdyna_go.m')
    catch
    error('Configure by running startup_EBdyna_go...')
    end
end

if ~isfield(par,'FILENAME')
    warning('Filename is undefined! Using a default name to save summary data...')
    par.FILENAME='EBdyna_summary_output_default.mat';
end




try
load('../../data_tokamak/physics_constants.mat');
load('../../data_tokamak/XZsmall_fields_tokamak_pre_collapse.mat');
load('../../data_tokamak/motions_map_dimensions.mat');
catch
load('./data_tokamak/physics_constants.mat');
load('./data_tokamak/XZsmall_fields_tokamak_pre_collapse.mat');
load('./data_tokamak/motions_map_dimensions.mat');
end

if strcmp(par.MACHINE_NAME,'MAST-U')
    try
        WALL_TOK=load('../../../data_common/MASTU_vessel_limits.mat');
    catch
        WALL_TOK=load('../data_common/MASTU_vessel_limits.mat');
    end
elseif strcmp(par.MACHINE_NAME,'COMPASS-U')
    try
        WALL_TOK=load('../../../data_common/CU_vessel_limits_standard_no_rescale.mat');
    catch
        WALL_TOK=load('../data_common/CU_vessel_limits_standard_no_rescale.mat');
    end
elseif strcmp(par.MACHINE_NAME,'COMPASS')
    try
        WALL_TOK=load('../../../data_common/COMPASS_vessel_limits.mat');
    catch
        WALL_TOK=load('../data_common/COMPASS_vessel_limits.mat');
    end
end


NB_TS_AVG=par.NB_TS_RECORDED-1;


%%
if PLOT_BASICS
hf=figure;
hold on
grid on
surf(summary.Ekin_values,summary.pitch_values,summary.dist_vspace'); shading interp; view(2)
xlabel('Ekin [eV]')
ylabel('pitch')

set(gca,'fontsize',16)

title('EBdyna')

end



%%


if PLOT_SD
    
    figure(11);
    %
    
    set(gca,'fontsize',20)
    hold on
    grid on
    plot(summary.Ekin_values, summary.dist_Ekin,'linewidth',2)
    xlabel('Ekin [eV]')
end



%%
disp('statistics for the total number of time stamps recoreded in this file');



disp('---------------------------------------------------------------')


disp(['POWER ION_WALL_LOSS   = ' num2str(1e-3*round(pfc_losses.POWER_ION_LOSS_TOT))]);
disp(['POWER_OWL_LOSS  = ' num2str(1e-3*round(pfc_losses.POWER_OWL_LOSS))]);
disp(['POWER_HFS_MIS_LOSS  = ' num2str(1e-3*round(pfc_losses.POWER_HFS_MID_LOSS))]);
disp(['POWER_TOP_LOSS  = ' num2str(1e-3*round(pfc_losses.POWER_TOP_LOSS))]);
disp(['POWER_DOWN_LOSS = ' num2str(1e-3*round(pfc_losses.POWER_DOWN_LOSS))]);
disp(['POWER_CX_LOSS   = ' num2str(1e-3*round(pfc_losses.POWER_CX_LOSS))]);
disp(['DID CX (POWER) = ' num2str(1e-3*round(pfc_losses.POW_DID_CX_DT))]);






%%
if par.CALCULATE_CX

    
    disp('---------------------------------------------------------------')
    
    disp(['CX_POWER_OWL_LOSS  = ' num2str(1e-3*round(pfc_losses.CX_POWER_OWL_LOSS))]);
    disp(['CX_POWER_TOP_LOSS  = ' num2str(1e-3*round(pfc_losses.CX_POWER_TOP_LOSS))]);
    disp(['CX_POWER_DOWN_LOSS = ' num2str(1e-3*round(pfc_losses.CX_POWER_DOWN_LOSS))]);
    
    
    
    
end

if PLOT_2D_LOSSES 
    if par.CALCULATE_CX
        figure
        hold on
        imagesc(pfc_losses.Rvals_cxloss,pfc_losses.Zvals_cxloss,pfc_losses.cxloss_RZ')
        plot_TOK_background_par; 
%         ylim([-0.6 0.6])
%         xlim([min(scale_X)+R0+0.02 1.21])
        
        set(gca,'fontsize',20)
        title('CU NBI-CX-losses (1MW) @wall W/m^{2} ')
        
        figure
        hold on
        imagesc(pfc_losses.Rvals_cxloss,pfc_losses.Zvals_cxloss,pfc_losses.ionloss_RZ')
        plot_TOK_background_par; 
%         ylim([-0.6 0.6])
%         xlim([min(scale_X)+R0+0.02 1.21])
        
        set(gca,'fontsize',20)
        title('CU NBI-ions-losses (1MW) @wall W/m^{2} ')
    end
    
end
    

if PLOT_NDD
    figure
    hold on
    imagesc(summary.Rvals_ndd,summary.Zvals_ndd,summary.ndd_RZ')
    plot_TOK_background_par; 
%     ylim([-0.45 0.45])
%     xlim([0.6 1.18])
    
    set(gca,'fontsize',16)
    title('Beam [1MW] => target neutrons /s/m^{3} ')
    
end
   

%%
if PLOT_CX_POWER_2D & par.CALCULATE_CX

    imagesc(summary.Rvals_cx_pow,summary.Zvals_cx_pow,summary.neutrals_power_map_2D);
    plot_TOK_background_par; ;
    set(gca,'fontsize',20)
    title('NBI neutral re-absorbed power [W/m^{3}]','fontsize',20)
%     ylim([-0.62 0.6])
%     xlim([0.55 1.22])
    
    
    figure
    set(gca,'fontsize',24)
    hold on
    grid on
    plot(summary.Rvals_cx_pow,mean(summary.CX_Ekin_power_map_2D,1),'linewidth',3);
    plot(summary.Rvals_cx_pow,mean(summary.neutrals_power_map_2D,1),'linestyle','-.','linewidth',3);
    xlabel('R [m]')
    ylabel('[W/m^{3}]')
%     xlim([0.56 1.23])
    
    legend('CX neutralization','re-ionizations')

end
%%
    
%         figure
%         hold on
%         imagesc(summary.Rvals_power,summary.Zvals_power,(summary.power_map)')
%         plot_COMPASS_vessel;
%         ylim([-0.4 0.4])
%         xlim([0.32 0.76])
%         
%         set(gca,'fontsize',16)
%         title('deposited power W/m^{3} ')
%         
% contour(scale_X+R0,scale_Z,psi_norm1_XZsmall_map',20,'g')
% contour(scale_X+R0,scale_Z,psi_norm1_XZsmall_map',[1 1],'r','linewidth',2)

figure
    imagesc(summary.Rvals_cx_pow,summary.Zvals_cx_pow,summary.NBI_power_map_2D);
    plot_TOK_background_par; ;
    set(gca,'fontsize',20)
    title('NBI deposited power [W/m^{3}]','fontsize',20)
%     ylim([-0.62 0.6])
%     xlim([0.55 1.22])


%%
figure
hold on
imagesc(pfc_losses.Rvals_cxloss,pfc_losses.Zvals_cxloss,(pfc_losses.ionloss_RZ+pfc_losses.cxloss_RZ)')
plot_TOK_background_par
ylim([-0.68 0.62])
xlim([min(scale_X)+R0+0.02 max(scale_X)+R0])
xlim([0.56 1.24])

set(gca,'fontsize',20)
title([TITLESTRING 'Fast particles [W/m^{2}]'],'fontsize',18);

colorbar


%%
view_radial_profile_from_summary;

%%