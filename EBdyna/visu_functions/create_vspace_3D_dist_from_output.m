
% clear all

function create_vspace_3D_dist_from_output(FNAME,SAVE_NEW_PREC_INPUT,PLOT_FI_DENSITY)

%% test for correct number of arguments


if nargin==0
    disp('Please specifty a filename first')
elseif nargin==1
    SAVE_NEW_PREC_INPUT=1
    PLOT_FI_DENSITY=0
end
%% Load simulation outputs and params

try
    load(FNAME)
catch
    error('invalid filename!')
end

%% init most relevant data in memory

try
    load('../../../data_common/physics_constants.mat')    
catch
    load('../data_common/physics_constants.mat')
end

%
% if strcmp(par.MACHINE_NAME,'MAST-U')
%     try
%         WALL_TOK=load('../../../data_common/MASTU_vessel_limits.mat');
%     catch
%         WALL_TOK=load('../data_common/MASTU_vessel_limits.mat');
%     end
% elseif strcmp(par.MACHINE_NAME,'COMPASS-U')
%     try
%         WALL_TOK=load('../../../data_common/CU_vessel_limits_standard_no_rescale.mat');
%     catch
%         WALL_TOK=load('../data_common/CU_vessel_limits_standard_no_rescale.mat');
%     end
% elseif strcmp(par.MACHINE_NAME,'COMPASS')
%     try
%         WALL_TOK=load('../../../data_common/COMPASS_vessel_limits.mat');
%     catch
%         WALL_TOK=load('../data_common/COMPASS_vessel_limits.mat');
%     end
% end


if par.NB_TS_RECORDED>0
    DT_SIM_RECORD=(par.time_scale(end)-par.time_scale(1));
else
    DT_SIM_RECORD=par.dt*par.NR_FUND_IN_LOOP*par.TIME_STAMP_PRECISION;
end

disp(['DT_SIM_RECORD = ' num2str(DT_SIM_RECORD)]);


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



if ~isfield(output,'birth_matrix')
    birth_time_scale=par.birth_time_scale;    
    if ~isfield('par','NB_TS_RECORDED')
        par.NB_TS_RECORDED=length(par.time_scale);
    end   
    birth_matrix=zeros(length(input.Ekin),par.NB_TS_RECORDED);    
    vR=squeeze(output.v(:,1,:));
    vR0=squeeze(input.v_ini(:,1));
    
    for TS=1:par.NB_TS_RECORDED
        birth_matrix(:,TS)=logical(vR(:,TS)~=vR0);
        %     vR0=squeeze(vR(:,TS));
    end   
    birth_matrix=birth_matrix';   
    save(par.FILENAME,'-append','birth_matrix');
else
    birth_matrix=output.birth_matrix;
end
LOSS_DT=output.time_stamp_loss>=(par.NB_STAMPS_saved-par.NB_TS_RECORDED);

if par.NB_TS_RECORDED>1
    BORN_MARKERS_SIMDT=logical(~birth_matrix(1,:)');
else
    BORN_MARKERS_SIMDT=logical(~birth_matrix(end,:)');
end

try
%     load('../../data_tokamak/physics_constants.mat');
    load('../../data_tokamak/XZsmall_fields_tokamak_pre_collapse.mat');
    load('../../data_tokamak/motions_map_dimensions.mat');
    load('../../data_tokamak/tokamak_parameters.mat');
    load('../../data_tokamak/flux_geometry.mat');
catch
%     load('./data_tokamak/physics_constants.mat');
    load('./data_tokamak/XZsmall_fields_tokamak_pre_collapse.mat');
    load('./data_tokamak/motions_map_dimensions.mat');
    load('./data_tokamak/tokamak_parameters.mat');
    load('./data_tokamak/flux_geometry.mat');
end


v2=squeeze(sum(output.v.^2,2));
output.Ekin=0.5*(input.m/eV).*v2;
output.pitch=output.vpll(:,:)./sqrt(v2(:,:));

% normalized Magnetic Moment
output.mm_norm=output.mm*Bphi0./output.Ekin;

NB_TS_AVG=par.NB_TS_RECORDED-1;
%%
warning('re-normalizing pphi values to standard convention')
% output.pphi_kin=output.pphi_kin.*(psi_global.*input.Z);
pphi_kin_overall=reshape(output.pphi_kin,1,(size(output.pphi_kin,1)).*(size(output.pphi_kin,2)));
% hist(pphi_kin_overall,linspace(-2,2,101),'r')
% hold on

% output.psi=output.pphi_kin*0;
for TS=1:size(output.pphi_kin,2)
    output.psi(:,TS)=interp2(scale_X,scale_Z,psi_XZsmall_map',squeeze(output.x(:,1,TS))-R0,squeeze(output.x(:,2,TS)),'*cubic');
end
output.psi=output.psi-psi_scale(1);

output.pphi_kin_recalc=-(input.m/eV).*squeeze(output.x(:,1,:)).*squeeze(output.v(:,3,:))+input.Z.*output.psi;
output.pphi_kin_recalc=-output.pphi_kin_recalc./(psi_global.*input.Z);

% output.pphi_kin=-(output.pphi_kin+input.Z.*output.psi); % change sign of vphi
% figure;title('old -vphi values');hist(-output.v(:,3,:),40)
% figure;title('new vphi values');hist((eV/mD).*output.pphi_kin(:,:)./squeeze(output.x(:,1,:))./input.Z+squeeze(output.v(:,3,:)),40)
%
% new convention for psi
% output.pphi_kin=(output.pphi_kin+input.Z.*(output.psi))./(psi_global.*input.Z);
% output.pphi_kin=-output.pphi_kin;
% output.pphi_kin=-(output.pphi_kin-input.Z.*(-psi_scale(1)))./psi_global;

pphi_kin_overall=reshape(output.pphi_kin_recalc,1,(size(output.pphi_kin,1)).*(size(output.pphi_kin,2)));
hist(pphi_kin_overall,linspace(-2,2,101))
% hist(pphi_kin_overall,40)

% output.pphi_kin=output.pphi_kin_recalc;

%%
disp('lost markers during overall simulation :')
disp(length(find(ejected)))

output.Ekin_end=output.Ekin(:,end);
output.mmn_end=output.mm_norm(:,end);
output.pphi_end=output.pphi_kin_recalc(:,end);
output.pitch_end=output.pitch(:,end);
output.R_end=output.x(:,1,end);
output.Z_end=output.x(:,2,end);

if strcmp(par.MACHINE_NAME,'COMPASS-U')
    Ekin_bins=(0:10:100)*1e3;
    pphi_bins=(-1.4:0.1:1.4);
elseif strcmp(par.MACHINE_NAME,'MAST-U')
    Ekin_bins=(0:10:90)*1e3;
    pphi_bins=(-1.4:0.1:1.4);
else
    Ekin_bins=(0:20:100)*1e3;
    pphi_bins=(-1.4:0.1:1.4);
end

mmn_bins=0:0.1:1.5;

Ekin_values=Ekin_bins(1:end-1)+0.5*(Ekin_bins(2)-Ekin_bins(1));
mmn_values=mmn_bins(1:end-1)+0.5*(mmn_bins(2)-mmn_bins(1));
pphi_values=pphi_bins(1:end-1)+0.5*(pphi_bins(2)-pphi_bins(1));

n_ej=logical(~ejected);
dist_vspace_3D=zeros(length(pphi_values),length(mmn_values),length(Ekin_values));

BORN_MARKERS_DT=logical(birth_matrix(end,:)');
POP_TS=logical(n_ej.*BORN_MARKERS_DT)&~output.CX_NEUTRALS;

output.POP_TS_matrix=zeros(size(output.Ekin));
output.POP_TS_matrix(:,end)=POP_TS;

% scanning through kinetic energies
for ekin_ind=1:length(Ekin_values)
    POP_EKIN_TS=POP_TS&(output.Ekin_end>Ekin_bins(ekin_ind))&(output.Ekin_end<=Ekin_bins(ekin_ind+1));
    dist_vspace=hist2d([output.pphi_end(POP_EKIN_TS) output.mmn_end(POP_EKIN_TS)],pphi_bins,mmn_bins);
    dist_vspace_3D(:,:,ekin_ind)=dist_vspace;
end

NB_IONS_SIM=sum(POP_TS);
% disp(['Number of ions in simulation at TS = ' num2str(NB_IONS_SIM)]);



% figure
% imagesc(pphi_values,mmn_values,squeeze(dist_vspace_3D(:,:,8.))');
% title(['Ekin=' num2str(Ekin_values(8)) ' eV']);
% axis xy
% xlabel('p_{\phi}')
% ylabel('\mu*B_0/E_{kin}')

%%

pphi_overall_values=output.pphi_end(POP_TS);
mmn_overall_values=output.mmn_end(POP_TS);
Ekin_overall_values=output.Ekin_end(POP_TS);
pitch_overall_values=output.pitch_end(POP_TS);
R_overall_values=output.R_end(POP_TS);
Z_overall_values=output.Z_end(POP_TS);

TS=1;
if NB_TS_AVG>0
    for TS=1:NB_TS_AVG
        output.Ekin_end=squeeze(output.Ekin(:,end-TS));
        output.mmn_end=squeeze(output.mm_norm(:,end-TS));
        output.pphi_end=squeeze(output.pphi_kin_recalc(:,end-TS));
        output.pitch_end=squeeze(output.pitch(:,end-TS));
        output.R_end=squeeze(output.x(:,1,end-TS));
        output.Z_end=squeeze(output.x(:,2,end-TS));
        
        n_ej=logical(output.time_stamp_loss>=par.NB_TIME_STAMPS-TS)|(isnan(output.time_stamp_loss));
        
        BORN_MARKERS_DT=logical(birth_matrix(end-TS,:)');
        POP_TS=logical(n_ej.*BORN_MARKERS_DT); %&squeeze(output.x_gc(:,1,end-TS))>1.1;
        output.POP_TS_matrix(:,end-TS)=POP_TS;
        
        pphi_overall_values=[pphi_overall_values ; output.pphi_end(POP_TS)];
        mmn_overall_values=[mmn_overall_values ; output.mmn_end(POP_TS)];
        Ekin_overall_values=[Ekin_overall_values ; output.Ekin_end(POP_TS)];
        pitch_overall_values=[pitch_overall_values ; output.pitch_end(POP_TS)];
        R_overall_values=[R_overall_values ; output.R_end(POP_TS)];
        Z_overall_values=[Z_overall_values ; output.Z_end(POP_TS)];

        
        dist_vspace_3D_TS=zeros(length(pphi_values),length(mmn_values),length(Ekin_values));
        for ekin_ind=1:length(Ekin_values)
            POP_EKIN_TS=POP_TS&(output.Ekin_end>Ekin_bins(ekin_ind))&(output.Ekin_end<=Ekin_bins(ekin_ind+1));
            dist_vspace_TS=hist2d([output.pphi_end(POP_EKIN_TS) output.mmn_end(POP_EKIN_TS)],pphi_bins,mmn_bins);
            dist_vspace_3D_TS(:,:,ekin_ind)=dist_vspace_TS;
        end
        
        dist_vspace_3D=dist_vspace_3D+dist_vspace_3D_TS;
        NB_IONS_SIM=NB_IONS_SIM+sum(POP_TS);
        disp(['Number of ions in simulation = ' num2str(NB_IONS_SIM)]);

    end
    DT_SIM_RECORD=(par.time_scale(end)-par.time_scale(1));

else
    DT_SIM_RECORD=par.dt*par.NR_FUND_IN_LOOP*par.TIME_STAMP_PRECISION;
end

disp(['DT_SIM_RECORD = ' num2str(DT_SIM_RECORD)]);

output.POP_TS_matrix=logical(output.POP_TS_matrix);

TOTAL_NUMBER_OF_IONS=sum(~ejected)*par.MARKER_WEIGHT;

% normalizing the vspace distribution so that its representative of the
% total number of ions generated by the NBI

dist_vspace_3D=dist_vspace_3D.*TOTAL_NUMBER_OF_IONS./sum(sum(sum(dist_vspace_3D)));

distCOM_3D=struct();
distCOM_3D.dist_vspace_3D=dist_vspace_3D;
distCOM_3D.pphi_values=pphi_values;
distCOM_3D.mmn_values=mmn_values;
distCOM_3D.Ekin_values=Ekin_values;

fast_ions.pphi=pphi_overall_values;
fast_ions.Ekin=Ekin_overall_values;
fast_ions.mm_norm=mmn_overall_values;
fast_ions.pitch=pitch_overall_values;
fast_ions.R=R_overall_values;
fast_ions.Z=Z_overall_values;
fast_ions.MARKER_WEIGHT=par.MARKER_WEIGHT

save('FI_CU5400_Rt45_dist_vspace_3D_Mads_coords.mat','distCOM_3D','fast_ions');


%%
EKIN_BIN=3

figure
imagesc(distCOM_3D.pphi_values,distCOM_3D.mmn_values,squeeze(distCOM_3D.dist_vspace_3D(:,:,EKIN_BIN))');
title(['Ekin=' num2str(Ekin_values(EKIN_BIN)) ' eV']);
axis xy
xlabel('p_{\phi}')
ylabel('\mu B_0/E_{kin}')
set(gca,'fontsize',22)


EKIN_BIN=5

figure
imagesc(distCOM_3D.pphi_values,distCOM_3D.mmn_values,squeeze(distCOM_3D.dist_vspace_3D(:,:,EKIN_BIN))');
title(['Ekin=' num2str(Ekin_values(EKIN_BIN)) ' eV']);
axis xy
xlabel('p_{\phi}')
ylabel('\mu B_0/E_{kin}')
set(gca,'fontsize',22)


EKIN_BIN=7

figure
imagesc(distCOM_3D.pphi_values,distCOM_3D.mmn_values,squeeze(distCOM_3D.dist_vspace_3D(:,:,EKIN_BIN))');
title(['Ekin=' num2str(Ekin_values(EKIN_BIN)) ' eV']);
axis xy
xlabel('p_{\phi}')
ylabel('\mu B_0/E_{kin}')
set(gca,'fontsize',22)


%%

if PLOT_FI_DENSITY
    
    Rpos_vspace=squeeze(output.x(~ejected,1,end));
    Zpos_vspace=squeeze(output.x(~ejected,2,end));
    Ekin_vspace=output.Ekin(~ejected,end);
    pitch_vspace=output.pitch(~ejected,end);
    
    Rbins = scale_X(1:25:end)+R0;
    Rbins = [Rbins(1)-(Rbins(2)-Rbins(1)) Rbins Rbins(end)+(Rbins(2)-Rbins(1))]
    Zbins = scale_Z(1:30:end);
    DR_vspace=Rbins(2)-Rbins(1);
    DZ_vspace=Zbins(2)-Zbins(1);
    
    vspace_R_scale=Rbins(1:end-1)+0.5*(DR_vspace);
    vspace_Z_scale=Zbins(1:end-1)+0.5*(DZ_vspace);
    DV_map=2*pi*repmat(vspace_R_scale, length(vspace_Z_scale),1);
    DV_map=DV_map*(Rbins(2)-Rbins(1))*(Zbins(2)-Zbins(1));
    DV_map=DV_map'; % (R,Z) convention
    
    %     Rvals = Rbins(1:end-1)+0.5*(Rbins(2)-Rbins(1));
    %     Zvals = Zbins(1:end-1)+0.5*(Zbins(2)-Zbins(1));
    pitch_bins=-1.0:0.5:1.0;
    
    Ekin_bins=(0:10:100)*1e3;
    
    DE_vspace=Ekin_bins(2)-Ekin_bins(1);
    Dpitch_vspace=pitch_bins(2)-pitch_bins(1);
    
    vspace_pitch_scale=pitch_bins(1:end-1)+0.5*(Dpitch_vspace);
    vspace_Ekin_scale=Ekin_bins(1:end-1)+0.5*(DE_vspace);
    Zbins_scale=Zbins(1:end-1)+0.5*(DZ_vspace);
    
    dens_FI_RZ=zeros(length(vspace_R_scale),length(vspace_Z_scale));
    vspace_Rmp_map=zeros(length(vspace_R_scale),length(vspace_pitch_scale),length(vspace_Ekin_scale));
    DZmp=0.04
    for Rindex=1:length(vspace_R_scale)-1
        for p_index=1:length(vspace_pitch_scale)-1
            for E_index=1:length(vspace_Ekin_scale)-1
                POP_BIN=find((Rpos_vspace>=Rbins(Rindex))&(Rpos_vspace<Rbins(Rindex+1))&(Zpos_vspace>=-DZmp)&(Zpos_vspace<=DZmp)...
                    &(Ekin_vspace>=vspace_Ekin_scale(E_index))&(Ekin_vspace<vspace_Ekin_scale(E_index+1))&(pitch_vspace>=vspace_pitch_scale(p_index))&(pitch_vspace<vspace_pitch_scale(p_index+1)));
                vspace_Rmp_map(Rindex,p_index,E_index)=length(POP_BIN);
            end
        end
    end
    for Rindex=1:length(vspace_R_scale)
        for Zindex=1:length(vspace_Z_scale)
            POP_BIN=find((Rpos_vspace>=Rbins(Rindex))&(Rpos_vspace<Rbins(Rindex+1))&(Zpos_vspace>=Zbins(Zindex))&(Zpos_vspace<Zbins(Zindex+1)));
            dens_FI_RZ(Rindex,Zindex)=length(POP_BIN);
        end
    end
    
    
    TS=1;
    disp('Building vspace R,pitch,Ekin......')
    if NB_TS_AVG>0
        for TS=1:NB_TS_AVG
            TS
            n_ej=logical(output.time_stamp_loss>=par.NB_TIME_STAMPS-TS)|(isnan(output.time_stamp_loss));
            
            BORN_MARKERS_DT=logical(birth_matrix(end-TS,:)');
            POP_TS=logical(n_ej.*BORN_MARKERS_DT);
            sum(POP_TS);
            
            Rpos_vspace=squeeze(output.x(POP_TS,1,end-TS));
            Zpos_vspace=squeeze(output.x(POP_TS,2,end-TS));
            Ekin_vspace=output.Ekin(POP_TS,end-TS);
            pitch_vspace=output.pitch(POP_TS,end-TS);
            
            for Rindex=1:length(vspace_R_scale)-1
                for p_index=1:length(vspace_pitch_scale)-1
                    for E_index=1:length(vspace_Ekin_scale)-1
                        POP_BIN=find((Rpos_vspace>=Rbins(Rindex))&(Rpos_vspace<Rbins(Rindex+1))&(Zpos_vspace>=-DZmp)&(Zpos_vspace<=DZmp)...
                            &(Ekin_vspace>=Ekin_bins(E_index))&(Ekin_vspace<Ekin_bins(E_index+1))&(pitch_vspace>=pitch_bins(p_index))&(pitch_vspace<pitch_bins(p_index+1)));
                        vspace_Rmp_map(Rindex,p_index,E_index)=vspace_Rmp_map(Rindex,p_index,E_index)+length(POP_BIN);
                    end
                end
            end
            for Rindex=1:length(vspace_R_scale)
                for Zindex=1:length(vspace_Z_scale)
                    POP_BIN=find((Rpos_vspace>=Rbins(Rindex))&(Rpos_vspace<Rbins(Rindex+1))&(Zpos_vspace>=Zbins(Zindex))&(Zpos_vspace<Zbins(Zindex+1)));
                    dens_FI_RZ(Rindex,Zindex)=dens_FI_RZ(Rindex,Zindex)+length(POP_BIN);
                end
            end
        end
        DT_SIM_RECORD=(par.time_scale(end)-par.time_scale(1));
        
    else
        DT_SIM_RECORD=par.dt*par.NR_FUND_IN_LOOP*par.TIME_STAMP_PRECISION;
    end
    
    vspace_Rmp_map=par.MARKER_WEIGHT*vspace_Rmp_map/(NB_TS_AVG+1);
    dens_FI_RZ=par.MARKER_WEIGHT*dens_FI_RZ/(NB_TS_AVG+1);
    dens_FI_RZ=dens_FI_RZ./DV_map;
    
    
    
    
    
    
    %%
    
    hf=figure;
    set(gca,'fontsize',22)
    hold on
    grid on
    pcolor(vspace_R_scale,vspace_Z_scale,dens_FI_RZ'); shading interp;
    plot_TOK_background_par;
    xlabel('R [m]')
    ylabel('Z [m]')
    
    title(['CU#5400 density FI [m^{-3}]']);
    ylim([min(vspace_Z_scale)+0.5*(vspace_Z_scale(2)-vspace_Z_scale(1)) max(vspace_Z_scale)-0.5*(vspace_Z_scale(2)-vspace_Z_scale(1))])
    xlim([min(vspace_R_scale)+0.5*(vspace_R_scale(2)-vspace_R_scale(1)) max(vspace_R_scale)-0.5*(vspace_R_scale(2)-vspace_R_scale(1))])
    colorbar
    
    
    xlim([0.5 1.3])
    ylim([-0.7 0.7])
    set(gca,'fontsize',16)
end

%%

disp ('saving new input prec file for further orbit categorization!')

if SAVE_NEW_PREC_INPUT
    new_input=struct();
    old_input=input;
    
    new_input.m=old_input.m;
    new_input.Z=old_input.Z;
    POP_TS=logical(output.POP_TS_matrix(:,end));
    new_input.x=squeeze(output.x(POP_TS,:,end));
    new_input.v=squeeze(output.v(POP_TS,:,end));
    
    for TS=1:NB_TS_AVG
        POP_TS=logical(output.POP_TS_matrix(:,end-TS));
        new_input.x=[new_input.x ; squeeze(output.x(POP_TS,:,end-TS))];
        new_input.v=[new_input.v ; squeeze(output.v(POP_TS,:,end-TS))];
    end
    
    new_input.N_total=size(new_input.x,1)
    if exist('folder_save_string')
        FILENAME_PREC=[folder_save_string par.FILENAME(1:end-4) '_ini_prec.mat'];
    else
        FILENAME_PREC=['./' par.FILENAME(1:end-4) '_ini_prec.mat'];
    end
    
    input=new_input;
    save(FILENAME_PREC,'input')
    input=old_input;
end

