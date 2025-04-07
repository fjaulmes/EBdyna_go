
% clear all

function summary_vspace_Ekin_pitch_CX(FNAME,SAVE_NEW_PREC_INPUT,TITLESTRING,ISALPHADIST)

try
    load('../../../data_common/physics_constants.mat')    
catch
    load('../data_common/physics_constants.mat')
end

if nargin==0
    disp('Please specifty a filename first')
end
try
    load(FNAME)
catch
    error('invalid filename!')
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

if nargin<=1
    SAVE_NEW_PREC_INPUT=0;
    disp('No post-sim input file generated.')
end
if nargin<=2
    TITLESTRING='';
end
if nargin<=3
    ISALPHADIST=0;
end


if par.NB_TS_RECORDED>0
    DT_SIM_RECORD=(par.time_scale(end)-par.time_scale(1));
else
    DT_SIM_RECORD=par.dt*par.NR_FUND_IN_LOOP*par.TIME_STAMP_PRECISION;
end

disp(['DT_SIM_RECORD = ' num2str(DT_SIM_RECORD)]);


PLOT_SD=0;
PLOT_NEUTRALS_RELOC=0;
PLOT_CX_POWER_2D=0;
PLOT_DOTS=0;
PLOT_DOTS_WALL=1;
PLOT_HIST_LOC=0;
PLOT_HIST_IONS_NEUTRALS=0;
PLOT_CXLOSS=1;
PLOT_IONLOSS=1;
PLOT_NDD=1;
PLOT_NDD_DETAILS=1;

PLOT_HIST_CX_SUMMARY=0;

CALC_3D_VSPACE=1;

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






% eV=1.6e-19; % recalculated after this
output.Ekin=0.5*(input.m/eV).*squeeze(sum(output.v.^2,2));

ERROR_MARKERS=[];
ERROR_MARKERS=((abs(imag(output.pphi_kin(:,end)))+abs(imag(output.Ekin(:,end))))>0)|...
                isnan(output.Ekin(:,end));
%     |(output.Ekin(:,end)>3.6e6)|(abs(output.Etot_i(:,1))>3.6e6)|output.Pdep_i(:,1)<0;
ERROR_MARKERS=find(ERROR_MARKERS);
if ~isempty('ERROR_MARKERS')
    disp('removing some error MARKERS based on pphi_kin')
    disp([num2str(length(ERROR_MARKERS)) ' markers removed!']);
    output.x(ERROR_MARKERS,:,:)=[];
    output.v(ERROR_MARKERS,:,:)=[];
%     output.x_ej_prev3(ERROR_MARKERS,:)=[];
%     output.x_ej_prev2(ERROR_MARKERS,:)=[];
%     output.x_ej_prev(ERROR_MARKERS,:)=[];
    output.x_ej(ERROR_MARKERS,:)=[];
    output.x_ej_next(ERROR_MARKERS,:)=[];
    output.Ekin_ej(ERROR_MARKERS)=[];
    output.vpll_ej(ERROR_MARKERS)=[];
    output.ejected_wall(ERROR_MARKERS)=[];
    output.Pdep_e(ERROR_MARKERS,:)=[];
    output.Pdep_i(ERROR_MARKERS,:)=[];
    output.Etot_e(ERROR_MARKERS,:)=[];
    output.Etot_i(ERROR_MARKERS,:)=[];
    output.delta_Ekin(ERROR_MARKERS,:)=[];
    output.cum_tor_ang_mom(ERROR_MARKERS,:)=[];
    output.Edep_th(ERROR_MARKERS)=[];
    output.time_step_loss(ERROR_MARKERS)=[];
    output.ejected_sd(ERROR_MARKERS)=[];
    output.psi_value_avg(ERROR_MARKERS,:)=[];
    output.x_gc(ERROR_MARKERS,:,:)=[];
    output.vpll(ERROR_MARKERS,:)=[];
    output.pphi_kin(ERROR_MARKERS,:)=[];
    output.Delta_pphi(ERROR_MARKERS)=[];
    output.mm(ERROR_MARKERS,:)=[];
    output.dpphi_dt(ERROR_MARKERS)=[];
    output.time_stamp_loss(ERROR_MARKERS)=[];
    output.Ekin(ERROR_MARKERS,:)=[];
    output.birth_matrix(:,ERROR_MARKERS)=[];
    
    input.Ekin(ERROR_MARKERS)=[];
    input.particle_nr(ERROR_MARKERS)=[];
    input.mm(ERROR_MARKERS)=[];
    input.pphi_kin(ERROR_MARKERS)=[];
    input.x_ini(ERROR_MARKERS,:)=[];
    input.v_ini(ERROR_MARKERS,:)=[];
    input.vpll_ini(ERROR_MARKERS)=[];
    input.x_gc_end(ERROR_MARKERS,:)=[];
    input.x_end(ERROR_MARKERS,:)=[];
    input.v_end(ERROR_MARKERS,:)=[];
    
    ejected(ERROR_MARKERS)=[];
    
    save(FNAME,'-append','ejected','output','input')
    
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

    BUGFIX=0
    
%     par.MARKER_WEIGHT=par.INPUT_POWER./(sum(input.Ekin)*eV/par.time_scale(end));

if ~isfield(par,'MARKER_WEIGHT_OLD')|BUGFIX
    warning('recalculating par.MARKER_WEIGHT')
    par.MARKER_WEIGHT_OLD=par.MARKER_WEIGHT;
%     par.MARKER_WEIGHT=par.INPUT_POWER./(sum(input.Ekin)*eV/par.time_scale(end));
%     par.MARKER_WEIGHT=par.INPUT_POWER./(sum(input.Ekin(BORN_MARKERS_SIMDT))*eV/DT_SIM_RECORD);
    par.MARKER_WEIGHT=0.5.*(par.MARKER_WEIGHT_OLD+par.INPUT_POWER./(sum(input.Ekin(BORN_MARKERS_SIMDT))*eV/DT_SIM_RECORD));
    
    % if isfield(input,'weight')
    if BUGFIX
    warning('adjusting outputs according to recalculated overall weight...')
    for TS=1:size(output.Pdep_e,2)
%         output.Pdep_e(:,TS)=par.NR_FUND_IN_LOOP.*par.dt.*output.Pdep_e(:,TS).*par.MARKER_WEIGHT./par.MARKER_WEIGHT_OLD;
%         output.Pdep_i(:,TS)=par.NR_FUND_IN_LOOP.*par.dt.*output.Pdep_i(:,TS).*par.MARKER_WEIGHT./par.MARKER_WEIGHT_OLD;
%         output.Etot_e(:,TS)=par.NR_FUND_IN_LOOP.*par.dt.*output.Etot_e(:,TS).*par.MARKER_WEIGHT./par.MARKER_WEIGHT_OLD;
%         output.Etot_i(:,TS)=par.NR_FUND_IN_LOOP.*par.dt.*output.Etot_i(:,TS).*par.MARKER_WEIGHT./par.MARKER_WEIGHT_OLD;
%         if isfield(output,'Pdep_imp')
%             output.Pdep_imp(:,TS)=par.NR_FUND_IN_LOOP.*par.dt.*output.Pdep_imp(:,TS).*par.MARKER_WEIGHT./par.MARKER_WEIGHT_OLD;
%             output.Etot_imp(:,TS)=par.NR_FUND_IN_LOOP.*par.dt.*output.Etot_imp(:,TS).*par.MARKER_WEIGHT./par.MARKER_WEIGHT_OLD;
%         end
        
        output.Pdep_e(:,TS)=output.Pdep_e(:,TS).*par.MARKER_WEIGHT./par.MARKER_WEIGHT_OLD;
        output.Pdep_i(:,TS)=output.Pdep_i(:,TS).*par.MARKER_WEIGHT./par.MARKER_WEIGHT_OLD;
        output.Etot_e(:,TS)=output.Etot_e(:,TS).*par.MARKER_WEIGHT./par.MARKER_WEIGHT_OLD;
        output.Etot_i(:,TS)=output.Etot_i(:,TS).*par.MARKER_WEIGHT./par.MARKER_WEIGHT_OLD;
        if isfield(output,'Pdep_imp')
            output.Pdep_imp(:,TS)=output.Pdep_imp(:,TS).*par.MARKER_WEIGHT./par.MARKER_WEIGHT_OLD;
            output.Etot_imp(:,TS)=output.Etot_imp(:,TS).*par.MARKER_WEIGHT./par.MARKER_WEIGHT_OLD;
        end
        
        output.ndd_val(:,TS)=output.ndd_val(:,TS).*par.MARKER_WEIGHT./par.MARKER_WEIGHT_OLD;
        if isfield(output,'ndd_val_BB')
            output.ndd_val_BB(:,TS)=output.ndd_val_BB(:,TS).*(par.MARKER_WEIGHT./par.MARKER_WEIGHT_OLD).^2;
            output.ndd_val_BP(:,TS)=output.ndd_val_BP(:,TS).*par.MARKER_WEIGHT./par.MARKER_WEIGHT_OLD;
        end
    end
    output.Edep_th=output.Edep_th.*par.MARKER_WEIGHT./par.MARKER_WEIGHT_OLD;
    output.Ekin_ej=output.Ekin_ej.*par.MARKER_WEIGHT./par.MARKER_WEIGHT_OLD;
    end
end    
    WEIGHTFIX=0;
    
    % if isfield(input,'weight')
    if WEIGHTFIX
        par.MARKER_WEIGHT=par.INPUT_POWER./(sum(input.Ekin(BORN_MARKERS_SIMDT))*eV/DT_SIM_RECORD);
        if ~isfield(input,'weight')
            input.weight=input.Ekin*0+par.MARKER_WEIGHT;
        else
            %         par.MARKER_WEIGHT=mean(input.weight);
        end
        warning('adjusting outputs according to input.weight...')
        for TS=1:size(output.Pdep_e,2)
            output.Pdep_e(:,TS)=output.Pdep_e(:,TS).*input.weight./mean(input.weight);
            output.Pdep_i(:,TS)=output.Pdep_i(:,TS).*input.weight./mean(input.weight);
            output.Etot_e(:,TS)=output.Etot_e(:,TS).*input.weight./mean(input.weight);
            output.Etot_i(:,TS)=output.Etot_i(:,TS).*input.weight./mean(input.weight);
            if isfield(output,'Pdep_imp')
                output.Pdep_imp(:,TS)=output.Pdep_imp(:,TS).*input.weight./mean(input.weight);
                output.Etot_imp(:,TS)=output.Etot_imp(:,TS).*input.weight./mean(input.weight);
            end
            output.ndd_val(:,TS)=output.ndd_val(:,TS).*input.weight./mean(input.weight);
            if isfield(output,'ndd_val_BB')
                output.ndd_val_BB(:,TS)=output.ndd_val_BB(:,TS).*input.weight./mean(input.weight);;
                output.ndd_val_BP(:,TS)=output.ndd_val_BP(:,TS).*input.weight./mean(input.weight);;
            end
        end
        output.Edep_th=output.Edep_th.*input.weight./mean(input.weight);;
        output.Ekin_ej=output.Ekin_ej.*input.weight./mean(input.weight);;
    end



par.NB_TS_RECORDED=size(output.Pdep_e,2)
% end
if ~isfield(output,'ndd_val_BB')
    output.ndd_val_BB=output.ndd_val_thermal*0;
end
if ~isfield(output,'ndd_val_BP')
    output.ndd_val_BP=output.ndd_val_thermal*0;
end
if ~isfield(output,'ndd_val_th')
    output.ndd_val_th=output.ndd_val_thermal*0;
end
% output.Pdep_e=output.Pdep_e.*par.MARKER_WEIGHT./par.MARKER_WEIGHT_OLD;
% output.Pdep_i=output.Pdep_i.*par.MARKER_WEIGHT./par.MARKER_WEIGHT_OLD;
% output.Pdep_imp=output.Pdep_imp.*par.MARKER_WEIGHT./par.MARKER_WEIGHT_OLD;
% output.Etot_e=output.Etot_e.*par.MARKER_WEIGHT./par.MARKER_WEIGHT_OLD;
% output.Etot_i=output.Etot_i.*par.MARKER_WEIGHT./par.MARKER_WEIGHT_OLD;
if ~WEIGHTFIX
save(FNAME,'-append','par','output')
end

if ~isfield(output,'CX_NEUTRALS')
    output.CX_NEUTRALS=ejected*0;
    output.ejected_CX=ejected*0;
    warning('No CX data for this simulation!');
end


output.ejected_CX(isnan(output.ejected_CX))=0;
output.ejected_sd(isnan(output.ejected_sd))=0;
output.ejected_wall(isnan(output.ejected_wall))=0;


try
    load('../../data_tokamak/physics_constants.mat');
    load('../../data_tokamak/XZsmall_fields_tokamak_pre_collapse.mat');
    load('../../data_tokamak/motions_map_dimensions.mat');
catch
    load('./data_tokamak/physics_constants.mat');
    load('./data_tokamak/XZsmall_fields_tokamak_pre_collapse.mat');
    load('./data_tokamak/motions_map_dimensions.mat');
end

% if isfield(output,'v')
%     v2=squeeze(sum(output.v.^2,2));
% end 
output.Ekin=0.5*(input.m/eV).*squeeze(sum(output.v.^2,2));

output.delta_Ekin(isnan(output.delta_Ekin))=0;
output.Pdep_e(isnan(output.Etot_e))=0;
output.Pdep_i(isnan(output.Etot_i))=0;

output.Ekin_ej(isnan(output.Ekin_ej))=0;
output.v(isnan(output.v))=0;
output.Ekin(isnan(output.Ekin))=0;

v2=squeeze(sum(output.v.^2,2));

if ~isfield(output,'Ekin')
    
    ejected=logical(ejected);
    
    % moved to combine_output_eq
    % output.time_stamp_loss=ceil(output.time_step_loss.*par.NB_TIME_STAMPS/par.NB_TIME_STEPS);
    %
    % save(par.FILENAME,'-append','output');
    
end
%security on max Ekin for alphas
MAX_EKIN=4e6;
ejected(output.Ekin(:,end)>MAX_EKIN)=1;
output.Ekin(output.Ekin>MAX_EKIN)=0;
% output.Ekin(isnan(output.Ekin(:,end)),:)=0;

if ~isfield(output,'CX_flag')
    output.CX_flag=output.Ekin*0;
end
if ~isfield(output,'x_CXn')
    output.x_CXn=output.x*0;
end
if ~isfield(output,'ejected_CX')
    output.ejected_CX=ejected*0;
end
output.ejected_CX=logical(output.ejected_CX);



NB_TS_AVG=par.NB_TS_RECORDED-1;

%%
disp('lost markers during overall simulation :')
disp(length(find(ejected)))

output.Ekin_end=output.Ekin(:,end);
output.pitch_end=output.vpll(:,end)./sqrt(v2(:,end));

if ISALPHADIST
    Ekin_bins=(0:0.2:3.6)*1e6;
else
    if strcmp(par.MACHINE_NAME,'COMPASS-U')
        Ekin_bins=(0:4:100)*1e3;
    elseif strcmp(par.MACHINE_NAME,'MAST-U')
        Ekin_bins=(0:4:90)*1e3;
    else
        Ekin_bins=(0:4:90)*1e3;
    end
end
pitch_bins=-1.00:0.04:1.00;
Ekin_values=Ekin_bins(1:end-1)+0.5*(Ekin_bins(2)-Ekin_bins(1));
pitch_values=pitch_bins(1:end-1)+0.5*(pitch_bins(2)-pitch_bins(1));

n_ej=logical(~ejected);

BORN_MARKERS_DT=logical(birth_matrix(end,:)');
POP_TS=logical(n_ej.*BORN_MARKERS_DT)&~output.CX_NEUTRALS; %&squeeze(output.x_gc(:,1,1))>1.1;
dist_vspace=hist2d([output.Ekin_end(POP_TS) output.pitch_end(POP_TS)],Ekin_bins,pitch_bins);

%%
TS=1;
if NB_TS_AVG>0
    for TS=1:NB_TS_AVG
        output.Ekin_end=squeeze(output.Ekin(:,end-TS));
        output.pitch_end=output.vpll(:,end-TS)./sqrt(v2(:,end-TS));
        
        n_ej=logical(output.time_stamp_loss>=par.NB_TIME_STAMPS-TS)|(isnan(output.time_stamp_loss));
        
        BORN_MARKERS_DT=logical(birth_matrix(end-TS,:)');
        POP_TS=logical(n_ej.*BORN_MARKERS_DT); %&squeeze(output.x_gc(:,1,end-TS))>1.1;
        dist_vspace_TS=hist2d([output.Ekin_end(POP_TS) output.pitch_end(POP_TS)],Ekin_bins,pitch_bins);
        
        dist_vspace=dist_vspace+dist_vspace_TS;
    end
    DT_SIM_RECORD=(par.time_scale(end)-par.time_scale(1));

else
    DT_SIM_RECORD=par.dt*par.NR_FUND_IN_LOOP*par.TIME_STAMP_PRECISION;
end

disp(['DT_SIM_RECORD = ' num2str(DT_SIM_RECORD)]);




%%

summary=struct();

if CALC_3D_VSPACE
    
    output.pitch=output.vpll(:,:)./sqrt(v2(:,:));
    
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
    if ISALPHADIST
        Ekin_bins=(0:0.2:3.6)*1e6;
    else
        Ekin_bins=(0:10:100)*1e3;
    end

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
    
    disp(['DT_SIM_RECORD = ' num2str(DT_SIM_RECORD)]);
    
    summary.vspace_pitch_scale=vspace_pitch_scale;
    summary.vspace_Ekin_scale=vspace_Ekin_scale;
    summary.vspace_R_scale=vspace_R_scale;
    summary.vspace_Z_scale=vspace_Z_scale;
    summary.vspace_Rmp_map=vspace_Rmp_map;
    summary.vspace_dens_FI_RZp=dens_FI_RZ;
    

    
    
    hf=figure;
    set(gca,'fontsize',22)
    hold on
    grid on
    pcolor(vspace_R_scale,vspace_Z_scale,dens_FI_RZ'); shading interp; 
    plot_TOK_background_par;
    xlabel('R [m]')
    ylabel('Z [m]')
    
    title([TITLESTRING ' density FI [m^{-3}]']);
    ylim([min(summary.vspace_Z_scale)+0.5*(vspace_Z_scale(2)-vspace_Z_scale(1)) max(summary.vspace_Z_scale)-0.5*(vspace_Z_scale(2)-vspace_Z_scale(1))])
    xlim([min(summary.vspace_R_scale)+0.5*(vspace_R_scale(2)-vspace_R_scale(1)) max(summary.vspace_R_scale)-0.5*(vspace_R_scale(2)-vspace_R_scale(1))])
    colorbar

        
    if PLOT_SD
        vspace_Rmp_map_ekin1=squeeze(vspace_Rmp_map(:,:,1));
       hf=figure;
        set(gca,'fontsize',22)
        hold on
        grid on
        surf(vspace_R_scale,vspace_pitch_scale,vspace_Rmp_map_ekin1'); shading interp; view(2)
        xlabel('R [m]')
        ylabel('pitch')
        
        title([TITLESTRING 'Ekin~7.5 keV']);
        
    end
end



%%
Ekin_array=repmat(Ekin_values,1,length(pitch_values));
pitch_array=repmat(pitch_values,length(Ekin_values),1);
% scatter(Ekin_array(:),pitch_array(:),200,dist_vspace(:),'filled')
% hf=figure;imagesc(Ekin_values,pitch_values,dist_vspace');
% 
dist_norm=sum(sum(dist_vspace))*(Ekin_bins(2)-Ekin_bins(1))*(pitch_bins(2)-pitch_bins(1));
dist_fabien=dist_vspace/dist_norm;

% dist_fabien=min(dist_fabien,1.55e-5);

% dist_fabien(end,end)=1.9e-5;

hf=figure;
set(gca,'fontsize',22)
hold on
grid on
surf(Ekin_values,pitch_values,dist_fabien'); shading interp; view(2)
xlabel('Ekin [eV]')
ylabel('pitch')

title([TITLESTRING 'EBdyna'])

summary.Ekin_values=Ekin_values;
summary.pitch_values=pitch_values;
summary.dist_vspace=dist_fabien;



%%

summary.dist_Ekin=sum(dist_fabien,2);

if PLOT_SD
    
    figure(11);
    %
    
    set(gca,'fontsize',20)
    hold on
    grid on
    plot(Ekin_values, summary.dist_Ekin,'linewidth',2)
    xlabel('Ekin [eV]')
end



%%
disp('statistics for the total number of time stamps recorded in this file');

output.B_ej=interp2(scale_X+R0,scale_Z,Btot_XZ_map',output.x_ej(:,1),output.x_ej(:,2));
% output.Ekin_ej=0.5*(input.m/eV).*output.vpll_ej(:,end).^2+output.mm_ej.*output.B_ej;

try
    output.ejected_sd_dt=logical(output.ejected_sd)&output.time_stamp_loss>(par.NB_STAMPS_saved-par.NB_TS_RECORDED);
catch
    output.ejected_sd_dt=logical(ejected*0);
end

try
%     output.ejected_CX_dt=logical(squeeze(output.ejected_CX))&~output.ejected_sd&output.time_stamp_loss>par.NB_STAMPS_saved-par.NB_TS_RECORDED;
%     output.did_CX_dt=logical(squeeze(output.CX_flag(:,end)>0))&~output.ejected_sd&output.time_stamp_loss>(par.NB_STAMPS_saved-par.NB_TS_RECORDED);
    output.did_CX_dt=BORN_MARKERS_SIMDT&logical(squeeze(output.CX_flag(:,end)-output.CX_flag(:,1))>0);
%     output.ejected_CX_dt=logical(squeeze(output.ejected_CX))&output.did_CX_dt;
    output.ejected_CX_dt=(output.time_stamp_loss>(par.NB_STAMPS_saved-par.NB_TS_RECORDED)&output.ejected_CX);
catch
    output.ejected_CX_dt=logical(ejected*0);
end

% try
%     output.did_CX_dt=logical(squeeze(output.CX_flag(:,end)>0))&~output.ejected_sd&output.time_stamp_loss>par.NB_STAMPS_saved-par.NB_TS_RECORDED;
%     ISNAN=output.time_stamp_loss<par.NB_STAMPS_saved-par.NB_TS_RECORDED;
%     output.did_CX_dt=~ISNAN&logical(squeeze(output.CX_flag(:,end)-output.CX_flag(:,1)>0));
% %     for ts=2:par.NB_STAMPS_saved
% %         ISNAN=ISNAN&isnan(output.did_CX_dt);
% %         output.did_CX_dt=logical(squeeze(output.CX_flag(:,ts)-output.CX_flag(:,1)>0))|output.did_CX_dt;
% %      end
% catch
%     output.did_CX_dt=logical(ejected*0);
% end

if ~isfield(output,'ejected_CX')
    output.ejected_CX=logical(ejected*0);
end
% output.ejected_wall=ejected&~output.ejected_sd&~logical(squeeze(output.ejected_CX(:,1)));
if ~isfield(output,'ejected_wall')
    output.ejected_wall=ejected&~output.ejected_sd;
end

if ~isfield(output,'ions_ejected_wall')
    output.ions_ejected_wall=output.ejected_wall&~output.ejected_CX;
end

output.ejected_wall_dt=output.ejected_wall&output.time_stamp_loss>(par.NB_STAMPS_saved-par.NB_TS_RECORDED);
output.ions_ejected_wall_dt=output.ions_ejected_wall&output.time_stamp_loss>(par.NB_STAMPS_saved-par.NB_TS_RECORDED);

summary.DT_SS=par.time_scale(end)-par.time_scale(1);
summary.ejected_wall_dt=output.ejected_wall_dt;
summary.ion_ejected_wall_dt=output.ions_ejected_wall_dt;
summary.did_CX_dt=output.did_CX_dt;
summary.ejected_CX_dt=output.ejected_CX_dt;
summary.ejected_sd_dt=output.ejected_sd_dt;

summary.Pdep_e_evol=zeros(1,par.NB_TS_RECORDED);
summary.Pdep_i_evol=zeros(1,par.NB_TS_RECORDED);
summary.Pdep_imp_evol=zeros(1,par.NB_TS_RECORDED);
summary.Ploss_CX_evol=zeros(1,par.NB_TS_RECORDED);
summary.Ploss_ion_evol=zeros(1,par.NB_TS_RECORDED);
summary.Ploss_th_evol=zeros(1,par.NB_TS_RECORDED);

DT_stamp=par.time_scale(2)-par.time_scale(1);

summary.Pdep_e_evol(end)=sum(output.Pdep_e(~ejected,end));
summary.Pdep_i_evol(end)=sum(output.Pdep_i(~ejected,end));
if par.IMPURITY_ZAVG>0
    summary.Pdep_imp_evol(end)=sum(output.Pdep_imp(~ejected,end));
end

summary.Ploss_CX_evol(end)=eV.*par.MARKER_WEIGHT*sum(output.Ekin_ej(output.ejected_CX &output.time_stamp_loss==par.NB_TIME_STAMPS))/DT_stamp;
summary.Ploss_ion_evol(end)=eV.*par.MARKER_WEIGHT*sum(output.Ekin_ej(output.ejected_wall&~output.ejected_CX &output.time_stamp_loss==par.NB_TIME_STAMPS))/DT_stamp;
summary.Ploss_th_evol(end)=eV.*par.MARKER_WEIGHT*sum(output.Ekin_ej(output.ejected_sd &output.time_stamp_loss==par.NB_TIME_STAMPS))/DT_stamp;

TS=1;
if NB_TS_AVG>0
    for TS=1:NB_TS_AVG
        n_ej=logical(output.time_stamp_loss>=par.NB_TIME_STAMPS-TS)|(isnan(output.time_stamp_loss));
        
        BORN_MARKERS_DT=logical(birth_matrix(end-TS,:)');
        POP_TS=logical(n_ej.*BORN_MARKERS_DT);
        
        summary.Pdep_e_evol(end-TS)=sum(output.Pdep_e(POP_TS,end-TS));
        summary.Pdep_i_evol(end-TS)=sum(output.Pdep_i(POP_TS,end-TS));
        if par.IMPURITY_ZAVG>0
            summary.Pdep_imp_evol(end-TS)=sum(output.Pdep_imp(POP_TS,end-TS));
        end
        summary.Ploss_CX_evol(end-TS)=eV.*par.MARKER_WEIGHT*sum(output.Ekin_ej(output.ejected_CX &output.time_stamp_loss==par.NB_TIME_STAMPS-TS))/DT_stamp;
        summary.Ploss_ion_evol(end-TS)=eV.*par.MARKER_WEIGHT*sum(output.Ekin_ej(output.ejected_wall&~output.ejected_CX &output.time_stamp_loss==par.NB_TIME_STAMPS-TS))/DT_stamp;
        summary.Ploss_th_evol(end-TS)=eV.*par.MARKER_WEIGHT*sum(output.Ekin_ej(output.ejected_sd &output.time_stamp_loss==par.NB_TIME_STAMPS-TS))/DT_stamp;
    end
    DT_SIM_RECORD=(par.time_scale(end)-par.time_scale(1));
    
else
    DT_SIM_RECORD=par.dt*par.NR_FUND_IN_LOOP*par.TIME_STAMP_PRECISION;
end

output.ejected_wallCX_dt=logical(output.ejected_wall|output.ejected_CX)&(output.time_stamp_loss>(par.NB_STAMPS_saved-par.NB_TS_RECORDED));
summary.Power_sim_evol=summary.Pdep_e_evol+summary.Pdep_i_evol+summary.Pdep_imp_evol+...
    summary.Ploss_CX_evol+...
    summary.Ploss_ion_evol + summary.Ploss_th_evol ;

summary.TOTAL_WALL_POWER_LOSS=mean(summary.Ploss_ion_evol+summary.Ploss_CX_evol);
summary.POWER_ION_WALL_LOSS=mean(summary.Ploss_ion_evol);
summary.POWER_CX_LOSS=mean(summary.Ploss_CX_evol);

disp('summary.Power_sim_evol = ');
disp(summary.Power_sim_evol)

%%
% if exist('WALL_TOK')
%     R_OWL_UP=WALL_TOK.PFC1_OWL(1,1);
%     Z_OWL_UP=WALL_TOK.PFC1_OWL(1,2);
%     R_OWL_DOWN=WALL_TOK.PFC1_OWL(end,1);
%     Z_OWL_DOWN=WALL_TOK.PFC1_OWL(end,2);
%     
%     disp(['OWL limiter top  (R,Z) coordinates : ' num2str(R_OWL_UP) ' , ' num2str(Z_OWL_UP)])
%     disp(['OWL limiter down (R,Z) coordinates : ' num2str(R_OWL_DOWN) ' , ' num2str(Z_OWL_DOWN)])
% else
    R_OWL_UP=WALL_TOK.R_outer_vessel(end,1);
    Z_OWL_UP=max(WALL_TOK.Z_outer_vessel);
    R_OWL_DOWN=WALL_TOK.R_outer_vessel(1,1);
    Z_OWL_DOWN=min(WALL_TOK.Z_outer_vessel);
    
    disp(['OWL limiter top  (R,Z) coordinates : ' num2str(R_OWL_UP) ' , ' num2str(Z_OWL_UP)])
    disp(['OWL limiter down (R,Z) coordinates : ' num2str(R_OWL_DOWN) ' , ' num2str(Z_OWL_DOWN)])
    WALL_TOK;
% end
% else
%     load('./deta_tokamak/CU_vessel_limits.mat')
%     R_OWL_UP=R_OWL_(1,1);
%     Z_OWL_UP=wall_CU.PFC1_OWL(1,2);
%     R_OWL_DOWN=wall_CU.PFC1_OWL(end,1);
%     Z_OWL_DOWN=wall_CU.PFC1_OWL(end,2);
%     
%     disp(['OWL limiter top  (R,Z) coordinates : ' num2str(R_OWL_UP) ' , ' num2str(Z_OWL_UP)])
%     disp(['OWL limiter down (R,Z) coordinates : ' num2str(R_OWL_DOWN) ' , ' num2str(Z_OWL_DOWN)])
% 
% end
if PLOT_DOTS_WALL
    
    
    figure
    hold on
    disp('ejected ion markers to wall : [whole sim, not from slowing down or CX]')
    disp(length(find(ejected&~output.ejected_sd&~output.ejected_CX_dt)))
    try
        plot(output.x_ej(output.ions_ejected_wall_dt,1),output.x_ej(output.ions_ejected_wall_dt,2),'g.','markersize',8)
        if par.CALCULATE_CX
            plot(output.x_ej(output.ejected_CX_dt,1),output.x_ej(find(output.ejected_CX_dt),2),'m.','markersize',6)
        end
        legend('ion wall loss','CX neutral wall loss')
    end
    
    disp('Total number of NBI ions lost to wall:')
    sum(output.ions_ejected_wall_dt)*par.MARKER_WEIGHT/summary.DT_SS
    
    plot_TOK_background_par;
    
    plot(output.x_ej(output.ions_ejected_wall_dt,1),output.x_ej(output.ions_ejected_wall_dt,2),'g.','markersize',8)
    if par.CALCULATE_CX
        plot(output.x_ej(output.ejected_CX_dt,1),output.x_ej(find(output.ejected_CX_dt),2),'m.','markersize',6)
    end
    
end
%%
if PLOT_DOTS
    
    
    figure
    hold on
    disp('ejected markers to wall : [not from slowing down or CX]')
    disp(length(find(ejected&~output.ejected_sd&~output.ejected_CX_dt)))
    try
        plot(output.x_CXn(output.CX_NEUTRALS&output.CX_flag(:,end)>=1,1),output.x_CXn(output.CX_NEUTRALS&output.CX_flag(:,end)>=1,2),'color',[0 0 1],'linestyle','none','Marker','.','markersize',8)
        plot(output.x_ej(output.ejected_wall_dt&~output.ejected_CX_dt,1),output.x_ej(output.ejected_wall_dt&~output.ejected_CX_dt,2),'g.','markersize',15)
        plot(output.x_ej(output.ejected_CX_dt,1),output.x_ej(find(output.ejected_CX_dt),2),'m.','markersize',15)
        legend('neutralization','ion wall loss','CX neutral wall loss')
    end

    plot_TOK_background_par;
    plot(output.x_ej(output.ejected_wall_dt&~output.ejected_CX_dt,1),output.x_ej(output.ejected_wall_dt&~output.ejected_CX_dt,2),'g.','markersize',12)
    plot(output.x_ej(output.ejected_CX_dt,1),output.x_ej(find(output.ejected_CX_dt),2),'m.','markersize',12)
    
    % plot(output.x_CXn(output.CX_NEUTRALS&output.CX_flag(:,end)==1,1),output.x_CXn(output.CX_NEUTRALS&output.CX_flag(:,end)==1,2),'color',[0 0 1],'linestyle','none','Marker','.','markersize',10)
    plot(output.x_CXn(output.CX_NEUTRALS&output.CX_flag(:,end)==1,1),output.x_CXn(output.CX_NEUTRALS&output.CX_flag(:,end)==1,2),'color',[0 0 1],'linestyle','none','Marker','.','markersize',8)
    plot(output.x_CXn(output.CX_NEUTRALS&output.CX_flag(:,end)==2,1),output.x_CXn(output.CX_NEUTRALS&output.CX_flag(:,end)==2,2),'color',[0.2 0.0 0.7],'linestyle','none','Marker','.','markersize',12)
    plot(output.x_CXn(output.CX_NEUTRALS&output.CX_flag(:,end)==3,1),output.x_CXn(output.CX_NEUTRALS&output.CX_flag(:,end)==3,2),'color',[0.4 0.0 0.4],'linestyle','none','Marker','.','markersize',12)
    
end

MID_LOSS=find(output.ejected_wall_dt&output.x_ej(:,1)>R_OWL_UP&output.x_ej(:,2)>Z_OWL_DOWN&output.x_ej(:,2)<Z_OWL_UP);
HFS_MIS_LOSS=find(output.ejected_wall_dt&output.x_ej(:,1)<1&output.x_ej(:,2)>-0.3&output.x_ej(:,2)<0.3);
TOP_LOSS=find(output.ejected_wall_dt&output.x_ej(:,2)>Z_OWL_UP);
DOWN_LOSS=find(output.ejected_wall_dt&output.x_ej(:,2)<Z_OWL_DOWN);
CX_LOSS=find(output.ejected_CX_dt);
ION_LOSS=find(output.ions_ejected_wall_dt);
DID_CX_DT=find(output.did_CX_dt);

pfc_losses=struct();

pfc_losses.OWL_LOSS=MID_LOSS;
pfc_losses.HFS_MIS_LOSS=HFS_MIS_LOSS;
pfc_losses.TOP_LOSS=TOP_LOSS;
pfc_losses.OWN_LOSS=DOWN_LOSS;
pfc_losses.CX_LOSS_TOT=CX_LOSS;
pfc_losses.ION_LOSS=ION_LOSS;
summary.DID_CX_DT=DID_CX_DT;

% disp(['MID_LOSS  = ' num2str(length(MID_LOSS))]);
% disp(['HFS_MIS_LOSS  = ' num2str(length(HFS_MIS_LOSS))]);
% disp(['TOP_LOSS  = ' num2str(length(TOP_LOSS))]);
% disp(['DOWN_LOSS = ' num2str(length(DOWN_LOSS))]);
disp('during DT:')
disp(['CX_LOSS   = ' num2str(length(CX_LOSS))]);
disp(['DID CX   = ' num2str(length(DID_CX_DT))]);


POWER_OWL_LOSS=eV*par.MARKER_WEIGHT*sum(output.Ekin_ej(MID_LOSS))/DT_SIM_RECORD;
POWER_HFS_MIS_LOSS=eV*par.MARKER_WEIGHT*sum(output.Ekin_ej(HFS_MIS_LOSS))/DT_SIM_RECORD;
POWER_TOP_LOSS=eV*par.MARKER_WEIGHT*sum(output.Ekin_ej(TOP_LOSS))/DT_SIM_RECORD;
POWER_DOWN_LOSS=eV*par.MARKER_WEIGHT*sum(output.Ekin_ej(DOWN_LOSS))/DT_SIM_RECORD;
POWER_CX_LOSS=eV*par.MARKER_WEIGHT*sum(output.Ekin_ej(CX_LOSS))/DT_SIM_RECORD;
POWER_ION_LOSS=eV*par.MARKER_WEIGHT*sum(output.Ekin_ej(ION_LOSS))/DT_SIM_RECORD;
TOTAL_WALL_LOSS=eV*par.MARKER_WEIGHT*sum(output.Ekin_ej(output.ejected_wall_dt))/DT_SIM_RECORD;

POW_DID_CX_DT=eV*par.MARKER_WEIGHT*sum(mean(output.Ekin(DID_CX_DT,:),2))/DT_SIM_RECORD;

disp('---------------------------------------------------------------')


% disp(['POWER TOTAL_WALL_LOSS   = ' num2str(1e-3*round(TOTAL_WALL_LOSS))]);
% disp(['POWER_ION_LOSS   = ' num2str(1e-3*round(POWER_ION_LOSS))]);
% disp(['POWER_CX_neutrals_LOSS   = ' num2str(1e-3*round(POWER_CX_LOSS))]);
% disp(['POWER_OWL_LOSS  = ' num2str(1e-3*round(POWER_OWL_LOSS))]);
% disp(['POWER_HFS_MIS_LOSS  = ' num2str(1e-3*round(POWER_HFS_MIS_LOSS))]);
% disp(['POWER_TOP_LOSS  = ' num2str(1e-3*round(POWER_TOP_LOSS))]);
% disp(['POWER_DOWN_LOSS = ' num2str(1e-3*round(POWER_DOWN_LOSS))]);
% disp(['DID CX (POWER) = ' num2str(1e-3*round(POW_DID_CX_DT))]);

pfc_losses.TOTAL_WALL_POWER_LOSS=TOTAL_WALL_LOSS;
pfc_losses.POWER_ION_LOSS_TOT=POWER_ION_LOSS;
pfc_losses.POWER_OWL_LOSS=POWER_OWL_LOSS;
pfc_losses.POWER_HFS_MID_LOSS=POWER_HFS_MIS_LOSS;
pfc_losses.POWER_TOP_LOSS=POWER_TOP_LOSS;
pfc_losses.POWER_DOWN_LOSS=POWER_DOWN_LOSS;
pfc_losses.POWER_CX_LOSS=POWER_CX_LOSS;
pfc_losses.POW_DID_CX_DT=POW_DID_CX_DT;





IONS_RECYCLED_ONCE_DT=logical(squeeze(~ejected))&~output.ejected_sd&(output.CX_flag(:,end)-output.CX_flag(:,1))==1;
IONS_RECYCLED_ONCE_POP=find(IONS_RECYCLED_ONCE_DT);

EJ_CX_POP=find(output.ejected_CX_dt&output.CX_NEUTRALS);

output.EJ_CX_POP=EJ_CX_POP;
output.BORN_MARKERS_DT=BORN_MARKERS_DT;

%%
PLOT_QUIV=0;

if PLOT_QUIV

    for ii=1:20:length(EJ_CX_POP)
        
%         quiver(output.x_CXn(EJ_CX_POP(ii),1),output.x_CXn(EJ_CX_POP(ii),2),...
%             output.x_ej(EJ_CX_POP(ii),1)-output.x_CXn(EJ_CX_POP(ii),1), output.x_ej(EJ_CX_POP(ii),2)-output.x_CXn(EJ_CX_POP(ii),2),...
%             'color',[0.0 0.0 0.0],'linewidth',1)
    end
    
    plot(output.x_ej(output.ejected_wall_dt,1),output.x_ej(output.ejected_wall_dt,2),'g.','markersize',6)
    plot(output.x_ej(output.ejected_CX_dt,1),output.x_ej(find(output.ejected_CX_dt),2),'m.','markersize',4)
end
%%

Rpos_CX=output.x_CXn(DID_CX_DT,1);
Zpos_CX=output.x_CXn(DID_CX_DT,2);
Rpos_TOT=squeeze(output.x_gc(~ejected,1,end));
Zpos_TOT=squeeze(output.x_gc(~ejected,2,end));
Ekin_TOT=squeeze(output.Ekin(~ejected,1,end));
% if NB_TS_AVG>0
%     for TS=1:NB_TS_AVG
%         Rpos_TOT=[Rpos_TOT ; squeeze(output.x_gc(~ejected,1,TS+1))];
%         Zpos_TOT=[Zpos_TOT ; squeeze(output.x_gc(~ejected,2,TS+1))];
%         Ekin_TOT=[Ekin_TOT ; squeeze(output.Ekin(~ejected,TS+1))];
%     end
% end

Rbins = scale_X(1:10:end)+R0+0.02;
Zbins = scale_Z(1:30:end);
Rbins_scale = Rbins(1:end-1)+0.5*(Rbins(2)-Rbins(1));
Zbins_scale = Zbins(1:end-1)+0.5*(Zbins(2)-Zbins(1));

% Rbins=0.3:0.01:1.4;
% Zbins=-1.3:0.04:1.3;
% Rbins_scale=Rbins(1:end-1)+0.5*(Rbins(2)-Rbins(1));
% Zbins_scale=Zbins(1:end-1)+0.5*(Zbins(2)-Zbins(1));

NBI_power_map_2D=hist2d([Zpos_TOT Rpos_TOT ],Zbins,Rbins);
DV_map=2*pi*repmat(Rbins_scale, length(Zbins_scale),1);
DV_map=DV_map*(Rbins(2)-Rbins(1))*(Zbins(2)-Zbins(1));
for Rindex=1:length(Rbins_scale)
    for Zindex=1:length(Zbins_scale)
        POP_BIN_TOT=find((Rpos_TOT>=Rbins(Rindex))&(Rpos_TOT<Rbins(Rindex+1))&(Zpos_TOT>=Zbins(Zindex))&(Zpos_TOT<Zbins(Zindex+1)));
        NBI_power_map_2D(Zindex,Rindex)=sum(abs(Ekin_TOT(POP_BIN_TOT))); % total power gone in to CX
    end
end
NBI_power_map_2D=(eV)*(par.MARKER_WEIGHT)*NBI_power_map_2D./DV_map/DT_SIM_RECORD; % /(NB_TS_AVG+1);

summary.NBI_power_map_2D=NBI_power_map_2D;

if par.CALCULATE_CX
    CX_power_map_2D=hist2d([Zpos_CX Rpos_CX ],Zbins,Rbins);
%     NBI_power_map_2D=hist2d([Zpos_TOT Rpos_TOT ],Zbins,Rbins);
    CX_Ekin_power_map_2D=CX_power_map_2D;
%     DV_map=2*pi*repmat(Rbins_scale, length(Zbins_scale),1);
%     DV_map=DV_map*(Rbins(2)-Rbins(1))*(Zbins(2)-Zbins(1));
    
    for Rindex=1:length(Rbins_scale)
        for Zindex=1:length(Zbins_scale)
            POP_BIN=find((Rpos_CX>=Rbins(Rindex))&(Rpos_CX<Rbins(Rindex+1))&(Zpos_CX>=Zbins(Zindex))&(Zpos_CX<Zbins(Zindex+1)));
%             POP_BIN_TOT=find((Rpos_TOT>=Rbins(Rindex))&(Rpos_TOT<Rbins(Rindex+1))&(Zpos_TOT>=Zbins(Zindex))&(Zpos_TOT<Zbins(Zindex+1)));
            CX_Ekin_power_map_2D(Zindex,Rindex)=sum(mean(abs(output.Ekin(POP_BIN,:)),2));
%             CX_power_map_2D(Zindex,Rindex)=sum(mean(abs(output.Ekin(POP_BIN,:)),2)); % total power gone in to CX
%             NBI_power_map_2D(Zindex,Rindex)=sum(abs(Ekin_TOT(POP_BIN_TOT))); % total power gone in to CX
        end
    end
%     NBI_power_map_2D=(eV)*(par.MARKER_WEIGHT)*NBI_power_map_2D./DV_map/DT_SIM_RECORD; % /(NB_TS_AVG+1);
%     NBI_power_map_2D=max(NBI_power_map_2D,CX_Ekin_power_map_2D);
    CX_Ekin_power_map_2D=(eV)*(par.MARKER_WEIGHT)*CX_Ekin_power_map_2D./DV_map/DT_SIM_RECORD;
    CX_power_map_2D=CX_Ekin_power_map_2D./NBI_power_map_2D;  % local fraction of NBI power that went into CX
    CX_power_map_2D(NBI_power_map_2D<=1)=0;

    summary.Rvals_cx_pow=Rbins_scale;
    summary.Zvals_cx_pow=Zbins_scale;
    summary.DV_map_cx_pow=DV_map;
    summary.CX_power_map_2D=CX_power_map_2D;
    summary.CX_Ekin_power_map_2D=CX_Ekin_power_map_2D';

    if PLOT_CX_POWER_2D
        figure
        imagesc(Rbins_scale,Zbins_scale,CX_Ekin_power_map_2D);
        plot_TOK_background_par;
        set(gca,'fontsize',20)
        title([TITLESTRING 'neutralized to CX [W/m^{3}]'],'fontsize',18)
%         ylim([-0.62 0.6])
        xlim([min(Rbins_scale) max(Rbins_scale) ])
    end
    
    % end
    %
    %
    % if PLOT_CX_POWER_2D
    
    Rpos_neut=output.x_CXi(DID_CX_DT,1);
    Zpos_neut=output.x_CXi(DID_CX_DT,2);
    neutrals_power_map_2D=hist2d([Rpos_neut Zpos_neut ],Zbins,Rbins);
    
    for Rindex=1:length(Rbins_scale)
        for Zindex=1:length(Zbins_scale)
            POP_BIN=find((Rpos_neut>=Rbins(Rindex))&(Rpos_neut<Rbins(Rindex+1))&(Zpos_neut>=Zbins(Zindex))&(Zpos_neut<Zbins(Zindex+1)));
            neutrals_power_map_2D(Zindex,Rindex)=sum(mean(abs(output.Ekin(POP_BIN,:)),2));
        end
    end
    neutrals_power_map_2D=(eV)*(par.MARKER_WEIGHT)*neutrals_power_map_2D./DV_map/DT_SIM_RECORD;
    
    if PLOT_CX_POWER_2D
        figure
        imagesc(Rbins_scale,Zbins_scale,neutrals_power_map_2D);
        plot_TOK_background_par;
        set(gca,'fontsize',20)
        title([TITLESTRING 'neutrals re-absorbed [W/m^{3}]'],'fontsize',20)
        xlim([min(Rbins_scale) max(Rbins_scale) ])
    end
    
    summary.Rvals_cx_pow=Rbins_scale;
    summary.Zvals_cx_pow=Zbins_scale;
    summary.DV_map_cx_pow=DV_map;
    summary.neutrals_power_map_2D=neutrals_power_map_2D';
    
    Raxis=R0+a;
    Rsep=interp1(scale_psi_mp,scale_R_mp,psi_scale(end));

    figure
    set(gca,'fontsize',24)
    hold on
    grid on
    plot(Rbins_scale,sum(CX_Ekin_power_map_2D,1),'linewidth',3);
    plot(Rbins_scale,sum(neutrals_power_map_2D,1),'linestyle','-.','linewidth',3);
    xlabel('R [m]')
%     ylabel('[W/m^{3}]')
    ylabel('[W/m]')
    xlim([0.3 1.4])
%     title('<P_{loss CX}>')
    title('P_{loss CX}')
    
    legend('CX neutralization','re-ionizations')
    
    Rline=[Raxis Raxis];    
%     Z_line=[0 4e5];
    Z_line=[0 7e6];
    hl=line(Rline,Z_line);
    set(hl,'linewidth',3);
    set(hl,'color','k');
    set(hl,'linestyle','--');
    
    Rline=[Rsep Rsep];
    Z_line=[0 4e5];
    hl=line(Rline,Z_line);
    set(hl,'linewidth',3);
    set(hl,'color','k');
    set(hl,'linestyle','--')

    
end


%%

if PLOT_NEUTRALS_RELOC & par.CALCULATE_CX
    figure
    
    imagesc(Rbins_scale,Zbins_scale,hist2d([Zpos_CX Rpos_CX ],Zbins,Rbins));
    
    plot_TOK_background_par;
    hold on
    disp('ejected CX markers to wall : ')
    disp(length(find(EJ_CX_POP)))
    disp('Latest neutrals that got reionized : ')
    disp(length(find(IONS_RECYCLED_ONCE_DT)))
    
    scale_factor=1.0;
    POP_FRAC=42;
    
    plot(output.x_CXi(IONS_RECYCLED_ONCE_POP(1:POP_FRAC:end),1),output.x_CXi(IONS_RECYCLED_ONCE_POP(1:POP_FRAC:end),2),'r.','markersize',8)
    
    for ii=1:POP_FRAC:length(IONS_RECYCLED_ONCE_POP)
        
        quiver(output.x_CXn(IONS_RECYCLED_ONCE_POP(ii),1),output.x_CXn(IONS_RECYCLED_ONCE_POP(ii),2),...
            scale_factor*(output.x_CXi(IONS_RECYCLED_ONCE_POP(ii),1)-output.x_CXn(IONS_RECYCLED_ONCE_POP(ii),1)), scale_factor*(output.x_CXi(IONS_RECYCLED_ONCE_POP(ii),2)-output.x_CXn(IONS_RECYCLED_ONCE_POP(ii),2)),1,...
            'Marker','o',...
            'color',[0.0 1.0 0.1],'linewidth',1.3);
    end
    plot(output.x_CXn(IONS_RECYCLED_ONCE_POP(1:POP_FRAC:end),1),output.x_CXn(IONS_RECYCLED_ONCE_POP(1:POP_FRAC:end),2),'b.','markersize',8)
%     xlim([min(Rbins_scale) max(Rbins_scale)])
    
end




if PLOT_NEUTRALS_RELOC & par.CALCULATE_CX
    figure
    imagesc(Rbins_scale,Zbins_scale,hist2d([Zpos_CX Rpos_CX ],Zbins,Rbins));
    plot_TOK_background_par;
    hold on
    disp('ejected CX markers to wall : ')
    disp(length(find(EJ_CX_POP)))
    
    scale_factor=1.0;
    POP_FRAC=100;
    
    plot(output.x_ej(EJ_CX_POP(1:POP_FRAC:end),1),output.x_ej(EJ_CX_POP(1:POP_FRAC:end),2),'m.','markersize',11)
    
    for ii=1:POP_FRAC:length(EJ_CX_POP)
        
        quiver(output.x_CXn(EJ_CX_POP(ii),1),output.x_CXn(EJ_CX_POP(ii),2),...
            scale_factor*(output.x_ej(EJ_CX_POP(ii),1)-output.x_CXn(EJ_CX_POP(ii),1)), scale_factor*(output.x_ej(EJ_CX_POP(ii),2)-output.x_CXn(EJ_CX_POP(ii),2)),1,...
            'Marker','o',...
            'color',[1.0 0.0 0.1],'linewidth',1.3);
    end
    plot(output.x_CXn(EJ_CX_POP(1:POP_FRAC:end),1),output.x_CXn(EJ_CX_POP(1:POP_FRAC:end),2),'b.','markersize',8)
    xlim([min(Rbins_scale) max(Rbins_scale)])
    
end

% save(par.FILENAME,'-append','output');

%%
if par.CALCULATE_CX
    if PLOT_HIST_CX_SUMMARY
        figure
        % title('CX loss energies')
        hist(output.Ekin_ej(output.ejected_CX_dt)*1e-3,linspace(2,82,42));
        xlim([2 82])
        hold on
        set(gca,'fontsize',16)
        title([TITLESTRING ' ' num2str(round(summary.POWER_CX_LOSS*1e-3)) ' kW / 5MW NBI CX loss : Ekin dist '])
        plot([75 75],[0 2000],'k--','linewidth',2)
        plot(0.5.*[75 75],[0 2000],'k--','linewidth',2)
        plot(0.3333.*[75 75],[0 2000],'k--','linewidth',2)
        
        xlabel('Ekin [keV]')
        
    end
    
    CX_OWL_LOSS=find(output.ejected_CX_dt&output.x_ej(:,1)>R_OWL_UP&output.x_ej(:,2)>Z_OWL_DOWN&output.x_ej(:,2)<Z_OWL_UP);
    CX_HFS_MID_LOSS=find(output.ejected_CX_dt&output.x_ej(:,1)<1&output.x_ej(:,2)>-0.3&output.x_ej(:,2)<0.3);
    CX_TOP_LOSS=find(output.ejected_CX_dt&output.x_ej(:,2)>Z_OWL_UP);
    CX_DOWN_LOSS=find(output.ejected_CX_dt&output.x_ej(:,2)<Z_OWL_DOWN);
    
    % CX_MID_LOSS=find(output.ejected_CX_dt&output.x_ej(:,2)>-0.24&output.x_ej(:,2)<0.24);
    % CX_TOP_LOSS=find(output.ejected_CX_dt&output.x_ej(:,2)>0.4);
    % CX_DOWN_LOSS=find(output.ejected_CX_dt&output.x_ej(:,2)<-0.3);
    
    % disp(['CX_MID_LOSS  = ' num2str(length(CX_MID_LOSS))]);
    % disp(['CX_TOP_LOSS  = ' num2str(length(CX_TOP_LOSS))]);
    % disp(['CX_DOWN_LOSS = ' num2str(length(CX_DOWN_LOSS))]);
    
    CX_POWER_OWL_LOSS=eV*par.MARKER_WEIGHT*sum(output.Ekin_ej(CX_OWL_LOSS))/DT_SIM_RECORD;
    CX_POWER_TOP_LOSS=eV*par.MARKER_WEIGHT*sum(output.Ekin_ej(CX_TOP_LOSS))/DT_SIM_RECORD;
    CX_POWER_DOWN_LOSS=eV*par.MARKER_WEIGHT*sum(output.Ekin_ej(CX_DOWN_LOSS))/DT_SIM_RECORD;
    
    disp('---------------------------------------------------------------')
    
    disp(['CX_POWER_OWL_LOSS  = ' num2str(1e-3*round(CX_POWER_OWL_LOSS))]);
    disp(['CX_POWER_TOP_LOSS  = ' num2str(1e-3*round(CX_POWER_TOP_LOSS))]);
    disp(['CX_POWER_DOWN_LOSS = ' num2str(1e-3*round(CX_POWER_DOWN_LOSS))]);
    
    pfc_losses.CX_OWL_LOSS=CX_OWL_LOSS;
    pfc_losses.CX_HFS_MID_LOSS=CX_HFS_MID_LOSS;
    pfc_losses.CX_TOP_LOSS=CX_TOP_LOSS;
    pfc_losses.CX_DOWN_LOSS=CX_DOWN_LOSS;
    pfc_losses.CX_POWER_OWL_LOSS=CX_POWER_OWL_LOSS;
    pfc_losses.CX_POWER_TOP_LOSS=CX_POWER_TOP_LOSS;
    pfc_losses.CX_POWER_DOWN_LOSS=CX_POWER_DOWN_LOSS;
    
    
    Ekin_hist_values=(5:5:85)*1e3;
    
    AVG_TIME=round(10*DT_SIM_RECORD*1e3)/10;
    disp(['AVG_TIME = ' num2str(AVG_TIME)]);
    
    if PLOT_HIST_LOC
        figure
        subplot(3,1,1)
        hist(output.Ekin_ej(CX_TOP_LOSS),Ekin_hist_values)
        set(gca,'fontsize',16)
        title([TITLESTRING 'CX top loss Ekin [' num2str(AVG_TIME) 'ms]'])
        
        subplot(3,1,2)
        hist(output.Ekin_ej(CX_OWL_LOSS),Ekin_hist_values)
        set(gca,'fontsize',16)
        title([TITLESTRING 'CX OWL loss Ekin [' num2str(AVG_TIME) 'ms]'])
        
        subplot(3,1,3)
        hist(output.Ekin_ej(CX_DOWN_LOSS),Ekin_hist_values)
        set(gca,'fontsize',16)
        title([TITLESTRING 'CX lower loss Ekin [' num2str(AVG_TIME) 'ms]'])
        
    end
    
end
%%
Ekin_hist_values=(5:5:85)*1e3;

if PLOT_HIST_LOC & par.CALCULATE_CX
    
    figure
    subplot(3,1,1)
    hist(output.Ekin_ej(TOP_LOSS),Ekin_hist_values)
    set(gca,'fontsize',16)
    title([TITLESTRING 'top loss Ekin [' num2str(AVG_TIME) 'ms]'])
    
    subplot(3,1,2)
    hist(output.Ekin_ej(MID_LOSS),Ekin_hist_values)
    set(gca,'fontsize',16)
    title([TITLESTRING 'OWL mid-plane loss Ekin [' num2str(AVG_TIME) 'ms]'])
    
    subplot(3,1,3)
    hist(output.Ekin_ej(DOWN_LOSS),Ekin_hist_values)
    set(gca,'fontsize',16)
    title([TITLESTRING 'lower loss Ekin [' num2str(AVG_TIME) 'ms]'])
    
    
end


if PLOT_HIST_IONS_NEUTRALS & par.CALCULATE_CX
    MID_LOSS_IONS=find(~output.ejected_CX&output.ejected_wall_dt&output.x_ej(:,1)>R_OWL_UP&output.x_ej(:,2)>Z_OWL_DOWN&output.x_ej(:,2)<Z_OWL_UP);
    MID_LOSS_NEUTRALS=find(output.ejected_CX&output.ejected_wall_dt&output.x_ej(:,1)>R_OWL_UP&output.x_ej(:,2)>Z_OWL_DOWN&output.x_ej(:,2)<Z_OWL_UP);

    subplot(2,1,1)
    hist(output.Ekin_ej(MID_LOSS_IONS),Ekin_hist_values)
    set(gca,'fontsize',16)
    title([TITLESTRING 'Outer limiter ions loss'])
    xlabel('Ekin [eV]')
    xlim([2.5 82.5]*1e3)
    
    subplot(2,1,2)
    hist(output.Ekin_ej(MID_LOSS_NEUTRALS),Ekin_hist_values)
    set(gca,'fontsize',16)
    title([TITLESTRING 'Outer limiter neutrals loss'])
    xlabel('Ekin [eV]')   
    xlim([2.5 82.5]*1e3)
    
end


%%
% output.ejected_wall_dt=output.ejected_wall&output.time_stamp_loss>(par.NB_STAMPS_saved-par.NB_TS_RECORDED+1);

% Rmesh = repmat(scale_X+R0, length(scale_Z),1)';
% Zmesh = repmat(scale_Z, length(scale_X),1);
if isfield(output,'ndd_val')
    
    Rpos=squeeze(output.x(~ejected,1,end));
    Zpos=squeeze(output.x(~ejected,2,end));
    ndd_values=squeeze(output.ndd_val(:,end));
    output.ndd_val_thermal(isnan(output.ndd_val_thermal))=0;   % security for ndd th bug
     output.ndd_val_thermal_evol=sum(output.ndd_val_thermal,1);
    ndd_thermal_values=squeeze(output.ndd_val_thermal(~ejected,end));
    ndd_thermal_values=ndd_thermal_values.*(output.ndd_val_thermal_evol(end)./sum(ndd_thermal_values));

    
    ndd_BB_values=squeeze(output.ndd_val_BB(~ejected,end));
    ndd_BP_values=squeeze(output.ndd_val_BP(~ejected,end));

    ndd_tot_values=ndd_thermal_values+ndd_BB_values+ndd_BP_values;
    
    Rbins = scale_X(1:27:end)+R0;Rbins=[Rbins Rbins(end)+(Rbins(2)-Rbins(1))];
    Zbins = scale_Z(1:27:end);
    Rvals = Rbins(1:end-1)+0.5*(Rbins(2)-Rbins(1));
    Zvals = Zbins(1:end-1)+0.5*(Zbins(2)-Zbins(1));
    
    DV_map=2*pi*repmat(Rvals, length(Zvals),1)';
    DV_map=DV_map*(Rbins(2)-Rbins(1))*(Zbins(2)-Zbins(1));
    
    ndd_map=zeros(length(Rvals),length(Zvals));
    ndd_BP_map=zeros(length(Rvals),length(Zvals));
    ndd_BB_map=zeros(length(Rvals),length(Zvals));
    ndd_tot_map=zeros(length(Rvals),length(Zvals));
    
    for Rindex=1:length(Rvals)
        for Zindex=1:length(Zvals)
            POP_BIN=find((Rpos>=Rbins(Rindex))&(Rpos<Rbins(Rindex+1))&(Zpos>=Zbins(Zindex))&(Zpos<Zbins(Zindex+1)));
            ndd_map(Rindex,Zindex)=sum(ndd_values(POP_BIN));
            ndd_BP_map(Rindex,Zindex)=sum(ndd_BP_values(POP_BIN));
            ndd_BB_map(Rindex,Zindex)=sum(ndd_BB_values(POP_BIN));
            ndd_tot_map(Rindex,Zindex)=sum(ndd_tot_values(POP_BIN));
        end
    end
    
    disp('Building nDD maps......')
    if NB_TS_AVG>0
        for TS=1:NB_TS_AVG

            n_ej=logical(output.time_stamp_loss>=par.NB_TIME_STAMPS-TS)|(isnan(output.time_stamp_loss));
            
            BORN_MARKERS_DT=logical(birth_matrix(end-TS,:)');
            POP_TS=logical(n_ej.*BORN_MARKERS_DT);
            
            Rpos=squeeze(output.x(POP_TS,1,end-TS));
            Zpos=squeeze(output.x(POP_TS,2,end-TS));
            ndd_values=squeeze(output.ndd_val(POP_TS,end-TS));
            ndd_BB_values=squeeze(output.ndd_val_BB(POP_TS,end-TS));
            ndd_BP_values=squeeze(output.ndd_val_BP(POP_TS,end-TS));
            ndd_thermal_values=squeeze(output.ndd_val_thermal(POP_TS,end-TS));
            ndd_thermal_values=ndd_thermal_values.*(output.ndd_val_thermal_evol(end-TS)./sum(ndd_thermal_values));
            
            ndd_tot_values=ndd_thermal_values+ndd_BP_values+ndd_BB_values;
            

            for Rindex=1:length(Rvals)
                for Zindex=1:length(Zvals)
                    POP_BIN=find((Rpos>=Rbins(Rindex))&(Rpos<Rbins(Rindex+1))&(Zpos>=Zbins(Zindex))&(Zpos<Zbins(Zindex+1)));
                    ndd_map(Rindex,Zindex)=ndd_map(Rindex,Zindex)+sum(ndd_values(POP_BIN));
                    ndd_BP_map(Rindex,Zindex)=ndd_BP_map(Rindex,Zindex)+sum(ndd_BP_values(POP_BIN));
                    ndd_BB_map(Rindex,Zindex)=ndd_BB_map(Rindex,Zindex)+sum(ndd_BB_values(POP_BIN));
                    ndd_tot_map(Rindex,Zindex)=ndd_tot_map(Rindex,Zindex)+sum(ndd_tot_values(POP_BIN));
                end
            end
        end
        DT_SIM_RECORD=(par.time_scale(end)-par.time_scale(1));
        
    else
        DT_SIM_RECORD=par.dt*par.NR_FUND_IN_LOOP*par.TIME_STAMP_PRECISION;
    end

    
    ndd_map=ndd_map/(NB_TS_AVG+1);
    ndd_BP_map=ndd_BP_map/(NB_TS_AVG+1);
    ndd_BB_map=ndd_BB_map/(NB_TS_AVG+1);
    ndd_tot_map=ndd_tot_map/(NB_TS_AVG+1);
    % ndd_map=griddata(Rpos,Zpos,ndd_values,Zmesh,Rmesh);
    
    summary.Rvals_ndd=Rvals;
    summary.Zvals_ndd=Zvals;
    summary.DV_map_ndd=DV_map;
    summary.ndd_RZ=ndd_tot_map./DV_map;
    
    if PLOT_NDD
        if PLOT_NDD_DETAILS
            figure
            hold on
            imagesc(Rvals,Zvals,(ndd_map./DV_map)')
            plot_TOK_background_par
            ylim([min(Zvals) max(Zvals)])
            xlim([min(Rvals) max(Rvals)])
            
            set(gca,'fontsize',20)
            title([TITLESTRING 'Beam[' num2str(round(par.INPUT_POWER*1e-6)) ']MW=>target neut/s/m^{3} (old formula) '])
            colorbar
            %           title(['MAST-U Beam[5MW]=>target neut/s/m^{3} '])
            
            figure
            hold on
            imagesc(Rvals,Zvals,(ndd_BP_map./DV_map)')
            plot_TOK_background_par
            ylim([min(Zvals) max(Zvals)])
            xlim([min(Rvals) max(Rvals)])
            
            set(gca,'fontsize',20)
            title([TITLESTRING 'Beam[' num2str(round(par.INPUT_POWER*1e-6)) 'MW]=>target neut/s/m^{3} '])
            colorbar
            
            figure
            hold on
            imagesc(Rvals,Zvals,(ndd_BB_map./DV_map)')
            plot_TOK_background_par
            ylim([min(Zvals) max(Zvals)])
            xlim([min(Rvals) max(Rvals)])
            
            set(gca,'fontsize',20)
            title([TITLESTRING 'Beam[' num2str(round(par.INPUT_POWER*1e-6)) 'MW]=>Beam neut/s/m^{3} '])
            colorbar
        end
        
        fh=figure(99)
        hold on
        pcolor(summary.Rvals_ndd,summary.Zvals_ndd,summary.ndd_RZ')
        shading interp;
        plot_TOK_background_par
        ylim([min(summary.Zvals_ndd)+0.5*(Zvals(2)-Zvals(1)) max(summary.Zvals_ndd)-0.5*(Zvals(2)-Zvals(1))])
        xlim([min(summary.Rvals_ndd)+0.5*(Rvals(2)-Rvals(1)) max(summary.Rvals_ndd)-0.5*(Rvals(2)-Rvals(1))])
        xlabel('R [m]');
        ylabel('Z [m]')
        set(gca,'fontsize',22)
        title([TITLESTRING 'DD neut/s/m^{3} '])
        colorbar
%           title(['MAST-U Beam[5MW]=>target neut/s/m^{3} '])
        fh.Position = [1.8000   41.8000  600.0000  836.8000];
      
    end
    
end
%%

if isfield(output,'ejected_CX_dt')&par.CALCULATE_CX
    
    Rpos=squeeze(output.x_ej(output.ejected_CX_dt,1));
    Zpos=squeeze(output.x_ej(output.ejected_CX_dt,2));
    pow_values=squeeze(output.Ekin_ej(output.ejected_CX_dt))/DT_SIM_RECORD;
    
    Rbins = scale_X(1:20:end)+R0;
    Rbins = [Rbins(1)-(Rbins(2)-Rbins(1)) Rbins Rbins(end)+(Rbins(2)-Rbins(1))];
    Zbins = scale_Z(1:30:end);
    Rvals = Rbins(1:end-1)+0.5*(Rbins(2)-Rbins(1));
    Zvals = Zbins(1:end-1)+0.5*(Zbins(2)-Zbins(1));
    
    Rpos_map=2*pi*repmat(Rvals, length(Zvals),1)';
    DV_map=Rpos_map*(Rbins(2)-Rbins(1))*(Zbins(2)-Zbins(1));
    
    CXloss_map=zeros(length(Rvals),length(Zvals));
    
    for Rindex=1:length(Rvals)
        for Zindex=1:length(Zvals)
            POP_BIN=find((Rpos>=Rbins(Rindex))&(Rpos<Rbins(Rindex+1))&(Zpos>=Zbins(Zindex))&(Zpos<Zbins(Zindex+1)));
            CXloss_map(Rindex,Zindex)=sum(pow_values(POP_BIN));
        end
    end
    CXloss_map=CXloss_map*eV*par.MARKER_WEIGHT;
    % ndd_map=griddata(Rpos,Zpos,ndd_values,Zmesh,Rmesh);
    
    pfc_losses.Rvals_cxloss=Rvals;
    pfc_losses.Zvals_cxloss=Zvals;
    pfc_losses.DV_map_cxloss=DV_map;
    pfc_losses.DR=(Rbins(2)-Rbins(1));
    pfc_losses.cxloss_RZ=CXloss_map.*(Rbins(2)-Rbins(1))./DV_map;
    
    if PLOT_CXLOSS
        figure
        hold on
        imagesc(Rvals,Zvals,(CXloss_map.*(Rbins(2)-Rbins(1))./DV_map)')
        plot_TOK_background_par
%         ylim([-0.6 0.6])
        xlim([min(scale_X)+R0+0.02 max(scale_X)+R0])
        
        set(gca,'fontsize',20)
        title([TITLESTRING 'CX neutrals losses W/m^{2}'],'fontsize',18)
%         xlim([0.58 1.22])
        
    end
    
end

%%

if isfield(output,'ejected_wall_dt')
    EJ_ION_POP_DT=output.ejected_wall_dt&~output.ejected_CX_dt;
    
    Rpos=squeeze(output.x_ej(EJ_ION_POP_DT,1));
    Zpos=squeeze(output.x_ej(EJ_ION_POP_DT,2));
    pow_values=squeeze(output.Ekin_ej(EJ_ION_POP_DT))/DT_SIM_RECORD;

    %same bins as cx neutrals losses map
    Rbins = scale_X(1:20:end)+R0;
    Rbins = [Rbins(1)-(Rbins(2)-Rbins(1)) Rbins Rbins(end)+(Rbins(2)-Rbins(1))];
    Zbins = scale_Z(1:30:end);
    Rvals = Rbins(1:end-1)+0.5*(Rbins(2)-Rbins(1));
    Zvals = Zbins(1:end-1)+0.5*(Zbins(2)-Zbins(1));

%     Rbins = scale_X(1:25:end)+R0+0.02;
%     Zbins = scale_Z(1:30:end);
%     Rvals = Rbins(1:end-1)+0.5*(Rbins(2)-Rbins(1));
%     Zvals = Zbins(1:end-1)+0.5*(Zbins(2)-Zbins(1));
    
    Rpos_map=2*pi*repmat(Rvals, length(Zvals),1)';
    DV_map=Rpos_map*(Rbins(2)-Rbins(1))*(Zbins(2)-Zbins(1));
    
    ionloss_map=zeros(length(Rvals),length(Zvals));
    
    for Rindex=1:length(Rvals)
        for Zindex=1:length(Zvals)
            POP_BIN=find((Rpos>=Rbins(Rindex))&(Rpos<Rbins(Rindex+1))&(Zpos>=Zbins(Zindex))&(Zpos<Zbins(Zindex+1)));
            ionloss_map(Rindex,Zindex)=sum(pow_values(POP_BIN));
        end
    end
    ionloss_map=ionloss_map*eV*par.MARKER_WEIGHT;
    % ndd_map=griddata(Rpos,Zpos,ndd_values,Zmesh,Rmesh);
    
%     summary.Rvals_cxloss=Rvals;
%     summary.Zvals_cxloss=Zvals;
%     summary.DV_map_cxloss=DV_map;
    pfc_losses.ionloss_RZ=ionloss_map.*(Rbins(2)-Rbins(1))./DV_map;
    
    if PLOT_IONLOSS
        figure
        hold on
        imagesc(Rvals,Zvals,(ionloss_map.*(Rbins(2)-Rbins(1))./DV_map)')
        plot_TOK_background_par
        ylim([-0.6 0.6])
        xlim([min(scale_X)+R0+0.02 max(scale_X)+R0])
        
        set(gca,'fontsize',20)
        title([TITLESTRING 'ions losses W/m^{2}'],'fontsize',18)
        xlim([0.58 1.22])
        
    end
    
end


try
    save(FNAME,'-append','summary')
    save(FNAME,'-append','pfc_losses')
%     save(par.FILENAME,'-append','summary')
%     save(par.FILENAME,'-append','pfc_losses')
catch
    warning('could not save summary to file!!!')
end


%%
if SAVE_NEW_PREC_INPUT
    new_input=struct();
    old_input=input;
    
    new_input.m=old_input.m;
    new_input.Z=old_input.Z;
    new_input.x=zeros(size(squeeze(output.x(~ejected,:,end))));
    new_input.v=zeros(size(squeeze(output.v(~ejected,:,end))));
    
    new_input.x=squeeze(output.x(~ejected,:,end));
    new_input.v=squeeze(output.v(~ejected,:,end));
    
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

%%
extract_radial_profile_from_markers;

%%
% adjusted to be able to get particle BEFORE it ihits the wall

% old_input=input;
% 
% EJ_WALL_DT=find(output.ejected_wall_dt);
% 
% new_input=struct();
% input=struct();
% 
% input.m=old_input.m;
% input.Z=old_input.Z;
% input.x=zeros(length(EJ_WALL_DT),3);
% input.v=zeros(length(EJ_WALL_DT),3);
% 
% 
% disp('Creating sample file with inputs for trajectories of loss particles : ejected to wall during recorded TS in this file')
% 
% new_input.x=squeeze(output.x(output.ejected_wall_dt,:,:));
% new_input.v=squeeze(output.v(output.ejected_wall_dt,:,:));
% 
% for part_ej=1:length(find(output.ejected_wall_dt))
%     input.x(part_ej,:)=squeeze(new_input.x(part_ej,:,round(output.time_stamp_loss(EJ_WALL_DT(part_ej))-(par.NB_STAMPS_saved-par.NB_TS_RECORDED)-1)));
%     input.v(part_ej,:)=squeeze(new_input.v(part_ej,:,round(output.time_stamp_loss(EJ_WALL_DT(part_ej))-(par.NB_STAMPS_saved-par.NB_TS_RECORDED)-1)));
% end
% 
% input.N_total=size(input.x,1)
% 
% save('..\input\TF3D_losses_input.mat','input')
% 
% 
% input=old_input;
