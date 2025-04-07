
% clear all

try
load('../../../data_common/wallRZ_CU.mat')
catch
load('../data_common/wallRZ_CU.mat')
end

SAVE_NEW_PREC_INPUT=0;

PLOT_SD=1;
PLOT_NEUTRALS_RELOC=0;
PLOT_CX_POWER_2D=1;
PLOT_DOTS=0;
PLOT_HIST=0;
PLOT_CXLOSS=1;
PLOT_IONLOSS=1;
PLOT_NDD=0;

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

if isfield(output,'v')
    v2=squeeze(sum(output.v.^2,2));
end 
if ~isfield(output,'Ekin')
   
    ejected=logical(ejected);
    
    % moved to combine_output_eq
    % output.time_stamp_loss=ceil(output.time_step_loss.*par.NB_TIME_STAMPS/par.NB_TIME_STEPS);
    %
    % save(par.FILENAME,'-append','output');
    
    output.Ekin=0.5*(input.m/eV).*v2(:,:);
end

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

NB_TS_AVG=par.NB_TS_RECORDED-1;

%%
disp('lost markers during overall simulation :')
disp(length(find(ejected)))

output.Ekin_end=0.5*(input.m/eV).*v2(:,end);
output.pitch_end=output.vpll(:,end)./sqrt(v2(:,end));

Ekin_bins=(0:4:100)*1e3;
pitch_bins=-1.00:0.04:1.00;
Ekin_values=Ekin_bins(1:end-1)+0.5*(Ekin_bins(2)-Ekin_bins(1));
pitch_values=pitch_bins(1:end-1)+0.5*(pitch_bins(2)-pitch_bins(1));

n_ej=logical(~ejected);

BORN_MARKERS=logical(birth_matrix(end,:)');
POP_TS=logical(n_ej.*BORN_MARKERS)&~output.CX_NEUTRALS; %&squeeze(output.x_gc(:,1,1))>1.1;
dist_vspace=hist2d([output.Ekin_end(POP_TS) output.pitch_end(POP_TS)],Ekin_bins,pitch_bins);

%%
TS=1;
if NB_TS_AVG>0
    for TS=1:NB_TS_AVG
        output.Ekin_end=squeeze(output.Ekin(:,end-TS));
        output.pitch_end=output.vpll(:,end-TS)./sqrt(v2(:,end-TS));
        
        n_ej=logical(logical(output.time_stamp_loss>par.NB_TIME_STAMPS-TS)+(isnan(output.time_stamp_loss))&~output.CX_NEUTRALS);
        
        BORN_MARKERS=logical(birth_matrix(end-TS,:)');
        POP_TS=logical(n_ej.*BORN_MARKERS); %&squeeze(output.x_gc(:,1,end-TS))>1.1;
        dist_vspace_TS=hist2d([output.Ekin_end(POP_TS) output.pitch_end(POP_TS)],Ekin_bins,pitch_bins);
        
        dist_vspace=dist_vspace+dist_vspace_TS;
    end
    DT_SIM_RECORD=(par.time_scale(end)-par.time_scale(1));

else
    DT_SIM_RECORD=par.dt*par.NR_FUND_IN_LOOP*par.TIME_STAMP_PRECISION;
end

disp(['DT_SIM_RECORD = ' num2str(DT_SIM_RECORD)]);

%%
Ekin_array=repmat(Ekin_values,1,length(pitch_values));
pitch_array=repmat(pitch_values,length(Ekin_values),1);
% scatter(Ekin_array(:),pitch_array(:),200,dist_vspace(:),'filled')
% hf=figure;imagesc(Ekin_values,pitch_values,dist_vspace');
% 
dist_norm=sum(sum(dist_vspace))*(Ekin_bins(2)-Ekin_bins(1))*(pitch_bins(2)-pitch_bins(1));
dist_fabien=dist_vspace/dist_norm;


% dist_fabien(end,end)=1.9e-5;

hf=figure;
set(gca,'fontsize',16)
hold on
grid on
surf(Ekin_values,pitch_values,dist_fabien'); shading interp; view(2)
xlabel('Ekin [eV]')
ylabel('pitch')


summary=struct();
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
disp('statistics for the total number of time stamps recoreded in this file');

output.B_ej=interp2(scale_X+R0,scale_Z,Btot_XZ_map',output.x_ej(:,1),output.x_ej(:,2));
% output.Ekin_ej=0.5*(input.m/eV).*output.vpll_ej(:,end).^2+output.mm_ej.*output.B_ej;

try
    output.ejected_sd_dt=logical(output.ejected_sd)&output.time_stamp_loss>(par.NB_STAMPS_saved-par.NB_TS_RECORDED);
catch
    output.ejected_sd_dt=logical(ejected*0);
end

try
%     output.ejected_CX_dt=logical(squeeze(output.ejected_CX))&~output.ejected_sd&output.time_stamp_loss>par.NB_STAMPS_saved-par.NB_TS_RECORDED;
    output.did_CX_dt=logical(squeeze(output.CX_flag(:,end)>0))&~output.ejected_sd&output.time_stamp_loss>(par.NB_STAMPS_saved-par.NB_TS_RECORDED);
    ISNAN=output.time_stamp_loss<=(par.NB_STAMPS_saved-par.NB_TS_RECORDED);
    output.did_CX_dt=~ISNAN&logical(squeeze(output.CX_flag(:,end)-output.CX_flag(:,1)>0));
%     output.ejected_CX_dt=logical(squeeze(output.ejected_CX))&output.did_CX_dt;
    output.ejected_CX_dt=~ISNAN&output.ejected_CX;
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


%%
R_OWL_UP=wall_CU.PFC1_OWL(1,1);
Z_OWL_UP=wall_CU.PFC1_OWL(1,2);
R_OWL_DOWN=wall_CU.PFC1_OWL(end,1);
Z_OWL_DOWN=wall_CU.PFC1_OWL(end,2);

disp(['OWL limiter top  (R,Z) coordinates : ' num2str(R_OWL_UP) ' , ' num2str(Z_OWL_UP)])
disp(['OWL limiter down (R,Z) coordinates : ' num2str(R_OWL_DOWN) ' , ' num2str(Z_OWL_DOWN)])

if PLOT_DOTS
    
    
    figure
    hold on
    disp('ejected markers to wall : [not from slowing down or CX]')
    disp(length(find(ejected&~output.ejected_sd&~output.ejected_CX_dt)))
    try
        plot(output.x_CXn(output.CX_NEUTRALS&output.CX_flag(:,end)>=1,1),output.x_CXn(output.CX_NEUTRALS&output.CX_flag(:,end)>=1,2),'color',[0 0 1],'linestyle','none','Marker','.','markersize',8)
        plot(output.x_ej(output.ejected_wall_dt&~output.ejected_CX_dt,1),output.x_ej(output.ejected_wall_dt&~output.ejected_CX_dt,2),'g.','markersize',6)
        plot(output.x_ej(output.ejected_CX_dt,1),output.x_ej(find(output.ejected_CX_dt),2),'m.','markersize',4)
        legend('neutralization','ion wall loss','CX wall loss')
    end
    plot_CU_TOK_background;
    
    % plot(output.x_CXn(output.CX_NEUTRALS&output.CX_flag(:,end)==1,1),output.x_CXn(output.CX_NEUTRALS&output.CX_flag(:,end)==1,2),'color',[0 0 1],'linestyle','none','Marker','.','markersize',10)
    plot(output.x_CXn(output.CX_NEUTRALS&output.CX_flag(:,end)==1,1),output.x_CXn(output.CX_NEUTRALS&output.CX_flag(:,end)==1,2),'color',[0 0 1],'linestyle','none','Marker','.','markersize',8)
    plot(output.x_CXn(output.CX_NEUTRALS&output.CX_flag(:,end)==2,1),output.x_CXn(output.CX_NEUTRALS&output.CX_flag(:,end)==2,2),'color',[0.2 0.0 0.7],'linestyle','none','Marker','.','markersize',6)
    plot(output.x_CXn(output.CX_NEUTRALS&output.CX_flag(:,end)==3,1),output.x_CXn(output.CX_NEUTRALS&output.CX_flag(:,end)==3,2),'color',[0.4 0.0 0.4],'linestyle','none','Marker','.','markersize',4)
    
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


disp(['POWER TOTAL_WALL_LOSS   = ' num2str(1e-3*round(TOTAL_WALL_LOSS))]);
disp(['POWER_ION_LOSS   = ' num2str(1e-3*round(POWER_ION_LOSS))]);
disp(['POWER_CX_LOSS   = ' num2str(1e-3*round(POWER_CX_LOSS))]);
disp(['POWER_OWL_LOSS  = ' num2str(1e-3*round(POWER_OWL_LOSS))]);
disp(['POWER_HFS_MIS_LOSS  = ' num2str(1e-3*round(POWER_HFS_MIS_LOSS))]);
disp(['POWER_TOP_LOSS  = ' num2str(1e-3*round(POWER_TOP_LOSS))]);
disp(['POWER_DOWN_LOSS = ' num2str(1e-3*round(POWER_DOWN_LOSS))]);
disp(['DID CX (POWER) = ' num2str(1e-3*round(POW_DID_CX_DT))]);

pfc_losses.TOTAL_WALL_POWER_LOSS=TOTAL_WALL_LOSS;
pfc_losses.POWER_ION_LOSS_TOT=POWER_ION_LOSS;
pfc_losses.POWER_OWL_LOSS=POWER_OWL_LOSS;
pfc_losses.POWER_HFS_MID_LOSS=POWER_HFS_MIS_LOSS;
pfc_losses.POWER_TOP_LOSS=POWER_TOP_LOSS;
pfc_losses.POWER_DOWN_LOSS=POWER_DOWN_LOSS;
pfc_losses.POWER_CX_LOSS=POWER_CX_LOSS;
pfc_losses.POW_DID_CX_DT=POW_DID_CX_DT;

summary.TOTAL_WALL_POWER_LOSS=TOTAL_WALL_LOSS;
summary.POWER_ION_WALL_LOSS=POWER_ION_LOSS;
summary.POWER_CX_LOSS=POWER_CX_LOSS;
summary.POWER_DID_CX=POW_DID_CX_DT;

% output.POWER_HFS_MID_LOSS=POWER_HFS_MID_LOSS;
% output.POWER_TOP_LOSS=POWER_TOP_LOSS;
% output.POWER_DOWN_LOSS=POWER_DOWN_LOSS;
output.POWER_CX_LOSS=POWER_CX_LOSS;

if par.NB_TS_RECORDED>1
    BORN_MARKERS=logical(birth_matrix(end-par.NB_TS_RECORDED+1,:)');
else
    BORN_MARKERS=logical(birth_matrix(end,:)');
end
EJ_CX_POP=find(output.ejected_CX_dt&output.CX_NEUTRALS);

IONS_RECYCLED_ONCE_DT=logical(squeeze(~ejected))&~output.ejected_sd&(output.CX_flag(:,end)-output.CX_flag(:,1))==1;
IONS_RECYCLED_ONCE_POP=find(IONS_RECYCLED_ONCE_DT);

output.EJ_CX_POP=EJ_CX_POP;
output.BORN_MARKERS_DT=BORN_MARKERS;

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

Rbins=0.54:0.01:1.24;
Zbins=-0.65:0.04:0.65;
Rbins_scale=Rbins(1:end-1)+0.5*(Rbins(2)-Rbins(1));
Zbins_scale=Zbins(1:end-1)+0.5*(Zbins(2)-Zbins(1));

if PLOT_CX_POWER_2D
    figure
    CX_power_map_2D=hist2d([Zpos_CX Rpos_CX ],Zbins,Rbins);
    DV_map=2*pi*repmat(Rbins_scale, length(Zbins_scale),1);
    DV_map=DV_map*(Rbins(2)-Rbins(1))*(Zbins(2)-Zbins(1));
    
   for Rindex=1:length(Rbins_scale)
        for Zindex=1:length(Zbins_scale)
            POP_BIN=find((Rpos_CX>=Rbins(Rindex))&(Rpos_CX<Rbins(Rindex+1))&(Zpos_CX>=Zbins(Zindex))&(Zpos_CX<Zbins(Zindex+1)));
            CX_power_map_2D(Zindex,Rindex)=sum(mean(abs(output.Ekin(POP_BIN,:)),2));
        end
    end
    CX_power_map_2D=(eV)*(par.MARKER_WEIGHT)*CX_power_map_2D./DV_map/DT_SIM_RECORD;

    imagesc(Rbins_scale,Zbins_scale,CX_power_map_2D);
    plot_CU_TOK_background;
    set(gca,'fontsize',18)
    title('NBI Charge-Exchange power [W/m^{3}]','fontsize',14)
    ylim([-0.62 0.6])
    xlim([0.55 1.22])
end

%%

if PLOT_NEUTRALS_RELOC
    figure
    
    imagesc(Rbins_scale,Zbins_scale,hist2d([Zpos_CX Rpos_CX ],Zbins,Rbins));
    
    plot_CU_TOK_background;
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
    
end




if PLOT_NEUTRALS_RELOC
    figure
    imagesc(Rbins_scale,Zbins_scale,hist2d([Zpos_CX Rpos_CX ],Zbins,Rbins));
    plot_CU_TOK_background;
    hold on
    disp('ejected CX markers to wall : ')
    disp(length(find(EJ_CX_POP)))
    
    scale_factor=1.0;
    POP_FRAC=42;
    
    plot(output.x_ej(EJ_CX_POP(1:POP_FRAC:end),1),output.x_ej(EJ_CX_POP(1:POP_FRAC:end),2),'g.','markersize',8)
    
    for ii=1:POP_FRAC:length(EJ_CX_POP)
        
        quiver(output.x_CXn(EJ_CX_POP(ii),1),output.x_CXn(EJ_CX_POP(ii),2),...
            scale_factor*(output.x_ej(EJ_CX_POP(ii),1)-output.x_CXn(EJ_CX_POP(ii),1)), scale_factor*(output.x_ej(EJ_CX_POP(ii),2)-output.x_CXn(EJ_CX_POP(ii),2)),1,...
            'Marker','o',...
            'color',[1.0 0.0 0.1],'linewidth',1.3);
    end
    plot(output.x_CXn(EJ_CX_POP(1:POP_FRAC:end),1),output.x_CXn(EJ_CX_POP(1:POP_FRAC:end),2),'b.','markersize',8)
    
end

% save(par.FILENAME,'-append','output');

%%
if par.CALCULATE_CX
    if PLOT_HIST_CX_SUMMARY
        figure
        % title('CX loss energies')
        hist(output.Ekin_ej(output.ejected_CX_dt),40)
        set(gca,'fontsize',16)
        title([par.ID ' : ' num2str(round(POWER_CX_LOSS*1e-3)) ' kW / MW CX loss : Ekin dist '])
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
    
    if PLOT_HIST
        figure
        subplot(3,1,1)
        hist(output.Ekin_ej(CX_TOP_LOSS),Ekin_hist_values)
        set(gca,'fontsize',16)
        title(['CX top loss Ekin [' num2str(AVG_TIME) 'ms]'])
        
        subplot(3,1,2)
        hist(output.Ekin_ej(CX_OWL_LOSS),Ekin_hist_values)
        set(gca,'fontsize',16)
        title(['CX OWL loss Ekin [' num2str(AVG_TIME) 'ms]'])
        
        subplot(3,1,3)
        hist(output.Ekin_ej(CX_DOWN_LOSS),Ekin_hist_values)
        set(gca,'fontsize',16)
        title(['CX lower loss Ekin [' num2str(AVG_TIME) 'ms]'])
        
    end
    
end
%%
Ekin_hist_values=(5:5:85)*1e3;

if PLOT_HIST
    
    figure
    subplot(3,1,1)
    hist(output.Ekin_ej(TOP_LOSS),Ekin_hist_values)
    set(gca,'fontsize',16)
    title(['top loss Ekin [' num2str(AVG_TIME) 'ms]'])
    
    subplot(3,1,2)
    hist(output.Ekin_ej(MID_LOSS),Ekin_hist_values)
    set(gca,'fontsize',16)
    title(['OWL mid-plane loss Ekin [' num2str(AVG_TIME) 'ms]'])
    
    subplot(3,1,3)
    hist(output.Ekin_ej(DOWN_LOSS),Ekin_hist_values)
    set(gca,'fontsize',16)
    title(['lower loss Ekin [' num2str(AVG_TIME) 'ms]'])
    
    
end

%%
% output.ejected_wall_dt=output.ejected_wall&output.time_stamp_loss>(par.NB_STAMPS_saved-par.NB_TS_RECORDED+1);

% Rmesh = repmat(scale_X+R0, length(scale_Z),1)';
% Zmesh = repmat(scale_Z, length(scale_X),1);
if isfield(output,'ndd_val')
    
    Rpos=squeeze(output.x(~output.ejected_wall_dt,1,end));
    Zpos=squeeze(output.x(~output.ejected_wall_dt,2,end));
    ndd_values=squeeze(output.ndd_val(~output.ejected_wall_dt,end));
    
    Rbins = scale_X(1:15:end)+R0;
    Zbins = scale_Z(1:15:end);
    Rvals = Rbins(1:end-1)+0.5*(Rbins(2)-Rbins(1));
    Zvals = Zbins(1:end-1)+0.5*(Zbins(2)-Zbins(1));
    
    DV_map=2*pi*repmat(Rvals, length(Zvals),1)';
    DV_map=DV_map*(Rbins(2)-Rbins(1))*(Zbins(2)-Zbins(1));
    
    ndd_map=zeros(length(Rvals),length(Zvals));
    
    for Rindex=1:length(Rvals)
        for Zindex=1:length(Zvals)
            POP_BIN=find((Rpos>=Rbins(Rindex))&(Rpos<Rbins(Rindex+1))&(Zpos>=Zbins(Zindex))&(Zpos<Zbins(Zindex+1)));
            ndd_map(Rindex,Zindex)=sum(ndd_values(POP_BIN));
        end
    end
    % ndd_map=griddata(Rpos,Zpos,ndd_values,Zmesh,Rmesh);
    
    summary.Rvals_ndd=Rvals;
    summary.Zvals_ndd=Zvals;
    summary.DV_map_ndd=DV_map;
    summary.ndd_RZ=ndd_map./DV_map;
    
    if PLOT_NDD
        figure
        hold on
        imagesc(Rvals,Zvals,(ndd_map./DV_map)')
        plot_CU_TOK_background
        ylim([-0.45 0.45])
        xlim([0.6 1.18])
        
        set(gca,'fontsize',16)
        title('Beam [1MW] => target neutrons /s/m^{3} ')
        
    end
    
end
%%

if isfield(output,'ejected_CX_dt')&par.CALCULATE_CX
    
    Rpos=squeeze(output.x_ej(output.ejected_CX_dt,1));
    Zpos=squeeze(output.x_ej(output.ejected_CX_dt,2));
    pow_values=squeeze(output.Ekin_ej(output.ejected_CX_dt))/DT_SIM_RECORD;
    
    Rbins = scale_X(1:25:end)+R0+0.02;
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
        plot_CU_TOK_background
        ylim([-0.6 0.6])
        xlim([min(scale_X)+R0+0.02 1.21])
        
        set(gca,'fontsize',20)
        title('CU NBI-CX (1MW) @wall W/m^{2} ')
        
    end
    
end

%%

if isfield(output,'ejected_wall_dt')
    EJ_ION_POP_DT=output.ejected_wall_dt&~output.ejected_CX_dt;
    
    Rpos=squeeze(output.x_ej(EJ_ION_POP_DT,1));
    Zpos=squeeze(output.x_ej(EJ_ION_POP_DT,2));
    pow_values=squeeze(output.Ekin_ej(EJ_ION_POP_DT))/DT_SIM_RECORD;
    
    Rbins = scale_X(1:25:end)+R0+0.02;
    Zbins = scale_Z(1:30:end);
    Rvals = Rbins(1:end-1)+0.5*(Rbins(2)-Rbins(1));
    Zvals = Zbins(1:end-1)+0.5*(Zbins(2)-Zbins(1));
    
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
        plot_CU_TOK_background
        ylim([-0.6 0.6])
        xlim([min(scale_X)+R0+0.02 1.21])
        
        set(gca,'fontsize',20)
        title('CU NBI-ions-losses (1MW) @wall W/m^{2} ')
        
    end
    
end


try
    save(par.FILENAME,'-append','summary')
    save(par.FILENAME,'-append','pfc_losses')
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
        FILENAME_PREC=[folder_save_string par.FILENAME(1:end-4) '_prec.mat'];
    else
        FILENAME_PREC=['./' par.FILENAME(1:end-4) '_prec.mat'];
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
