
USE_GC_PSIPOS_FOR_PROFILES=0;
ANALYZE_BB_INTERACTION=0;
SAVE_NFAST_PROF_FOR_BB_INTERACTION=1;


if ~exist('psi_norm1_XZsmall_map')
    load('../../data_tokamak/physics_constants.mat');
    load('../../data_tokamak/XZsmall_fields_tokamak_pre_collapse.mat');
    load('../../data_tokamak/motions_map_dimensions.mat');
end
try
    load('../../data_tokamak/flux_geometry.mat');
catch
    load('./data_tokamak/flux_geometry.mat');
end

% larger interval because of Beam-beam data calculations
Ekin_bins_ndd=0:2e3:240*1e3;
Ekin_values_ndd=Ekin_bins_ndd(1:end-1)+0.5.*(Ekin_bins_ndd(2)-Ekin_bins_ndd(1));

output.psi=interp2(scale_X+R0,scale_Z,psi_XZsmall_map',squeeze(output.x(:,1,:)),squeeze((output.x(:,2,:))));
psi_max_orbits=max(max(output.psi));
psi_min_orbits=min(min(output.psi));
if psi_scale(end)>psi_scale(1)
    psi_scale_extend=[psi_scale psi_max_orbits];
else
    psi_scale_extend=[psi_scale psi_min_orbits];
end
rho_tor_extend=[rho_tor_scale (rho_tor_scale(end)-rho_tor_scale(end-1))*20+rho_tor_scale(end)];

% ITER like case
if psi_scale(end)>psi_scale(1)
    rho_tor_extension=interp1(psi_scale_extend,rho_tor_extend,psi_max_orbits)-1
    volume_flux_extension=volume_flux(end)-interp1(psi_scale,volume_flux,psi_max_orbits-psi_scale(end)+psi_scale(1))
else
    % MAST-U or AUG case
    rho_tor_extension=interp1(psi_scale_extend,rho_tor_extend,psi_min_orbits)-1
    volume_flux_extension=volume_flux(end)-volume_flux(end-20)
end

% arbitraray extension of scale to avoid nans
psi_max_scale=max(psi_scale);
psi_min_scale=min(psi_scale);


PLOT_DETAILS=0;
PLOT_TORQUE=0;


if ~exist('birth_matrix')
    try
        if ~isfield(output,'birth_matrix')
            error('you need to run summary script to generate birth_matrix first!')
        else
            birth_matrix=output.birth_matrix;
        end
    catch
        error('load an output file first!')
    end
end

rho_tor_scale_adjust=[0 rho_tor_scale(3:3:10) rho_tor_scale(20:10:end-2) rho_tor_scale(end)];
volume_flux_adjust=[0 volume_flux(3:3:10) volume_flux(20:10:end-2) volume_flux(end)];
% rho_tor_scale_adjust=rho_tor_scale(1:8:end);
% volume_flux_adjust=volume_flux(1:8:end);
% dv_flux(1:11)=volume_flux_adjust(2:12)-volume_flux_adjust(1:11);
% dv_flux(12:12:end)=volume_flux_adjust(24:12:end)-volume_flux_adjust(12:12:end-12);


% extending the scales to SOL : simple linear extrapolation
% ITER like case
if psi_scale(end)>psi_scale(1)
    psi_scale_ext=[psi_scale psi_max_orbits];
    psi_scale_norm=psi_scale_ext-psi_max_scale;
    psi_scale_norm=1-psi_scale_norm/psi_scale_norm(1);
else
    % MAST-U or AUG case
    psi_scale_ext=[psi_scale psi_min_orbits];
    psi_scale_norm=psi_scale_ext-psi_min_scale;
    psi_scale_norm=1-psi_scale_norm/psi_scale_norm(1);
end

rho_tor_scale_adjust=[rho_tor_scale_adjust 1+rho_tor_extension];
volume_flux_adjust=[volume_flux_adjust volume_flux_adjust(end)+volume_flux_extension];
rho_tor_scale_ext=[rho_tor_scale 1+rho_tor_extension];
volume_flux_ext=[volume_flux volume_flux(end)+volume_flux_extension];

dv_flux=volume_flux_adjust(2:end)-volume_flux_adjust(1:end-1);
summary.volume_flux=volume_flux_adjust;

% summary.rho_tor_scale_adjust=rho_tor_scale_adjust;
% summary.volume_flux_adjust=volume_flux_adjust;
summary.dv_flux=dv_flux;


% psi_scale_norm=1-psi_scale_norm/psi_scale_norm(1);

%%
output.psi_final=ejected*0;
output.rhot_final=ejected*0;

output.psi_final(~ejected|output.ejected_sd)=interp2(scale_X+R0,scale_Z,psi_XZsmall_map',input.x_end(~ejected|output.ejected_sd,1),input.x_end(~ejected|output.ejected_sd,2));
output.psi_ini=interp2(scale_X+R0,scale_Z,psi_XZsmall_map',input.x_ini(:,1),input.x_ini(:,2));
output.psi_ini_nej=interp2(scale_X+R0,scale_Z,psi_XZsmall_map',input.x_ini(~ejected,1),input.x_ini(~ejected,2));
output.psi_gc=interp2(scale_X+R0,scale_Z,psi_XZsmall_map',squeeze(output.x_gc(:,1,:)),squeeze((output.x_gc(:,2,:))));
output.psi=interp2(scale_X+R0,scale_Z,psi_XZsmall_map',squeeze(output.x(:,1,:)),squeeze((output.x(:,2,:))));

output.rhot_final=interp1(psi_scale_ext,rho_tor_scale_ext,output.psi_final);
output.rhot_ini=interp1(psi_scale_ext,rho_tor_scale_ext,output.psi_ini);
output.rhot_ini_nej=interp1(psi_scale_ext,rho_tor_scale_ext,output.psi_ini_nej);
% output.rhot=interp1(psi_scale,rho_tor_scale,output.psi_gc(:,:));
if ~USE_GC_PSIPOS_FOR_PROFILES
    try
        output.rhot=interp1(psi_scale_ext,rho_tor_scale_ext,output.psi(:,:));
    catch
        output.rhot=interp1(psi_scale_ext,rho_tor_scale_ext,output.psi_gc(:,:));
    end
else
    output.rhot=interp1(psi_scale_ext,rho_tor_scale_ext,output.psi_gc(:,:));
end
% output.rhot=interp1(psi_scale,rho_tor_scale,output.psi(:,:));

NB_TS_AVG=par.NB_TS_RECORDED-1;

%%
% PROFILES_POP=~ejected|output.ejected_sd;
PROFILES_POP=~ejected;
% ej_before=logical(output.time_stamp_loss<=(par.NB_TIME_STAMPS-par.NB_TS_RECORDED))&ejected;
% BORN_MARKERS=~logical(birth_matrix(1,:)');

% we can take all markers : energy evolution stops at loss
PROFILES_ETOT_POP=logical(ejected*0+1);
% PROFILES_ETOT_POP=logical(~ejected.*BORN_MARKERS);

if isfield(output,'ndd_val_BP')
    values_markers_ndd=squeeze(output.ndd_val_BP(PROFILES_POP,end));
elseif isfield(output,'ndd_val')
    values_markers_ndd=squeeze(output.ndd_val(PROFILES_POP,end));
end
if isfield(output,'cum_tor_ang_mom')
    values_markers_torq_end=squeeze(output.cum_tor_ang_mom(:,end));
    values_markers_torq_ini=squeeze(output.cum_tor_ang_mom(:,1));
end
values_markers_Pdep_e=squeeze(output.Pdep_e(PROFILES_POP,end));
values_markers_Pdep_i=squeeze(output.Pdep_i(PROFILES_POP,end));
if par.IMPURITY_ZAVG>0
    values_markers_Pdep_imp=squeeze(output.Pdep_imp(PROFILES_POP,end));
end
values_markers_Edep_tot=squeeze(output.Ekin(PROFILES_ETOT_POP,1)-output.Ekin(PROFILES_ETOT_POP,end));
rho_tor_bins=rho_tor_scale_adjust(1:end-1)+0.5*(rho_tor_scale_adjust(2:end)-rho_tor_scale_adjust(1:end-1));

%arbitrary adjustment
rho_tor_bins(end)=1+(rho_tor_bins(end)-1)*0.5;

[NMARKERS,BINVALS] = histc(output.rhot(PROFILES_POP,end),rho_tor_scale_adjust);
[NMARKERS_TOT,BINVALS_TOT] = histc(output.rhot(PROFILES_ETOT_POP,end),rho_tor_scale_adjust);

values_markers_Pdep_th=squeeze(output.Edep_th(output.ejected_sd_dt));
[NMARKERS_TH,BINVALS_TH] = histc(mean(output.rhot(output.ejected_sd_dt,:),2),rho_tor_scale_adjust);

values_markers_Ploss_CX=squeeze(output.Ekin_ej(output.ejected_CX_dt));
[NMARKERS_CX,BINVALS_CX] = histc(mean(output.rhot(output.ejected_CX_dt,:),2),rho_tor_scale_adjust);


ndd_profile=rho_tor_bins*0;
if isfield(output,'ndd_val')
    for ii=1:length(rho_tor_bins)
        ndd_profile(ii)=sum(values_markers_ndd(BINVALS==ii));
    end
end
torque_profile=rho_tor_bins*0;
if isfield(output,'cum_tor_ang_mom')
    for ii=1:length(rho_tor_bins)
        torque_profile(ii)=sum(values_markers_torq_end(BINVALS==ii))-sum(values_markers_torq_ini(BINVALS==ii));
    end
end
torque_profile=torque_profile/summary.DT_SS; % in N.m
summary.torque_profile=torque_profile;
summary.torque_density_profile=torque_profile./dv_flux;

pdep_e_profile=rho_tor_bins*0;
for ii=1:length(rho_tor_bins)
pdep_e_profile(ii)=sum(values_markers_Pdep_e(BINVALS==ii));
end

pdep_i_profile=rho_tor_bins*0;
for ii=1:length(rho_tor_bins)
pdep_i_profile(ii)=sum(values_markers_Pdep_i(BINVALS==ii));
end

pdep_imp_profile=rho_tor_bins*0;
if par.IMPURITY_ZAVG>0
    for ii=1:length(rho_tor_bins)
        pdep_imp_profile(ii)=sum(values_markers_Pdep_imp(BINVALS==ii));
    end
end

pdep_th_profile=rho_tor_bins*0;
for ii=1:length(rho_tor_bins)
pdep_th_profile(ii)=sum(values_markers_Pdep_th(BINVALS_TH==ii));
end
pdep_th_profile=pdep_th_profile*par.MARKER_WEIGHT*eV/DT_SIM_RECORD;

ploss_cx_profile=rho_tor_bins*0;
for ii=1:length(rho_tor_bins)
ploss_cx_profile(ii)=sum(values_markers_Ploss_CX(BINVALS_CX==ii));
end
ploss_cx_profile=ploss_cx_profile*par.MARKER_WEIGHT*eV/DT_SIM_RECORD;

Edep_tot_profile=rho_tor_bins*0;
for ii=1:length(rho_tor_bins)
Edep_tot_profile(ii)=sum(values_markers_Edep_tot(BINVALS_TOT==ii));
end
% Pdep_tot_profile=Edep_tot_profile*par.INPUT_POWER/sum(Edep_tot_profile);
Pdep_tot_profile=Edep_tot_profile*par.MARKER_WEIGHT*eV/DT_SIM_RECORD;

ndd_Ekin_dist=Ekin_bins_ndd*0;

if PLOT_NDD
    output.Ekin(ejected,end)=nan;
    [EKIN_DIST_vals EKIN_DIST_POP]=histc(output.Ekin(:,end),Ekin_bins_ndd);
    for ek=1:length(EKIN_DIST_vals)
        ndd_Ekin_dist(ek)=sum(output.ndd_val_BP(EKIN_DIST_POP==ek,end));
    end
    ndd_Ekin_dist_TS=ndd_Ekin_dist;
end

nb1D_markers=par.MARKER_WEIGHT*NMARKERS(1:end-1);

%
BORN_MARKERS=logical(birth_matrix(end,:)');
POP_TS=logical(~ejected.*BORN_MARKERS);

%%
if NB_TS_AVG>0
    for TS=1:NB_TS_AVG
        
        n_ej=logical(output.time_stamp_loss>=par.NB_TIME_STAMPS-TS)|(isnan(output.time_stamp_loss));
        % keeping thermalized markers in density profile
        % n_ej=n_ej|~output.ejected_sd;
        
        BORN_MARKERS=logical(birth_matrix(end-TS,:)');
        POP_TS=logical(n_ej.*BORN_MARKERS);
        disp(['at TS=' num2str(TS) ' POP_TS ' num2str(length(find(POP_TS)))]);
        

        if PLOT_NDD
            output.Ekin(~POP_TS,end-TS)=nan;
            [EKIN_DIST_vals EKIN_DIST_POP]=histc(output.Ekin(:,end-TS),Ekin_bins_ndd);
            for ek=1:length(EKIN_DIST_vals)
                ndd_Ekin_dist_TS(ek)=sum(output.ndd_val_BP(EKIN_DIST_POP==ek,end-TS));
            end
            ndd_Ekin_dist=ndd_Ekin_dist+ndd_Ekin_dist_TS;
        end
        
        [NMARKERS,BINVALS] = histc(output.rhot(POP_TS,end-TS),rho_tor_scale_adjust);
        
        if PLOT_NDD
            nb1D_markers=nb1D_markers+par.MARKER_WEIGHT*NMARKERS(1:end-1);
            if isfield(output,'ndd_val_BP')
                values_markers_ndd=squeeze(output.ndd_val_BP(POP_TS,end-TS));
            elseif isfield(output,'ndd_val')
                values_markers_ndd=squeeze(output.ndd_val(POP_TS,end-TS));
           end
        end
        values_markers_Pdep_e=squeeze(output.Pdep_e(POP_TS,end-TS));
        values_markers_Pdep_i=squeeze(output.Pdep_i(POP_TS,end-TS));
        if par.IMPURITY_ZAVG>0
            values_markers_Pdep_imp=squeeze(output.Pdep_imp(POP_TS,end-TS));
        end
%         values_markers_Edep_tot=squeeze(output.Ekin(POP_TS,1)-output.Ekin(POP_TS,end-TS));
        
        if PLOT_NDD
            ndd_profile_TS=rho_tor_bins*0;
            if isfield(output,'ndd_val')
                for ii=1:length(rho_tor_bins)
                    ndd_profile_TS(ii)=sum(values_markers_ndd(BINVALS==ii));
                end
            end
        end
        
        pdep_e_profile_TS=rho_tor_bins*0;
        for ii=1:length(rho_tor_bins)
            pdep_e_profile_TS(ii)=sum(values_markers_Pdep_e(BINVALS==ii));
        end
        
        pdep_i_profile_TS=rho_tor_bins*0;
        for ii=1:length(rho_tor_bins)
            pdep_i_profile_TS(ii)=sum(values_markers_Pdep_i(BINVALS==ii));
        end
        
        if par.IMPURITY_ZAVG>0
            pdep_imp_profile_TS=rho_tor_bins*0;
            for ii=1:length(rho_tor_bins)
                pdep_imp_profile_TS(ii)=sum(values_markers_Pdep_imp(BINVALS==ii));
            end
        end
%         Edep_tot_profile_TS=rho_tor_bins*0;
%         for ii=1:length(rho_tor_bins)
%             Edep_tot_profile_TS(ii)=sum(values_markers_Edep_tot(BINVALS==ii));
%         end

%         Edep_tot_profile=Edep_tot_profile+Edep_tot_profile_TS;
        if PLOT_NDD
            if isfield(output,'ndd_val')
                ndd_profile=ndd_profile+ndd_profile_TS;
            end
        end
        pdep_e_profile=pdep_e_profile+pdep_e_profile_TS;
        pdep_i_profile=pdep_i_profile+pdep_i_profile_TS;
        if par.IMPURITY_ZAVG>0
            pdep_imp_profile=pdep_imp_profile+pdep_imp_profile_TS;
        end
%         figure(21)
%         hold on
%         plot(pdep_e_profile_TS)
        
    end
    %     Edep_tot_profile=Edep_tot_profile/(NB_TS_AVG+1);
    if PLOT_NDD
        ndd_Ekin_dist=ndd_Ekin_dist/(NB_TS_AVG+1);
    end
    nb1D_markers=nb1D_markers/(NB_TS_AVG+1);
    ndd_profile=ndd_profile/(NB_TS_AVG+1);
    pdep_e_profile=pdep_e_profile/(NB_TS_AVG+1);
    pdep_i_profile=pdep_i_profile/(NB_TS_AVG+1);
    if par.IMPURITY_ZAVG>0
        pdep_imp_profile=pdep_imp_profile/(NB_TS_AVG+1);
    end
%     Pdep_tot_profile=Edep_tot_profile*par.INPUT_POWER/(eV*par.MARKER_WEIGHT*sum(Edep_tot_profile)/(par.time_scale(end)-par.time_scale(end-NB_TS_AVG)));
%     Pdep_tot_profile=Edep_tot_profile*par.INPUT_POWER/sum(Edep_tot_profile);

end
% Pdep_tot_profile=Edep_tot_profile*par.MARKER_WEIGHT*eV/DT_SIM_RECORD;

% increased pdep_i to take into account remaining power after v~1.5vthi
% pdep_i_profile=pdep_i_profile*(par.INPUT_POWER-sum(pdep_e_profile))/(sum(pdep_i_profile));

% if everything that is not given to electrons goes to ions
% pdep_i_profile_according_to_tot=Pdep_tot_profile-pdep_e_profile;

neutron_rate=sum(ndd_profile);
Power_to_electrons=sum(pdep_e_profile);
Power_to_ions=sum(pdep_i_profile);
Power_to_thermal=sum(pdep_th_profile);
Power_loss_CX=sum(ploss_cx_profile);
Power_tot_eiimp=sum(Pdep_tot_profile);
Power_tot_imp=sum(pdep_imp_profile);

disp(['neutron_rate [10^14] for ' num2str(par.INPUT_POWER) ' MW of NBI = ' num2str(1e-14*neutron_rate)])
disp('---------------------------------------------------------------')

disp(['Power to plasma (dEkin) = ' num2str(1e-3*Power_tot_eiimp)])
disp(['Power to electrons  = ' num2str(1e-3*Power_to_electrons)])
disp(['Power to ions       = ' num2str(1e-3*Power_to_ions)])
disp(['Power to impurity   = ' num2str(1e-3*Power_tot_imp)])
disp(['Power thermal       = ' num2str(1e-3*mean(summary.Ploss_th_evol))])
disp(['Power CX loss       = ' num2str(1e-3*mean(summary.Ploss_CX_evol))])
disp(['Power ion loss      = ' num2str(1e-3*mean(summary.Ploss_ion_evol))])
disp(['Power loss + plasma (e,i,imp) + th = ' num2str(1e-3*(summary.POWER_ION_WALL_LOSS+summary.POWER_CX_LOSS+Power_to_electrons+Power_to_ions+Power_tot_imp+Power_to_thermal))])

Pdep_tot_profile=Pdep_tot_profile./dv_flux;
pdep_e_profile=pdep_e_profile./dv_flux;
pdep_i_profile=pdep_i_profile./dv_flux;
pdep_imp_profile=pdep_imp_profile./dv_flux;
pdep_th_profile=pdep_th_profile./dv_flux;
ndd_profile=ndd_profile./dv_flux;

density_markers=nb1D_markers./dv_flux';
% density_markers=par.MARKER_WEIGHT*NMARKERS(1:end-1)./dv_flux';

summary.neutron_rate=neutron_rate;
summary.Power_to_electrons=Power_to_electrons;
summary.Power_to_ions=Power_to_ions;
summary.Power_to_thermal=Power_to_thermal;

summary.ndd_Ekin_dist=ndd_Ekin_dist(1:end-1);
summary.Ekin_values_ndd=Ekin_values_ndd;

summary.Pdep_tot_profile=Pdep_tot_profile;
summary.pdep_e_profile=pdep_e_profile;
summary.pdep_i_profile=pdep_i_profile;
summary.pdep_imp_profile=pdep_imp_profile;
summary.pdep_th_profile=pdep_th_profile;
summary.ndd_profile=ndd_profile;
summary.density_markers=density_markers;
summary.rho_tor_bins=rho_tor_bins;

if SAVE_NFAST_PROF_FOR_BB_INTERACTION
    load('../../data_tokamak/pressure_profile.mat');
    
    nfast_prof_rhot=interp1(summary.rho_tor_bins,summary.density_markers ,  rho_tor_extend,'pchip');
    nfast_prof_rhot(rho_tor_extend>summary.rho_tor_bins(end))=0;
    nfast_prof=interp1(psi_scale_norm,nfast_prof_rhot , psi_norm,'pchip');
    save('../../data_tokamak/pressure_profile.mat','-append','nfast_prof');
end

try
    save(FNAME,'-append','summary')
% 
% save(par.FILENAME,'-append','summary')
catch
    warning('could not save summary to file!!!')
end


if ANALYZE_BB_INTERACTION
    calculate_BeamBeam_nDD;
end

%%
figure(8)
subplot(2,1,1)
grid on
set(gca,'fontsize',20)
xlabel('\rho_t')
ylabel('Pdep [W/m^{3}]')
hold on

plot(rho_tor_bins,pdep_e_profile,'linewidth',2)
plot(rho_tor_bins,pdep_i_profile,'linewidth',2)
plot(rho_tor_bins,pdep_imp_profile,'linewidth',2)
% plot(rho_tor_bins,pdep_e_profile+pdep_i_profile);
% plot(rho_tor_bins,Pdep_tot_profile);
legend('Pe','Pi+Pth')
xlim([0 1.1])


%%
if PLOT_NDD
    
%     figure
%     plot(Ekin_values_ndd,  ndd_Ekin_dist(1:end-1))
%     grid on
%     xlabel('Ekin [eV]')
%     ylabel('nDD')
    
    figure
    grid on
    set(gca,'fontsize',20)
    xlabel('\rho _{pol}')
    ylabel('BP neutron rate [m^{-3} s^{-1}]')
    hold on
    
    if psi_scale(end)>psi_scale(1)
        psi_norm_scale=(psi_scale - psi_scale(1))/psi_global;
    else
        psi_norm_scale=-(psi_scale - psi_scale(1))/psi_global;
    end
    summary.rho_pol_bins=interp1(rho_tor_scale,psi_norm_scale,summary.rho_tor_bins,'pchip');
    summary.rho_pol_bins=sqrt(summary.rho_pol_bins);

    plot(summary.rho_pol_bins,ndd_profile,'linewidth',2)
    xlim([0 1.1])

end

%%

figure(8)
subplot(2,1,2)
grid on
set(gca,'fontsize',20)
xlabel('\rho_t')
ylabel('density of markers')
hold on

plot(rho_tor_bins,density_markers,'linewidth',2)

xlim([0 1.1])


%%
if PLOT_TORQUE
    figure(9)
    grid on
    set(gca,'fontsize',20)
    xlabel('\rho')
    ylabel('Nm^{-2}')
    hold on
    
    plot(rho_tor_bins,summary.torque_density_profile,'linewidth',2)
    
    xlim([0 1.1])
end

%%
%Edep_th

if isfield(output,'Edep_th')
    if PLOT_DETAILS

    figure
    plot(output.rhot(find(output.Edep_th>0&output.ejected_sd_dt)),output.Edep_th(find(output.Edep_th>0&output.ejected_sd_dt)),'k.')
    xlim([0 1.1])

    end
    
end

%%
if PLOT_DETAILS
    
    figure
    grid on
    set(gca,'fontsize',20)
    xlabel('\rho')
    ylabel('(#ions)/d\rho')
    hold on
    
    % plot(rho_tor_bins,values_profile)
    [NMARKERS,BINVALS] = histc(output.rhot_final,rho_tor_scale_adjust);
    plot(rho_tor_bins,par.MARKER_WEIGHT*NMARKERS(1:end-1))
    
    
    [NMARKERS,BINVALS] = histc(output.rhot_ini,rho_tor_scale_adjust);
    plot(rho_tor_bins,par.MARKER_WEIGHT*NMARKERS(1:end-1))
    
    
    [NMARKERS,BINVALS] = histc(output.rhot_ini_nej,rho_tor_scale_adjust);
    plot(rho_tor_bins,par.MARKER_WEIGHT*NMARKERS(1:end-1),'linestyle','--')
    
    legend('final (without ejected)','initial (without ejected)','initial (with ejected)')

end

%%

% if PLOT_DETAILS
%     
%     Rbins=0.5:0.025:1.2;
%     Zbins=-0.7:0.025:0.7;
%     
%     figure;
%     hold on
%     % imagesc(Rbins,Zbins,hist2d([input.x_ini(~ejected,2) input.x_ini(~ejected,1)],Zbins,Rbins));
%     
%     imagesc(Rbins,Zbins,hist2d([input.x_end(~ejected,2) input.x_end(~ejected,1)],Zbins,Rbins));
%     plot_CU_TOK_background;
%     
%     
% end

