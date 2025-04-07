

PLOT_RHO_POLOIDAL=0


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

if psi_scale(end)>psi_scale(1)
    psi_norm_scale=(psi_scale - psi_scale(1))/psi_global;
else
    psi_norm_scale=-(psi_scale - psi_scale(1))/psi_global;
end
summary.rho_pol_bins=interp1(rho_tor_scale,psi_norm_scale,summary.rho_tor_bins);
summary.rho_pol_bins=sqrt(summary.rho_pol_bins);

% if ~exist('birth_matrix')
%     try
%         if ~isfield(output,'birth_matrix')
%             error('you need to run summary script to generate birth_matrix first!')
%         else
%             birth_matrix=output.birth_matrix;
%         end
%     catch
%         error('load an output file first!')
%     end
% end



PLOT_NDD=0;


%%



%%
figure(9)
grid on
set(gca,'fontsize',22)
xlabel('\rho')
ylabel('Pdep [W/m^3]')
hold on

% plot(summary.rho_tor_bins,summary.pdep_e_profile+summary.pdep_i_profile+summary.pdep_th_profile,'linewidth',2)
if PLOT_RHO_POLOIDAL
    plot(summary.rho_pol_bins,summary.pdep_e_profile+summary.pdep_i_profile+summary.pdep_th_profile,'linewidth',2)
else
    plot(summary.rho_tor_bins,summary.pdep_e_profile+summary.pdep_i_profile+summary.pdep_th_profile,'linewidth',2)
%     plot(summary.rho_tor_bins,summary.pdep_e_profile+summary.pdep_i_profile+summary.pdep_th_profile,'color',[1 0.6 0.1],'linestyle','--','linewidth',2)
    xlabel('\rho_t')
end
% plot(summary.rho_tor_bins,summary.pdep_i_profile)
% plot(rho_tor_bins,pdep_e_profile+pdep_i_profile);
% plot(rho_tor_bins,Pdep_tot_profile);
% legend('Pe','Pi')
xlim([0.05 1.05])


%%
if PLOT_NDD
    
    figure
    grid on
    set(gca,'fontsize',20)
    xlabel('\rho')
    ylabel('neutron rate [m^{-3} s^{-1}]')
    hold on
    
    if PLOT_RHO_POLOIDAL
        plot(summary.rho_pol_bins,summary.ndd_profile)
    else
        plot(summary.rho_tor_bins,summary.ndd_profile)
    end
    xlim([0 1.1])
    
end

%%

figure(7)
grid on
set(gca,'fontsize',20)
xlabel('\rho')
ylabel('density of Fast Ions')
hold on

%
if PLOT_RHO_POLOIDAL
    plot(summary.rho_pol_bins,summary.density_markers)
else
    plot(summary.rho_tor_bins,summary.density_markers)
end
xlim([0 1.1])

%%
figure(8)
subplot(2,1,2)
grid on
set(gca,'fontsize',16)
xlabel('\rho_t')
ylabel('density of markers')
hold on

% plot(summary.rho_tor_bins,summary.density_markers,'linewidth',2)
if PLOT_RHO_POLOIDAL
    plot(summary.rho_pol_bins,summary.density_markers,'linewidth',2)
    xlabel('\rho_p')
else
    pl=plot(summary.rho_tor_bins,summary.density_markers,'linewidth',2)
    xlabel('\rho_t')
end

xlim([0 1.1])




figure(8)
subplot(2,1,1)
grid on
set(gca,'fontsize',16)
xlabel('\rho')
ylabel('Pdep [W/m^{3}]')
hold on
% 
% plot(summary.rho_tor_bins,summary.pdep_e_profile,'linewidth',2)
% plot(summary.rho_tor_bins,summary.pdep_i_profile+summary.pdep_th_profile,'linewidth',2)
if PLOT_RHO_POLOIDAL
    plot(summary.rho_pol_bins,summary.pdep_e_profile,'linewidth',2)
    plot(summary.rho_pol_bins,summary.pdep_i_profile+summary.pdep_th_profile,'linewidth',2)
    plot(summary.rho_pol_bins,summary.pdep_e_profile+summary.pdep_i_profile+summary.pdep_th_profile,'k--','linewidth',2)
    xlabel('\rho_p')
else
    pl2=plot(summary.rho_tor_bins,summary.pdep_e_profile,'linewidth',2)
    pl2.Color=pl.Color;
    pl3=plot(summary.rho_tor_bins,summary.pdep_i_profile+summary.pdep_th_profile,'linewidth',2,'linestyle','--')
    pl3.Color=pl.Color;
    xlabel('\rho_t')
end

xlim([0 1.1])

% subplot(2,1,2)
% grid on
% set(gca,'fontsize',20)
% xlabel('\rho')
% ylabel('density of NBI [10^{19} m^{-3}]')
% hold on
% 
% plot(summary.rho_tor_bins,summary.density_markers*1e-19,'linewidth',2)
% 
% xlim([0 1.1])




