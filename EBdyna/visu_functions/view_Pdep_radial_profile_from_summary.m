

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








