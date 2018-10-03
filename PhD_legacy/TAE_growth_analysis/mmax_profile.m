run('reset_data_analysis_environment.m')
rho_scale=radial_r_value_flux/max(radial_r_value_flux);
Btot_PR_map=sqrt(Bpol_PR_map.^2+Btor_PR_map.^2);
Bavg_radial=mean(Btot_PR_map(1:end-1,:),1);
PSI_BIN_SIZE=12

load('initialG_DT_MB_all_pre_collapse.mat')
load('initialG_DT_MB_all_precession_stats.mat')

psi_scale_binned=(1+0.5*PSI_BIN_SIZE:PSI_BIN_SIZE:Nradial);

for psi_pos=1:length(psi_scale_binned)
    Fpsi_pop=(alphas_psi>=psi_scale_binned(psi_pos)-0.5*PSI_BIN_SIZE).*(alphas_psi<=psi_scale_binned(psi_pos)+0.5*PSI_BIN_SIZE);
    POP_SLICE=find(Fpsi_pop);
    orbit_width_binned(psi_pos)=mean(delta_r_avg(POP_SLICE))
end
orbit_width_profile=zeros(1,Nradial);
orbit_width_profile(1:Nradial-PSI_BIN_SIZE)=interp1(psi_scale_binned,orbit_width_binned,1:Nradial-PSI_BIN_SIZE);
rho_scale_ITER=rho_scale;

%%
figure(3)
set(gca,'FontSize',26)
hold on; grid on
% plot(rho_scale,m_max_profile,'b','LineWidth',3)
plot(rho_scale,orbit_width_profile,'r','LineWidth',3)
xlabel('\rho')
ylabel('m_{max}')
legend('JET','ITER')

%%
Ni_XZsmall_map=interp1(1:257,Ne_profile,psi_norm_XZsmall_map);
omega_ci_profile=eV*Bavg_radial/mDT;
rhoLi_profile=sqrt(Te_profile/mDT)./omega_ci_profile;
m_max_profile_ITER=radial_r_value_flux./(orbit_width_profile);

vA_map=(Btot_XZ_map./sqrt(mu0*Ni_XZsmall_map*mDT));
omegaTAE_map=vA_map./(2*(q_initial_XZsmall_map).*R0);
omegaA_map=vA_map.*(6-6./q_initial_XZsmall_map)./R0;

%%
figure(2)
set(gca,'FontSize',26)
hold on; grid on
% plot(rho_scale,m_max_profile,'b','LineWidth',3)
plot(rho_scale_ITER,m_max_profile_ITER,'r','LineWidth',3)
xlabel('\rho')
ylabel('m_{max}')
ylim([0 220])
legend('JET','ITER')

%%
figure(1)
title('\omega_{TAE}')
hold on; grid on

contour(scale_X+R0,-scale_Z,omegaTAE_map','LineWidth',3);
colorbar

% contour(X_scale+R0,-Z_scale,psi_XZ_map',psi_scale(2:28:end),'k')

% contour(scale_X+R0,-scale_Z,q_initial_XZsmall_map',((1:30)+0.5)/5,'k--');

contour(X_scale+R0,-Z_scale,psi_XZ_map',[0 0],'k','linewidth',4)
set(gca,'FontSize',26)
xlabel('R (m)')
ylabel('Z (m)')
% colorbar
axis square
% title('Contour of poloidal flux surfaces of AUG')

nTAE=6
q_TAE_m=(nTAE)/nTAE
psi_TAE=interp1(q_initial_profile,psi_scale,q_TAE_m)
contour(X_scale+R0,-Z_scale,psi_XZ_map',[psi_TAE psi_TAE],'color',[0.1 0.1 0.5],'linewidth',2)

q_TAE_m=(nTAE+1)/nTAE
psi_TAE=interp1(q_initial_profile,psi_scale,q_TAE_m)
contour(X_scale+R0,-Z_scale,psi_XZ_map',[psi_TAE psi_TAE],'color',[0.5 0.1 0.1],'linewidth',2)

contour(scale_X+R0,-scale_Z,omegaTAE_map','LineWidth',3);

axis equal
xlabel('X (m)')
ylabel('Z (m)')