
FILENAME=strcat(DATA_FOLDER,'tokamak_PR_map.mat')
save (FILENAME,'psiH_PR_map','psi_PR_map','NP','Nradial','R0','a','elongation','theta_PR_map','Rpos_PR_map','BX_PR_map','BZ_PR_map','X_axis','Z_axis','radial_r_value_flux','X_scale','Z_scale','finesse_data_X','finesse_data_Z');

FILENAME=strcat(DATA_FOLDER,'pressure_profile.mat')
save (FILENAME,'P_initial_profile','plasma_beta_tot','P0','Ne_profile','Te_profile','Ne0','Te0');

Iaxis=I_flux(end)
FILENAME=strcat(DATA_FOLDER,'plasma_current.mat')
save (FILENAME,'Iaxis','I_flat_flux','I_flux','J_PHI_profile','J_PHI_PR_map','JPHI_XZ_map','F2_profile','F2_prime')

psiH_pol_profile=mean(psiH_PR_map(1:NP-1,:),1);
psi_star_initial_profile=psi_star_Nradial_avg;
psi_pol_initial_profile=psi_profile;
FILENAME=strcat(DATA_FOLDER,'psi_profiles.mat')
save (FILENAME,'psi_star_initial_profile','psiH_pol_profile','psi_pol_initial_profile');

q_initial_profile=q_profile;
FILENAME=strcat(DATA_FOLDER,'q_profile.mat')
save (FILENAME,'q_initial_profile','psi_q1','r_value_q1_mean','psi_rank_q1','q_XZ_map');

BpolX_initial_map=BR_XZ_map;
BpolZ_initial_map=BZ_XZ_map;
FILENAME=strcat(DATA_FOLDER,'B_fields.mat')
save (FILENAME,'BpolX_initial_map','BpolZ_initial_map','Bpol_PR_map','Btor_PR_map','Bstar_PR_map','Bphi0','BHpol_PR_map','BHpol_X_map','BHpol_Z_map','BHpol_X_PR_map','BHpol_Z_PR_map','Bphi_XZ_map');

xi_profile=xi_psi_fit; % NPSI points
%volume_flux(2)=0.5*(volume_flux(3)+volume_flux(1));
FILENAME=strcat(DATA_FOLDER,'flux_geometry.mat')
save (FILENAME,'psi_scale','tor_flux_scale','rho_tor_scale','theta_ki_psi_map','X_PR_map','psiH_XZ_map','psi_XZ_map','volume_tor_diff','scale_tor_x',...
'Z_PR_map','radial_XZ_map','theta_XZ_map','dr_PR_map','dr_X_PR_map','dr_Z_PR_map','cos_ki_PR_map',...
'dl_PR_map','dl_X_PR_map','dl_Z_PR_map','psi_global','xi_profile','volume_flux','surf_flux','dist_surf_PR_map','radial_r_value_flux');

Z_psi_fit_up=Z_psi_fit;
FILENAME=strcat(DATA_FOLDER,'volume_flux_geometry.mat')
save (FILENAME,'tor_flux_profile','volume_tor','volume_tor_diff','Z_psi_fit_up','Z_psi_fit_down','X1_Nradial','X2_Nradial');

FILENAME=strcat(DATA_FOLDER,'grad_flux_geometry.mat')
save (FILENAME,'grad_theta_PR_map_X','grad_theta_PR_map_Z','g_polX_XZ_map','g_polZ_XZ_map');