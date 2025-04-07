if ADD_FIESTA_TO_FILENAMES==1
    FILENAME=strcat(DATA_FOLDER,'tokamak_PR_map_FIESTA.mat')
else
    FILENAME=strcat(DATA_FOLDER,'tokamak_PR_map.mat')
end
save (FILENAME,'Nradial','R0','a','elongation','X_axis','Z_axis','radial_r_value_flux','X_scale','Z_scale');

if ADD_FIESTA_TO_FILENAMES==1
    FILENAME=strcat(DATA_FOLDER,'pressure_profile_FIESTA.mat')
else
    FILENAME=strcat(DATA_FOLDER,'pressure_profile.mat')
end
save (FILENAME,'P_initial_profile','plasma_beta_tot','P0');

if ADD_FIESTA_TO_FILENAMES==1
    FILENAME=strcat(DATA_FOLDER,'plasma_current_FIESTA.mat')
else
    FILENAME=strcat(DATA_FOLDER,'plasma_current.mat')
end
save (FILENAME,'Iaxis','JPHI_XZ_map','F2_profile')

if ADD_FIESTA_TO_FILENAMES==1
    FILENAME=strcat(DATA_FOLDER,'q_profile_FIESTA.mat')
else
    FILENAME=strcat(DATA_FOLDER,'q_profile.mat')
end
save (FILENAME,'q_initial_profile','q_XZ_map');

BpolX_initial_map=BR_XZ_map;
BpolZ_initial_map=BZ_XZ_map;
if ADD_FIESTA_TO_FILENAMES==1
    FILENAME=strcat(DATA_FOLDER,'B_fields_FIESTA.mat')
else
    FILENAME=strcat(DATA_FOLDER,'B_fields.mat')
end
save (FILENAME,'BpolX_initial_map','BpolZ_initial_map','Bphi0','Bphi_XZ_map');

if ADD_FIESTA_TO_FILENAMES==1
    FILENAME=strcat(DATA_FOLDER,'flux_geometry_FIESTA.mat')
else
    FILENAME=strcat(DATA_FOLDER,'flux_geometry.mat')
end
save (FILENAME,'psi_XZ_map','radial_XZ_map','psi_scale','tor_flux_scale','rho_tor_scale',...
    'theta_XZ_map','psi_global','volume_flux','surf_flux','radial_r_value_flux');


