save('XZsmall_fields_tokamak.mat','Btot_XZ_map','Bphi_XZsmall_map','BpolX_initial_XZsmall_map','BpolZ_initial_XZsmall_map','BpolX_final_XZsmall_map','BpolZ_final_XZsmall_map','Rpos_XZsmall_map','radial_XZsmall_map','psi_global');

Raxis=R0+X_axis;
save('map_dimensions.mat','mid_Xaxis_large','mid_Xzero','DPHI','DX','NP','NB_PHI','NB_PHI_DATA_HALF','R0','a','XX_small','ZZ_small','RR','Z_PR_map','scale_X','scale_Z','size_r','X_axis','Z_axis','mid_X','mid_Z','finesse_mesh','finesse_mesh_dtri','IDTRI','XZ_mesh','XZ_mesh_dtri','Raxis','radial_r_value_flux');

ZHe=2;
save('physics_constants.mat','kb','mD','mDT','mHe','mT','me','mu0','epsilon0','ZHe','eV');

Psih_initial=Psih;
save('reconnection_description.mat','tau_cr','size_r','Nomega','time_scale_precise','ksi0_evol_interp','rx_evol_interp','NB_TIME_STEP','TIME_STEP','Psih_initial','Psih_final','xPsih_zero','xPsih_max');

%pre-collapse-maps
save pre_collapse_XZsmall_maps.mat bX_XZ_map bZ_XZ_map bphi_XZ_map vD_X_XZ_map vD_Z_XZ_map vD_phi_XZ_map Fmirror_XZ_map
