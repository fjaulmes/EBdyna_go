mcc -m ./G_evolution_eq_RMP.m ./GT_many_evolution_eq_pll_RMP.m ../initialize_folder_names.m ./initialize_eq_sim_parameters_RMP.m ./initialize_eq_sim_maps.m ./initialize_eq_simulation_arrays.m ./interpolate_theta_psi_fromXZ.m ./time_step_integration_GT_eq.m ./time_step_integration_GT_eq_centrifug.m  ./adapt_speed_pphi_G.m ./adapt_speed_Ekin_G.m  ./update_GT_3D_eq.m ./update_GT_3D_eq_centrifug.m ./reposition_lost_particles_3DG.m ./wrap2pi.m ./interp2_XZ.m ./lininterp3.m ./build_3Dinterp_indexarrays.m ./record_time_stamp.m ./save_data_file_RMP.m ./reset_data_analysis_environment.m
cp ./G_evolution_eq_RMP ./G_evolution_eq
rm ./run_G_evolution_eq_RMP.sh