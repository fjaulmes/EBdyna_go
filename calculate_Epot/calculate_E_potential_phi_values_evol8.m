
close all;

REINIT_ALL_TOKAMAK_DATA=1;


if REINIT_ALL_TOKAMAK_DATA==1
    clear all;
    warning off MATLAB:griddata:DuplicateDataPoints;
    initialize_folder_names;
    initialize_Epot_calculation_context;
    
    finesse_mesh=[finesse_data_X  finesse_data_Z];
    [finesse_mesh, IDTRI, J] = unique(finesse_mesh,'last','rows');
    finesse_mesh_dtri=DelaunayTri(finesse_mesh);
    
    XZ_mesh=[X_scale_data Z_scale_data];
    XZ_mesh_dtri=DelaunayTri(XZ_mesh);
    
    clc
    
else
%     load Epot_calculation_context
end


%simulation options
initialize_Epot_calculation_parameters;


for frame_rank=43:48
    clear E_potential_PR_map_phi BstarX_PR_map_phi BstarZ_PR_map_phi psi_star_dot_PR_map_phi psi_star_PR_map_phi;
    
    disp('frame_rank = ');
    disp(frame_rank);
    
	% had been done in case we need more than 101 frames......
    reference_frame=min((frame_rank-1)*10+1,1001)
    calculate_f_Epot_phi_values_evol;

%     if (reference_frame>EXPULSION_FRAME_NUMBER_INF) && (reference_frame<EXPULSION_FRAME_NUMBER_SUP)
%         % intermediate frame for transition at core expulsion
%         reference_frame=reference_frame+5
%         calculate_f_Epot_phi_values_evol;
%     end
end

if SAVE_LOG_FILE==1
    fclose(fid);
    disp('****** log file closed properly ****');
end
if SAVE_DATA_FILE==1
    disp('****** data is available in \reconnection_maps\ ****');
else
    disp('****** size of data ****');
    size(E_potential_PR_map_phi)
end


exit
