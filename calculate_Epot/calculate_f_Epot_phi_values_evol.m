    %f=round(0.1*(reference_frame-1));
    f=frame_rank;
	f0=f;
    
    rx_value=interp1(time_evol,rx_evol,time_scale_lin(f))
    ksi0_value=interp1(time_evol,ksi0_evol,time_scale_lin(f));
    rx_precise=interp1(radial_r_value_flux,1:Nradial, rx_value);
    rx_rank=round(rx_precise);
    minX=max(mid_X_zoom-round(rx_value/DX_zoom)-X_axis_offset,3);
    maxX=min(mid_X_zoom+round(rx_value/DX_zoom)+X_axis_offset,NZ_zoom-2);
    minZ=max(mid_Z_zoom-round(elongation*rx_value/DZ_zoom)-Z_axis_offset,3);
    maxZ=min(mid_Z_zoom+round(elongation*rx_value/DZ_zoom)+Z_axis_offset,NZ_zoom-2);
    
    psi_limit13=max(interp1(1:Nradial,psi_star_initial,rx_precise),0)
    
    x_pos_sep_final=interp1(psi_star_final,1:Nradial,psi_limit13);
    

    delta_rx=0;
    define_RZ_mask_rx_rank;
    mask_XZ_zoom_map_default=mask_XZ_zoom_map_reconnection;
    
    if (f<10)
        frame_name='00';
    elseif (f<100)
        frame_name='0';
    else
        frame_name='';
    end
    frame_name=strcat(frame_name,num2str(f));
    filename=strcat('./',RECONNECTION_MAPS_FOLDER);
    filename=strcat(filename,'t0');
    filename=strcat(filename,frame_name,'.mat');
           
    % Looping around in the toroidal direction
%    for (phi_rank=1:NB_PHI)
   for (phi_rank=1:PHI_STEP_SIZE:NB_PHI)
       avg_phi_rank=((phi_rank-1)/PHI_STEP_SIZE)+1
%     for (phi_rank=NB_PHI_DATA_HALF-1:NB_PHI_DATA_HALF+1)
        
%         disp('**************************************************************');
        
        clear psi_star_dot_XZ_zoom_map psi_star_XZ_zoom_map;
        
        if SAVE_LOG_FILE==1
            fprintf(fid, '****** phi_rank = %d ********************************\n',phi_rank);
            fprintf(fid, '********************************************************\n');
        end
        
% tic

        f0=f;
        phi=(1/(NB_PHI-1))*(phi_rank-1)*(2*pi);
        
        calculate_psi_star_dot_XZ_map_v2;
        f=frame_rank;
        f0=f;
        initialize_rotated_PR_maps;
        %disp('**************************************************************');
        
        psi_star_dot_data=reshape(psi_star_dot_PR_map(:,1:size_r),NP*size_r,1);
        psi_star_dot_XZ_zoom_map=griddata(finesse_data_X,finesse_data_Z,psi_star_dot_data,XX_zoom,ZZ_zoom,'cubic');
%         psi_star_dot_XZ_zoom_map=dtinterp(finesse_mesh,finesse_mesh_dtri,psi_star_dot_data(IDTRI),XX_zoom,ZZ_zoom,DT_INTERPOLATION_METHOD);
        psi_star_dot_XZ_zoom_map=psi_star_dot_XZ_zoom_map';
        psi_star_dot_XZ_zoom_map(isnan(psi_star_dot_XZ_zoom_map))=0;
%         psi_star_dot_XZ_zoom_map=smooth_small_map(psi_star_dot_XZ_zoom_map);
        
        psi_star_data=reshape(psi_star_PR_map(:,1:size_r),NP*size_r,1);
        psi_star_XZ_zoom_map=griddata(finesse_data_X,finesse_data_Z,psi_star_data,XX_zoom,ZZ_zoom,'cubic');
%         psi_star_XZ_zoom_map=dtinterp(finesse_mesh,finesse_mesh_dtri,psi_star_data(IDTRI),XX_zoom,ZZ_zoom,DT_INTERPOLATION_METHOD);
        psi_star_XZ_zoom_map=psi_star_XZ_zoom_map';
        psi_star_XZ_zoom_map(isnan(psi_star_XZ_zoom_map))=0;

        psi_star_XZ_zoom_map_copy=psi_star_XZ_zoom_map;
        
        %disp('**************************************************************');
        
        mask_XZ_zoom_map_reconnection=mask_XZ_zoom_map_default;
        
        close all;
        calculate_Bstar_RZ;
        %smoothing required here
        Bstar_XZ_zoom_map=smooth_small_map(Bstar_XZ_zoom_map);
        
        %disp('**************************************************************');
        
        %region1_rank=round((rx_value-ksi0_value)/Dr)+1;
        
        integ_element_XZ_map=zeros(NZ_zoom,NZ_zoom);
        integ_element_XZ_map=B3_XZ_zoom_map.*psi_star_dot_XZ_zoom_map./Bstar_XZ_zoom_map;
        
        if ADDITIONAL_MAP_SMOOTHING==1
            integ_element_XZ_map=smooth_small_map(integ_element_XZ_map);
        end

% toc
% tic

        
        %disp('**************************************************************');
        
        define_reference_potential_axis;
        
        calculate_XZ_psi_star_contours_halves;
        
        if DISPLAY_OUTPUTS==1
            display_integration_map_for_E_potential;
            pause(0.1);
        end
% toc
% tic

        mask_XZ_zoom_map_reconnection=mask_XZ_zoom_map_default;
        
        if INTEGRATION_PRECISION==0
            smooth_XZ_psi_star_contour_integration_v2;
        else
            smooth_XZ_psi_star_contour_integration;
        end
        
        
        %run('display_E_potential_RZ')
        close all
        
        if ADDITIONAL_MAP_SMOOTHING==1
            E_potential_RZ_map=smooth_small_map(E_potential_RZ_map);
        end
        
        % working on smaller maps for performance issues
        %NZ=4*NX;
        smooth_potential_data=reshape(E_potential_RZ_map(:,:)',NZ*NZ,1);
        
        %Recording the potential map and Bstar fields
        %E_potential_PR_map=griddata(X_scale_data,Z_scale_data,smooth_potential_data,RR,Z_PR_map(:,1:size_r),'cubic');
        E_potential_PR_map=dtinterp(XZ_mesh,XZ_mesh_dtri,smooth_potential_data,RR,Z_PR_map(:,1:size_r),PR_DT_INTERPOLATION_METHOD);
        
        %if DISPLAY_OUTPUTS==1
            imagesc(X_scale(NX:3*NX),Z_scale(NX:3*NX),E_potential_RZ_map(NX:3*NX,NX:3*NX)');
            axis xy;
            colorbar;
            pause(0.1);
        %end
        disp('****************E_potential_PR_map DONE!**********************');
        
%         Bstar_data=reshape(Bstar_X_RZ_map(:,:)',NZ*NZ,1);
% %                  BstarX_PR_map=griddata(X_scale_data,Z_scale_data,Bstar_data,RR,Z_PR_map(:,1:size_r),'cubic');
%         BstarX_PR_map=dtinterp(XZ_mesh,XZ_mesh_dtri,Bstar_data,RR,Z_PR_map(:,1:size_r),PR_DT_INTERPOLATION_METHOD);
%         
%         Bstar_data=reshape(Bstar_Z_RZ_map(:,:)',NZ*NZ,1);
% %                  BstarZ_PR_map=griddata(X_scale_data,Z_scale_data,Bstar_data,RR,Z_PR_map(:,1:size_r),'cubic');
%         BstarZ_PR_map=dtinterp(XZ_mesh,XZ_mesh_dtri,Bstar_data,RR,Z_PR_map(:,1:size_r),PR_DT_INTERPOLATION_METHOD);
        
        %     E_potential_PR_map=E_potential_PR_map';
        %     BstarX_PR_map=BstarX_PR_map';
        %     BstarZ_PR_map=BstarZ_PR_map';
        %     psi_star_dot_PR_map=psi_star_dot_PR_map';
        
        E_potential_PR_map_phi(avg_phi_rank,:,:)=E_potential_PR_map;
%         BstarX_PR_map_phi(phi_rank,:,:)=BstarX_PR_map;
%         BstarZ_PR_map_phi(phi_rank,:,:)=BstarZ_PR_map;
        psi_star_dot_PR_map_phi(avg_phi_rank,:,1:size_r)=psi_star_dot_PR_map(:,1:size_r);
        %psi_star_PR_map_phi(phi_rank,:,1:size_r)=psi_star_PR_map(:,1:size_r);
        
        
        % E_potential_PR_data=reshape(E_potential_PR_map(1:Nradial,:),NP*Nradial,1);
        % E_potential_XZ_zoom_map=griddata(finesse_data_X,finesse_data_Z,E_potential_PR_data,XX_zoom,ZZ_zoom,'cubic');
        
        if SAVE_LOG_FILE==1
            fprintf(fid, '********************************************************\n');
        end

% toc

   end
    
   if SIGN_PSI0==-1
    %E_potential_PR_map_phi=-E_potential_PR_map_phi;
    psi_star_dot_PR_map_phi=-psi_star_dot_PR_map_phi;
   end
   
    if (SAVE_DATA_FILE==1)
%        save(filename,'E_potential_PR_map_phi','psi_star_PR_map_phi','psi_star_dot_PR_map_phi');
        save(filename,'E_potential_PR_map_phi','psi_star_dot_PR_map_phi');
        disp('****** SAVING FILE in \reconnection_maps\ ****');
        disp(filename);
    end
    disp('**************************************************************');