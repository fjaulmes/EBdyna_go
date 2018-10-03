
% Bpol_XZ_map=sqrt(BpolX_XZ_map.^2+BpolZ_XZ_map.^2);

        
        Fmirror_X_map=gradB_X.*BpolX_XZ_map;
        Fmirror_Z_map=gradB_Z.*BpolZ_XZ_map;
        %Fmirror_phi_map=gradB_phi_XZ_map.*Bphi_XZsmall_map;
        
%         Fmirror_X_map=Fmirror_X_map./(Btot_XZ_map);
%         Fmirror_Z_map=Fmirror_Z_map./(Btot_XZ_map);
        %Fmirror_phi_map=Fmirror_phi_map./(Btot_XZ_map);

        Fmirror_tot_map=(Fmirror_X_map+Fmirror_Z_map)./(Btot_XZ_map);
        
%         Fmirror_tot_map=smooth_small_map(Fmirror_tot_map);
        
        run('calculate_rotB_corr');
        
        rotB_corr_map=rotB_X.*BpolX_XZ_map+rotB_Z.*BpolZ_XZ_map+rotB_phi.*Bphi_XZsmall_map;
        % mapping the small correction b x (rot(b)) 
        rotB_corr_map=rotB_corr_map./(Btot_XZ_map.^2);
        rotB_corr_map=smooth_small_map(rotB_corr_map);
       
        

        if MAPS_IN_PR_COORDINATES==1
            
            Fmirror_tot_PR_map=(gradB_X_PR_map.*BpolX_PR_map+gradB_Z_PR_map.*BpolZ_PR_map)./(Btot_PR_map);
            
%             v_data=reshape(Fmirror_tot_map(:,:)',sizeX*sizeZ,1);
%             Fmirror_tot_PR_map=dtinterp(XZ_mesh,XZ_mesh_dtri,v_data,X_PR_map(:,1:size_r),Z_PR_map(:,1:size_r),PR_DT_INTERPOLATION_METHOD);
% 
%             v_data=reshape(rotB_corr_map(:,:)',sizeX*sizeZ,1);
%             rotB_corr_PR_map=dtinterp(XZ_mesh,XZ_mesh_dtri,v_data,RR,Z_PR_map(:,1:size_r),PR_DT_INTERPOLATION_METHOD);
% 
%             v_data=reshape(Fmirror_X_map(:,:)',sizeX*sizeZ,1);
%             Fmirror_X_PR_map=dtinterp(XZ_mesh,XZ_mesh_dtri,v_data,RR,Z_PR_map(:,1:size_r),PR_DT_INTERPOLATION_METHOD);
%             
%             v_data=reshape(Fmirror_Z_map(:,:)',sizeX*sizeZ,1);
%             Fmirror_Z_PR_map=dtinterp(XZ_mesh,XZ_mesh_dtri,v_data,RR,Z_PR_map(:,1:size_r),PR_DT_INTERPOLATION_METHOD);
            
%             v_data=reshape(Fmirror_phi_map(:,:)',sizeX*sizeZ,1);
%             Fmirror_phi_PR_map=dtinterp(XZ_mesh,XZ_mesh_dtri,v_data,RR,Z_PR_map(:,1:size_r),PR_DT_INTERPOLATION_METHOD);
        end
        