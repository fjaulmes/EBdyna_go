

%         psi_star_dot_data=reshape(psi_star_dot_PR_map(:,1:size_r),NP*size_r,1);
%         psi_star_dot_XZ_zoom_map=dtinterp(finesse_mesh,finesse_mesh_dtri,psi_star_dot_data(IDTRI),XX_zoom,ZZ_zoom,DT_INTERPOLATION_METHOD);
%         psi_star_dot_XZ_zoom_map=psi_star_dot_XZ_zoom_map';
%         psi_star_dot_XZ_zoom_map(isnan(psi_star_dot_XZ_zoom_map))=0;
%         psi_star_dot_XZ_zoom_map=smooth_small_map(psi_star_dot_XZ_zoom_map);
        
%         Epot_data=reshape(E_potential_PR_map(:,1:size_r),NP*size_r,1);
%         Epot_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,Epot_data(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
%         Epot_XZ_map=Epot_XZ_map';
%         Epot_XZ_map(isnan(Epot_XZ_map))=0;
		
Epot_data=reshape(E_potential_PR_map(:,1:size_r),NP*size_r,1);
Epot_XZ_map=gridfit(finesse_data_X,finesse_data_Z,Epot_data,scale_X,scale_Z,'smoothness',INTERP_SMOOTHNESS);
Epot_XZ_map=Epot_XZ_map'; 
 
        calculate_Efield_XZ;
%         EX_XZ_map=smooth_small_map(EX_XZ_map);
%         EZ_XZ_map=smooth_small_map(EZ_XZ_map);
        
        E_potential_PR_map_DPHI=E_potential_PR_map_next-E_potential_PR_map_prev;
        grad_Phi_tor_PR_map=(E_potential_PR_map_DPHI./Rpos_PR_map(:,1:size_r))/(2*DOMEGA);
        
%         %use saved data
%          grad_Phi_tor_PR_map=squeeze(grad_Phi_tor_map_phi(phi_index,:,:));

%         grad_Phi_tor_PR_map=smooth_small_map(grad_phi_tor_num);
%         grad_Phi_tor_PR_map=grad_phi_tor_num;
        Ephi_PR_map=psi_star_dot_PR_map./Rpos_PR_map(:,1:size_r)-grad_Phi_tor_PR_map;

% Ephi_data=reshape(Ephi_PR_map(:,1:size_r),NP*size_r,1);
% Ephi_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,Ephi_data(IDTRI),XX_small,ZZ_small,PR_DT_INTERPOLATION_METHOD);
% Ephi_XZ_map=Ephi_XZ_map';

%         E_data=reshape(EX_XZ_map(:,:)',sizeX*sizeZ,1);
%         EX_PR_map=dtinterp(XZ_mesh,XZ_mesh_dtri,E_data,RR,Z_PR_map(:,1:size_r),PR_DT_INTERPOLATION_METHOD);
%         E_data=reshape(EZ_XZ_map(:,:)',sizeX*sizeZ,1);
%         EZ_PR_map=dtinterp(XZ_mesh,XZ_mesh_dtri,E_data,RR,Z_PR_map(:,1:size_r),PR_DT_INTERPOLATION_METHOD);
% 
        
        % converting the PR Efield on phi direction to an XZ map

        %if MAPS_IN_PR_COORDINATES==0
        %end
        
%         Ephi_XZ_map=smooth_small_map(Ephi_XZ_map);