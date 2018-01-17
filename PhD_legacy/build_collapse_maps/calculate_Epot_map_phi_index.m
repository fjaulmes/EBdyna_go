
        
        

        
        if PHI_OMEGA_RATIO==1
            phi_rank=phi_index
        else
            phi_rank=PHI_OMEGA_RATIO*phi_index-1
        end
        Epot_omega_map(:,:)=Epot_evol(frame_rank,:,:);
        
        E_potential_PR_map_phi=zeros(NB_THETA,NB_THETA,size_r);
        for index=1:NB_THETA
            %         phi_rank=round(2*phi_index-1);
            E_potential_PR_map_phi(index,:,:)=rotate_map_phi(Epot_omega_map,index)';
        end
                
        E_potential_PR_map(:,:)=E_potential_PR_map_phi(phi_rank,:,:);
        if phi_rank==1
            E_potential_PR_map_prev(:,:)=E_potential_PR_map_phi(end-1,:,:);
            E_potential_PR_map_next(:,:)=E_potential_PR_map_phi(phi_rank+1,:,:);
        elseif phi_rank==257
            E_potential_PR_map_prev(:,:)=E_potential_PR_map_phi(phi_rank-1,:,:);
            E_potential_PR_map_next(:,:)=E_potential_PR_map_phi(2,:,:);
        else
            E_potential_PR_map_prev(:,:)=E_potential_PR_map_phi(phi_rank-1,:,:);
            E_potential_PR_map_next(:,:)=E_potential_PR_map_phi(phi_rank+1,:,:);
        end
        if phi_rank==1
            E_potential_PR_map_prev_prev(:,:)=E_potential_PR_map_phi(end-2,:,:);
            E_potential_PR_map_next_next(:,:)=E_potential_PR_map_phi(phi_rank+2,:,:);
        elseif phi_rank==2
            E_potential_PR_map_prev_prev(:,:)=E_potential_PR_map_phi(end-1,:,:);
            E_potential_PR_map_next_next(:,:)=E_potential_PR_map_phi(phi_rank+2,:,:);
        elseif phi_rank==256
            E_potential_PR_map_prev_prev(:,:)=E_potential_PR_map_phi(phi_rank-2,:,:);
            E_potential_PR_map_next_next(:,:)=E_potential_PR_map_phi(2,:,:);
        elseif phi_rank==257
            E_potential_PR_map_prev_prev(:,:)=E_potential_PR_map_phi(phi_rank-2,:,:);
            E_potential_PR_map_next_next(:,:)=E_potential_PR_map_phi(3,:,:);
        else
            E_potential_PR_map_prev_prev(:,:)=E_potential_PR_map_phi(phi_rank-2,:,:);
            E_potential_PR_map_next_next(:,:)=E_potential_PR_map_phi(phi_rank+2,:,:);
        end        
        %     calculate_Efield_phi_rank;
        Epot_data=reshape(E_potential_PR_map(:,1:size_r),NP*size_r,1);
        
        Epot_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,Epot_data(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
        Epot_XZ_map=Epot_XZ_map';
        Epot_XZ_map(isnan(Epot_XZ_map))=0;
        
        calculate_lap_Phi;
        
        
        EX_map(:,:)=Efield_X_map_phi(phi_index,:,:);
        EZ_map(:,:)=Efield_Z_map_phi(phi_index,:,:);
        
        
        E_data=reshape(lap_Phi_PR_map(:,1:size_r),NP*size_r,1);
        lap_Phi_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,E_data(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
        lap_Phi_XZ_map=lap_Phi_XZ_map';
 
