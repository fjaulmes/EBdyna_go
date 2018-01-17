
        
        
        psi_star_omega_map_half(:,:)=psi_star_2D_evol_lin(frame_rank,:,:);
        % Using symmetry to reconstruct a poloidal turn
        psi_star_omega_map_rank=zeros(size_r,NB_THETA);
        psi_star_omega_map_rank(:,1:round(0.5*NB_THETA))=psi_star_omega_map_half(:,:);
        psi_star_omega_map_rank(:,round(0.5*NB_THETA):NB_THETA)=psi_star_omega_map_half(:,round(0.5*NB_THETA):-1:1);
        
        psi_star_omega_map=zeros(size_r,NB_THETA);
        psi_star_PR_map=psi_star_omega_map';
        
        if PHI_OMEGA_RATIO==1
            phi_rank=phi_index
        else
            phi_rank=PHI_OMEGA_RATIO*phi_index-1
        end
        
        % if phi_index==1
        %     phi_rank=1
        % end
        
        psi_star_omega_map(:,:)=rotate_map_phi(psi_star_omega_map_rank,phi_rank);
        psi_star_PR_map=psi_star_omega_map';
        
        delta_phi_rank=phi_rank+1;
        psi_star_omega_map_rank_next(:,:)=rotate_map_phi(psi_star_omega_map_rank,delta_phi_rank);
        psi_star_PR_map_rank_next=psi_star_omega_map_rank_next';
        
        delta_phi_rank=phi_rank-1;
        psi_star_omega_map_rank_prev(:,:)=rotate_map_phi(psi_star_omega_map_rank,delta_phi_rank);
        psi_star_PR_map_rank_prev=psi_star_omega_map_rank_prev';
        
        calculate_Btot_phi_rank;
        
        
        Epot_omega_map(:,:)=Epot_evol(frame_rank,:,:);
        
        E_potential_PR_map_phi=zeros(NB_THETA,NB_THETA,size_r);
        for index=1:NB_THETA
            %         phi_rank=round(2*phi_index-1);
            E_potential_PR_map_phi(index,:,:)=rotate_map_phi(Epot_omega_map,index)';
        end
        
        psi_star_dot_omega_map(:,:)=psi_star_dot_evol(frame_rank,:,:);
        
        psi_star_dot_PR_map_phi=zeros(NB_THETA,NB_THETA,size_r);
        for index=1:NB_THETA
            %         phi_rank=round(2*phi_index-1);
            psi_star_dot_PR_map_phi(index,:,:)=rotate_map_phi(psi_star_dot_omega_map,index)';
        end
        
        psi_star_dot_PR_map(:,:)=psi_star_dot_PR_map_phi(phi_rank,:,:);
        % load psi_star_dot_PR_map_ref3.mat
        % psi_star_dot_PR_map(:,:)=psi_star_dot_PR_map_ref3(:,1:size_r);
        
        %     E_data=reshape(psi_star_dot_PR_map(:,1:size_r),NP*size_r,1);
        %     % Efield_3_XZ_map=griddata(finesse_data_X(1:size_r*NP),finesse_data_Z(1:size_r*NP),E_data,XX_small,ZZ_small,'cubic');
        %     Efield_3_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,E_data(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
        %     Efield_3_XZ_map=Efield_3_XZ_map';
        %     Efield_3_XZ_map(isnan(Efield_3_XZ_map))=0;
        
        %     Ephi_XZ_map=Efield_3_XZ_map./Rpos_XZsmall_map;
        % %
        % Bfield_3_XZmap=Bphi_XZsmall_map./Rpos_XZsmall_map;
        %
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
        
        %     calculate_Efield_phi_rank;
        Epot_data=reshape(E_potential_PR_map(:,1:size_r),NP*size_r,1);
        
        Epot_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,Epot_data(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
        Epot_XZ_map=Epot_XZ_map';
        Epot_XZ_map(isnan(Epot_XZ_map))=0;
        
        calculate_Efield_XZ;
        
        E_potential_PR_map_DPHI=E_potential_PR_map_next-E_potential_PR_map_prev;
        grad_Phi_tor_PR_map=(E_potential_PR_map_DPHI./Rpos_PR_map(:,1:size_r))/(2*DOMEGA);
        Ephi_PR_map=psi_star_dot_PR_map./Rpos_PR_map(:,1:size_r)-grad_Phi_tor_PR_map;
        
        EX_map(:,:)=Efield_X_map_phi(phi_index,:,:);
        EZ_map(:,:)=Efield_Z_map_phi(phi_index,:,:);
        
        
        %     vExB_X_map=EZ_XZ_map.*Bphi_XZsmall_map-Ephi_XZ_map.*BpolZ_XZ_map;
        %     vExB_Z_map=Ephi_XZ_map.*BpolX_XZ_map-Bphi_XZsmall_map.*EX_XZ_map;
        %     vExB_phi_map=EX_XZ_map.*BpolZ_XZ_map-BpolX_XZ_map.*EZ_XZ_map;
        %     vExB_X_map=vExB_X_map./(Btot_XZ_map.^2);
        %     vExB_Z_map=vExB_Z_map./(Btot_XZ_map.^2);
        %     vExB_phi_map=vExB_phi_map./(Btot_XZ_map.^2);
        
        
        % dr_data=reshape(dr_X_PR_map(:,1:size_r),NP*size_r,1);
        % drX_map=dtinterp(finesse_mesh,finesse_mesh_dtri,dr_data(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
        % drX_map=drX_map';
        % dr_data=reshape(dr_Z_PR_map(:,1:size_r),NP*size_r,1);
        % drZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,dr_data(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
        % drZ_map=drZ_map';
        %
        % vExB_radial=vExB_X_map.*drX_map+vExB_Z_map.*drZ_map;
        
        %     Btot_PR_map=sqrt(BpolX_PR_map.^2+BpolZ_PR_map.^2+Btor_PR_map(:,1:size_r).^2);
        
        vExB_X_PR_map=EZ_PR_map.*Btor_PR_map(:,1:size_r)-Ephi_PR_map.*BpolZ_PR_map;
        vExB_Z_PR_map=Ephi_PR_map.*BpolX_PR_map-Btor_PR_map(:,1:size_r).*EX_PR_map;
        vExB_phi_PR_map=EX_PR_map.*BpolZ_PR_map-BpolX_PR_map.*EZ_PR_map;
        vExB_X_PR_map=vExB_X_PR_map./(Btot_PR_map(:,1:size_r).^2);
        vExB_Z_PR_map=vExB_Z_PR_map./(Btot_PR_map(:,1:size_r).^2);
        vExB_phi_PR_map=vExB_phi_PR_map./(Btot_PR_map(:,1:size_r).^2);
        
        vExB_radial_PR_map=vExB_X_PR_map.*dr_X_PR_map(:,1:size_r)+vExB_Z_PR_map.*dr_Z_PR_map(:,1:size_r);
        vExB_poloidal_PR_map=vExB_X_PR_map.*grad_theta_PR_map_X(:,1:size_r)+vExB_Z_PR_map.*grad_theta_PR_map_Z(:,1:size_r);
        vExB_toroidal_PR_map=vExB_phi_PR_map;
        
        vE_data=reshape(vExB_radial_PR_map(:,1:size_r),NP*size_r,1);
        vEr_map=dtinterp(finesse_mesh,finesse_mesh_dtri,vE_data(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
        vEr_map=vEr_map';
        vEr_map(isnan(vEr_map))=0;
%         vEr_map(mask_XZ_small==0)=0;

        vE_data=reshape(vExB_X_PR_map(:,1:size_r),NP*size_r,1);
        vEX_map=dtinterp(finesse_mesh,finesse_mesh_dtri,vE_data(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
        vEX_map=vEX_map';
%         vEX_map(mask_XZ_small==0)=0;

        vE_data=reshape(vExB_Z_PR_map(:,1:size_r),NP*size_r,1);
        vEZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,vE_data(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
        vEZ_map=vEZ_map';
%         vEZ_map(mask_XZ_small==0)=0;

        
        E_data=reshape(EZ_PR_map(:,1:size_r),NP*size_r,1);
        EZ_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,E_data(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
        EZ_XZ_map=EZ_XZ_map';
        
        E_data=reshape(EX_PR_map(:,1:size_r),NP*size_r,1);
        EX_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,E_data(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
        EX_XZ_map=EX_XZ_map';

        E_data=reshape(Ephi_PR_map(:,1:size_r),NP*size_r,1);
        Ephi_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,E_data(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
        Ephi_XZ_map=Ephi_XZ_map';
 
        
        B_data=reshape(BpolZ_PR_map(:,1:size_r),NP*size_r,1);
        BZ_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,B_data(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
        BZ_XZ_map=BZ_XZ_map';
        
        B_data=reshape(BpolX_PR_map(:,1:size_r),NP*size_r,1);
        BX_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,B_data(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
        BX_XZ_map=BX_XZ_map';

        B_data=reshape(Btor_PR_map(:,1:size_r),NP*size_r,1);
        Bphi_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,B_data(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
        Bphi_XZ_map=Bphi_XZ_map';
               