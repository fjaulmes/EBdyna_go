
close all

REINIT_ALL_TOKAMAK_DATA=0;

if REINIT_ALL_TOKAMAK_DATA==1
    clear all
    initialize_folder_names;
    initialize_collapse_map_calculation_context
    rescaling_to_XZsmall_maps
    DT_INTERPOLATION_METHOD='quadratic'     % by default
    CALCULATE_VD_DATA_FILE=0;
end
% Epot_evol=-Epot_evol;
phi_index=12;

for f=1:101
    
    frame_rank=f
    
    psi_star_omega_map_half(:,:)=psi_star_2D_evol_lin(frame_rank,:,:);
    % Using symmetry to reconstruct a poloidal turn
    psi_star_omega_map_rank=zeros(size_r,NB_THETA);
    psi_star_omega_map_rank(:,1:round(0.5*NB_THETA))=psi_star_omega_map_half(:,:);
    psi_star_omega_map_rank(:,round(0.5*NB_THETA):NB_THETA)=psi_star_omega_map_half(:,round(0.5*NB_THETA):-1:1);
    
    psi_star_omega_map=zeros(size_r,NB_THETA);
    psi_star_PR_map=psi_star_omega_map';
    
    if PHI_OMEGA_RATIO==1
        phi_rank=phi_index;
    else
        phi_rank=PHI_OMEGA_RATIO*phi_index-1;
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
    E_potential_PR_map(:,:)=E_potential_PR_map_phi(phi_rank,:,:);
    
    if phi_rank>1
        E_potential_PR_map_prev(:,:)=E_potential_PR_map_phi(phi_rank-1,:,:);
    else
        E_potential_PR_map_prev(:,:)=E_potential_PR_map_phi(end-1,:,:);
    end
    E_potential_PR_map_next(:,:)=E_potential_PR_map_phi(phi_rank+1,:,:);
    
    figure(1);
    imagesc(psi_star_PR_map);
    
    figure(2);
    imagesc(E_potential_PR_map);
    colorbar
    
    pause
    
end
