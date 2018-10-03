psi_star_dot_omega_map(:,:)=psi_star_dot_evol(FRAME_NUMBER,:,:);

Efield_phi_map_phi=Efield_X_map_phi*0;

for phi_rank=1:NB_PHI
    
    phi_rank_kinetic_energy=0;
    EX_PR_map=squeeze(Efield_X_map_phi(phi_rank,:,:));
    EZ_PR_map=squeeze(Efield_Z_map_phi(phi_rank,:,:));
    
    psi_star_dot_PR_map_phi=zeros(NB_THETA,NB_THETA,size_r);
    for index=1:NB_THETA
        psi_star_dot_PR_map_phi(index,:,:)=rotate_map_phi(psi_star_dot_omega_map,index)';
    end
    
    grad_Phi_tor_PR_map=squeeze(grad_Phi_tor_map_phi(phi_rank,:,:));
    psi_star_dot_PR_map=squeeze(psi_star_dot_PR_map_phi(phi_rank,:,:));
    
    Ephi_PR_map=psi_star_dot_PR_map./Rpos_PR_map(:,1:size_r)-grad_Phi_tor_PR_map;
    
    % Ephi_PR_map=squeeze(Efield_phi_map_phi(phi_rank,:,:));
    Etot_PR_map=sqrt(EX_PR_map.^2+EZ_PR_map.^2+Ephi_PR_map.^2);
    
    Efield_phi_map_phi(phi_rank,:,:)=Ephi_PR_map;
    
end

save (strcat('dvD_map_frame_rank',num2str(FRAME_NUMBER),'.mat'),'-append','Efield_phi_map_phi');