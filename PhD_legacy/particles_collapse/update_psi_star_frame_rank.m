

psi_star_omega_map_half(:,:)=psi_star_2D_evol_lin(round(frame_rank_next/10)+1,:,:);


% Using symmetry to reconstruct a poloidal turn
% psi_star_omega_map_rank=zeros(size_r,NB_THETA);

psi_star_omega_map_rank(:,1:round(0.5*NB_THETA))=psi_star_omega_map_half(:,:);
psi_star_omega_map_rank(:,round(0.5*NB_THETA):NB_THETA)=psi_star_omega_map_half(:,round(0.5*NB_THETA):-1:1);


