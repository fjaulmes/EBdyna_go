


psi_star_PR_map_DPHI=psi_star_PR_map_rank_next-psi_star_PR_map_rank_prev;
grad_psi_star_tor_PR_map=(psi_star_PR_map_DPHI./Rpos_PR_map(:,1:size_r))/(2*DOMEGA);

% converting the PR Efield on phi direction to an XZ map

% psi_data=reshape(grad_psi_star_tor_PR_map(:,1:size_r),NP*size_r,1);
% 
% grad_psi_star_tor_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,psi_data(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
% grad_psi_star_tor_XZ_map=grad_psi_star_tor_XZ_map';
% grad_psi_star_tor_XZ_map(isnan(grad_psi_star_tor_XZ_map))=0;
