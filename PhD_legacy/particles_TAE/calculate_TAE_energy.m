bsX_PR_map=volume_circ_approx_PR_map(:,1:size_r)*0;
bsZ_PR_map=volume_circ_approx_PR_map(:,1:size_r)*0;
bs_PR_map=volume_circ_approx_PR_map(:,1:size_r)*0;
W_PR_map=volume_circ_approx_PR_map(:,1:size_r)*0;
W_TAE=0;

for phi_rank=1:NB_PHI-1
    bsX_PR_map=squeeze(bsX_map_phi(phi_rank,:,:));
    bsZ_PR_map=squeeze(bsZ_map_phi(phi_rank,:,:));
    bs_PR_map=bsX_PR_map.^2+bsZ_PR_map.^2;
    W_PR_map=0.5*bs_PR_map/mu0;
    W_TAE=W_TAE+sum(sum(W_PR_map.*volume_circ_approx_PR_map(:,pTAE_inf:pTAE_sup)));
end

W_TAE=nTAE*W_TAE;