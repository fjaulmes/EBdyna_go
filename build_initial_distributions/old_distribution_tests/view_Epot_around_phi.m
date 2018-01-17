
E_potential_PR_map0(:,:)=E_potential_PR_map_phi(1,:,:);
E_potential_PR_map(:,:)=zeros(NP,size_r);
    
for (phi_rank=2:NB_PHI_DATA_HALF)
    E_potential_PR_map(:,:)=zeros(NP,size_r);
    E_potential_PR_map(:,:)=E_potential_PR_map_phi(phi_rank,:,:);
    OFFSET=NP-(phi_rank-1)*4;
    E_potential_PR_map_recalc(1:NP-OFFSET,:)=E_potential_PR_map0(OFFSET:NP-1,:);
    E_potential_PR_map_recalc(NP-OFFSET+1:NP,:)=E_potential_PR_map0(1:OFFSET,:);
    imagesc(abs(E_potential_PR_map-E_potential_PR_map_recalc)');
    %imagesc((E_potential_PR_map)');
    axis xy;
    title(num2str(phi_rank));
    colorbar;
    pause(0.5);
end