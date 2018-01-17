% for phi_rank=1:17
% imagesc(squeeze(Efield_X_map_phi(phi_rank,:,:)))
% pause(0.1)
% end

% for phi_rank=1:17
%     imagesc(squeeze(bsX_map_phi(phi_rank,:,:).^2+bsZ_map_phi(phi_rank,:,:).^2))
%     pause(0.1)
% end

% for phi_rank=1:17
%     imagesc(squeeze(psi_map_phi(phi_rank,:,:)))
%     pause(0.1)
% end

% for phi_rank=1:17
%     imagesc(squeeze(psi_star_map_phi(phi_rank,:,:)))
%     pause(0.1)
% end

for phi_rank=1:33
    imagesc(squeeze(Epot_map_phi(phi_rank,:,:)))
    pause(0.05)
end

%%
for phi_rank=1:33
    imagesc(squeeze(iEpot_map_phi(phi_rank,:,:)))
    pause(0.05)
end