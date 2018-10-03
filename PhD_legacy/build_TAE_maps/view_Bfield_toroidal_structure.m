% for phi_rank=1:17
% imagesc(squeeze(Efield_X_map_phi(phi_rank,:,:)))
% pause(0.1)
% end

% for phi_rank=1:17
%     imagesc(squeeze(bsX_map_phi(phi_rank,:,:).^2+bsZ_map_phi(phi_rank,:,:).^2))
%     pause(0.1)
% end

for phi_rank=1:17
    imagesc(squeeze(BsX_map_phi(phi_rank,:,:)))
    pause(0.1)
end