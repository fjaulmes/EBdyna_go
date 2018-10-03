% for phi_rank=1:17
% imagesc(squeeze(Efield_X_map_phi(phi_rank,:,:)))
% pause(0.1)
% end

for phi_rank=1:33
    imagesc(squeeze(Epot_map_phi(phi_rank,:,:)))
    pause(0.1)
end