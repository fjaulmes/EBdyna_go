% for phi_rank=1:17
% imagesc(squeeze(Efield_X_map_phi(phi_rank,:,:)))
% pause(0.1)
% end
phi_rank=1

for phi_rank=1:33
    imagesc(squeeze(Epot_map_phi(phi_rank,:,:)))
    pause(0.05)
end

for phi_rank=2:33
    imagesc(squeeze(Epot_map_phi(phi_rank,:,:)))
    pause(0.05)
end

% for phi_rank=1:33
%     imagesc(squeeze(-(kTAE/omega_TAE)*grad_psi_star_map_phi(phi_rank,:,:)))
%     pause(0.05)
% end