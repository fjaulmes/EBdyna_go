% for phi_rank=1:17
% imagesc(squeeze(Efield_X_map_phi(phi_rank,:,:)))
% pause(0.1)
% end

% for phi_rank=1:17
%     imagesc(squeeze(bsX_map_phi(phi_rank,:,:).^2+bsZ_map_phi(phi_rank,:,:).^2))
%     pause(0.1)
% end

for f=180:200
    if (f<10)
        frame_name='00';
    elseif (f<100)
        frame_name='0';
    else
        frame_name='';
    end
    frame_name=strcat(frame_name,num2str(f));
    filename='../B_maps/B0';
    filename=strcat(filename,frame_name,'.mat');   
    load(filename,'psi_map_phi')

    imagesc(squeeze(psi_map_phi(1,:,:)))
    title(frame_name)
    pause(0.05)
end