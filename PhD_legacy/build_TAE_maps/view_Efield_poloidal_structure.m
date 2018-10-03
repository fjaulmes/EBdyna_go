% for phi_rank=1:17
% imagesc(squeeze(Efield_X_map_phi(phi_rank,:,:)))
% pause(0.1)
% end

% for phi_rank=1:17
%     imagesc(squeeze(bsX_map_phi(phi_rank,:,:).^2+bsZ_map_phi(phi_rank,:,:).^2))
%     pause(0.1)
% end

for f=1:200
    if (f<10)
        frame_name='00';
    elseif (f<100)
        frame_name='0';
    else
        frame_name='';
    end
    frame_name=strcat(frame_name,num2str(f));
    filename='../E_maps/E0';
    filename=strcat(filename,frame_name,'.mat');   
    load(filename,'Epot_map_phi')

    imagesc(squeeze(Epot_map_phi(1,:,:)))
    title(frame_name)
    pause(0.05)
end