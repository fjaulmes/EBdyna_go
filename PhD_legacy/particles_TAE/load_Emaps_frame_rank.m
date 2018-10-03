
f=frame_rank;

if (f<10)
    frame_name='00';
elseif (f<100)
    frame_name='0';
else
    frame_name='';
end

filename=strcat(EMAPS_FOLDER,'E0');
frame_name=strcat(frame_name,num2str(f));
filename=strcat(filename,frame_name,'.mat');
load(filename);

Efield_X_map_phi=TAE_amplitude*Efield_X_map_phi;
Efield_Z_map_phi=TAE_amplitude*Efield_Z_map_phi;
% Epot_map_phi=TAE_amplitude*Epot_map_phi;

% psi_map_phi=-(kTAE/omega_TAE)*Epot_map_phi;

% be careful with the sign here
% the stored derivative may have a wrong sign

% iEpot_map_phi=TAE_amplitude*iEpot_map_phi;
% grad_psi_star_map_phi=-(kTAE/omega_TAE)*iEpot_map_phi;