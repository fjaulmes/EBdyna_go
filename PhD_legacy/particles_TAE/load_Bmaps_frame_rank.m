
f=frame_rank;

if (f<10)
    frame_name='00';
elseif (f<100)
    frame_name='0';
else
    frame_name='';
end

filename=strcat(BMAPS_FOLDER,'B0');
frame_name=strcat(frame_name,num2str(f));
filename=strcat(filename,frame_name,'.mat');
load(filename);

BsX_map_phi=TAE_amplitude*BsX_map_phi;
BsZ_map_phi=TAE_amplitude*BsZ_map_phi;
Bsphi_map_phi=TAE_amplitude*Bsphi_map_phi;
% psi_map_phi=TAE_amplitude*psi_map_phi;

