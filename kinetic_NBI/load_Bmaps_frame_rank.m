
% f=frame_rank;

if (frame_rank<10)
    frame_name='00';
elseif (frame_rank<100)
    frame_name='0';
else
    frame_name='';
end

filename=strcat(BMAPS_FOLDER,'B0');
frame_name=strcat(frame_name,num2str(frame_rank));
filename=strcat(filename,frame_name,'.mat')
load(filename);



