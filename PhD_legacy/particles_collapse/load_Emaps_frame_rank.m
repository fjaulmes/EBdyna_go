
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
