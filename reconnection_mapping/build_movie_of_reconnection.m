close all

writerObj = VideoWriter('reconnection_movie.avi');
writerObj.FrameRate=6
open(writerObj);

fig=figure;
set(fig, 'position', [150 150 w h]);

for (f=1:101)

    time_step=(f-1)
    
    
    if (time_step<100)
        frame_name='0';
    else
        frame_name='';
    end
    if (time_step<10)
        frame_name=strcat(frame_name,'0');
    end
    frame_name=strcat(frame_name,num2str(time_step));
    filename='cartoon\t';
    filename=strcat(filename,frame_name,'.jpg');
    disp(filename);
    
    A = imread(filename, 'JPG');
    image(A);
    axis off

    F = getframe(fig);
    writeVideo(writerObj,F);

end
close(fig);
% aviobj = close(aviobj);
close(writerObj);