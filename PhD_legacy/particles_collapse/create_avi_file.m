% aviobj = avifile('bulk_redistribution.avi','compression','none');

writerObj = VideoWriter('bulk_redistribution.avi');
writerObj.FrameRate=8
open(writerObj);
fig=figure;

    filename='cartoon\t0001';
    filename=strcat(filename,'.bmp');
    disp(filename);
    
    A = imread(filename, 'BMP');
    image(scale_X,scale_Z,A);
    F = getframe(fig);
    writeVideo(writerObj,F);


for k=1:400
    time_step=k*50;
    if (time_step<10000)
        frame_name='0';
    else
        frame_name='';
    end
    if (time_step<1000)
        frame_name=strcat(frame_name,'0');
    end
    if (time_step<100)
        frame_name=strcat(frame_name,'0');
    end
    frame_name=strcat(frame_name,num2str(time_step));
    filename='cartoon\t';
    filename=strcat(filename,frame_name,'.bmp');
    disp(filename);
    
    A = imread(filename, 'BMP');
    image(scale_X,scale_Z,A);
    F = getframe(fig);
    writeVideo(writerObj,F);

%     aviobj = addframe(aviobj,F);
end
close(fig);
% aviobj = close(aviobj);
close(writerObj);