close all
w=920
h=905

writerObj = VideoWriter('W_centrifug_fc1_AUG30809_2p25_Elim1000_movie_pi_div_8_mod_slice.avi');
writerObj.FrameRate=4.0
open(writerObj);

fig=figure;
set(fig, 'position', [150 150 w h]);


%%
for (f=1:99)
    close all
    fig=figure;
    set(fig, 'position', [150 150 w h]);

    time_step=(f)
    
    frame_name='t'
    if (time_step<10)
        frame_name=strcat(frame_name,'0');
    end
    frame_name=strcat(frame_name,num2str(time_step));
    filename='cartoon/';
    filename=strcat(filename,frame_name,'.png');
    disp(filename);
    
    A = imread(filename);
    image(A);
    axis off

    F = getframe(fig);
    writeVideo(writerObj,F);
    pause(0.1)

end

F = getframe(fig);
writeVideo(writerObj,F);
pause(0.1)

F = getframe(fig);
writeVideo(writerObj,F);
pause(0.1)

F = getframe(fig);
writeVideo(writerObj,F);
pause(0.1)

close(fig);
% aviobj = close(aviobj);
close(writerObj);