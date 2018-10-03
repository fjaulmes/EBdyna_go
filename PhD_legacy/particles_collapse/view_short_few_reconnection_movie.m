close all
figure(1);

for frame_movie=1:10:size(psipos_outputG,2)
    clf(1)
    grid on;
    hold on
    plot(1:size_r,1:size_r,'k--','linewidth',2)
    plot(psipos_outputG(:,1),psipos_outputG(:,frame_movie),'b.')
    pause
end