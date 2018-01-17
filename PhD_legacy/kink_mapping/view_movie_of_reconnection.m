f=1

xMixing_radius=size_r-3;

initialize_map_r_omega;

for (f=1:101)
    fig=figure(1)
    clear psi_star_2D
    clear F gcf
    psi_star_2D=zeros(Nomega,size_r)';
    psi_star_2D(:,:)=psi_star_2D_evol(f,:,:);
    % psi_star_2D(:,:)=psi_star_2D_evol(f+1,:,:)-psi_star_2D_evol(f-1,:,:);
    
    %     psi_star_2D=psi_star_2D';
    
    % Create grids and convert polar coordinates to rectangularfigure(3)
    subplot(2,1,1)
    [THETA,RR] = meshgrid(scale_omega,scale_r);
    [XX,YY] = pol2cart(THETA,RR);
    
    
    %figure(3)
    surfc(XX,YY,psi_star_2D,'edgecolor','none');view(0,90);
    xlim([-0.5 0.5])
    ylim([0 0.5])
    %     contour(XX,YY,psi_star_2D,psi_contour_values);
    
    max_psi_star=max(max(psi_star_2D));
    title(f)
    
    
    sp=subplot(2,1,2);
    cla(sp);
    grid on;hold on;
    psi_star_2D(:,:)=psi_star_2D_evol(f,:,:);
    plot([-radial_r_value_flux(size_r:-1:1) radial_r_value_flux(2:size_r)],[psi_star_2D(end:-1:1,end)' psi_star_2D(2:end,1)' ]);
    xlim([-0.5 0.5])
    time_step=(f-1)
    
    figure(1)
    pause(0.4);
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
    title(strcat('frame rank #',num2str(frame_name)));

    
    F=getframe(fig);
[h, w, p] = size(F.cdata);  % use 1st frame to get dimensions
% hf = figure; 
% % resize figure based on frame's w x h, and place at (150, 150)
% set(hf, 'position', [150 150 w h]);
% axis off
    [im,map] = frame2im(F);    %Return associated image data
    imwrite(im,filename,'jpg');
    saveas(gcf,filename,'jpg');%saving freq plot
    
    
    pause(0.4);

end
