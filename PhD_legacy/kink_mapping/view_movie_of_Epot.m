f=1
DATA_FOLDER='../data_tokamak/';

xMixing_radius=size_r-3;
filename=strcat(DATA_FOLDER,'physics_constants.mat');
load(filename);
filename=strcat(DATA_FOLDER,'psi_star_evol.mat');
load(filename);
filename=strcat(DATA_FOLDER,'Epot_psi_star_dot_evol.mat');
load(filename);
filename=strcat(DATA_FOLDER,'tokamak_map_dimensions.mat');
load(filename);
filename=strcat(DATA_FOLDER,'tokamak_PR_map.mat');
load(filename);
filename=strcat(DATA_FOLDER,'B_fields.mat');
load(filename);
filename=strcat(DATA_FOLDER,'flux_geometry.mat');
load(filename);
filename=strcat(DATA_FOLDER,'psi_profiles.mat');
load(filename);
filename=strcat(DATA_FOLDER,'psi_profiles_kadomtsev.mat');
load(filename);
filename=strcat(DATA_FOLDER,'grad_flux_geometry.mat');
load(filename);

% initialize_map_r_omega;
size_omega=pi;
Nomega=round((NP+1)/2);
Domega=size_omega/(Nomega-1);
half_Nomega=round(Nomega/2);
scale_omega=[0:Nomega-1]*Domega;
scale_r=radial_r_value_flux(1:size_r);

[THETA,RR] = meshgrid(scale_omega,scale_r);
[XX,YY] = pol2cart(THETA,RR);



[xxs yys ]=meshgrid((-0.5:0.001:0.5),(0:0.001:0.5));

for (f=2:100)
    time_step=f;
    fig=figure(1)
    clear psi_star_2D
    clear F gcf
    psi_star_2D=zeros(Nomega,size_r)';
    Epot_2D=zeros(Nomega,size_r)';
%     psi_star_2D(:,:)=psi_star_2D_evol(f,:,:);
%     psi_star_2D(:,:)=psi_star_2D_evol(f+1,:,:)-psi_star_2D_evol(f-1,:,:);
    psi_star_2D(:,:)=psi_star_2D_evol_lin(f,:,:);
    Epot_2D(:,:)=Epot_evol(f,1:size_r,1:Nomega);
    
    %     psi_star_2D=psi_star_2D';
    
    % Create grids and convert polar coordinates to rectangularfigure(3)
    sp=subplot(2,1,1)
    set(gca,'fontsize',20);
    
    %figure(3)
    psi_XY_map=griddata(XX,YY,psi_star_2D,xxs,yys,'linear');
    psi_XY_map=psi_XY_map';
    cla(sp);
    imagesc((-0.5:0.001:0.5),(0:0.001:0.5),psi_XY_map',[-0.6 0.001]*1e-2);
    axis xy;
%     surfc(XX,YY,psi_star_2D,'edgecolor','none');view(0,90);
    xlim([-0.3 0.3])
    ylim([0 0.3])
    %     contour(XX,YY,psi_star_2D,psi_contour_values);
    
    max_psi_star=max(max(psi_star_2D));
    title('\psi_*')
    
%     [lEpot wEpot]=size(Epot_2D);
%     DEpot=1/(round(0.5*lEpot)-1);
%     Cmapping=meshgrid([-1-DEpot:DEpot:1]*1e3,1:wEpot);
%     
    sp=subplot(2,1,2);
    set(gca,'fontsize',20);
    if (time_step<100)
        frame_name='0';
    else
        frame_name='';
    end
    if (time_step<10)
        frame_name=strcat(frame_name,'0');
    end
    frame_name=strcat(frame_name,num2str(time_step));
    filename='cartoon/t';
    filename=strcat(filename,frame_name,'.jpg');
    disp(filename);
    
    Epot_XY_map=griddata(XX,YY,Epot_2D,xxs,yys,'linear');
    Epot_XY_map=Epot_XY_map';
    cla(sp);
    imagesc((-0.5:0.001:0.5),(0:0.001:0.5),Epot_XY_map',[-0.5 0.5]*1e2);
    axis xy;
    title(strcat('Epot frame rank #',num2str(frame_name)));
%     surfc(XX,YY,Epot_2D,Cmapping,'edgecolor','none');view(0,90);

    ylim([0 0.3])
    xlim([-0.3 0.3])
    time_step=(f-1)
    
    figure(1)
    pause(0.4);


    
    F=getframe(fig);
[h, w, p] = size(F.cdata);  % use 1st frame to get dimensions
% hf = figure; 
% % resize figure based on frame's w x h, and place at (150, 150)
% set(hf, 'position', [150 150 w h]);
% axis off
    [im,map] = frame2im(F);    %Return associated image data
    imwrite(im,filename,'jpg');
    saveas(gcf,filename,'jpg');%saving freq plot
    
    
%     pause(0.4);

end
