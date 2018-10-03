function matlab2tikz_fast(hf,filename,format)
%matlab2tikz_fast Have matlab2tikz structure, but with pdf of image
%   Matlab2tikz is called with imageAsPng option. Instead of the PNG one
%   can use the produced PDF names FIGURE < date > 
%   This works for one axes only, namely 

debug_name='matlab2tikz_fast_figure.tikz';

if nargin==0
    hf=gcf;
    filename=debug_name;
    format='pdf';
elseif ishandle(hf)
switch nargin
    case 1
        filename=debug_name;
        format='pdf';
    case 2
        format='pdf';
    case 3
        if isempty(filename); filename=debug_name; end
end
else
switch nargin
    case 1
        format='pdf';
        filename=hf;
        hf=gcf;
    case 2
        format=filename;
        filename=hf;
        hf=gcf;
end
end         
%% Get proper handle and handle of parent
if ~ishandle(hf)
    filename=hf;
    hf=gcf;
elseif isempty(hf)
    hf=gcf;
end
has=findall(get(hf,'children'),'type','Axes');
nr_has=length(has);

%% Put all data in other figure to prevent it from being copied in matlab2tikz
% Make a zombie figure to copy the lines temperorely
hf_zombie{nr_has}=[];
lines{nr_has}=[];
ha_zombie{nr_has}=[];
h_image{nr_has}=[];
for i=1:nr_has
    % Axes of 'live' figure
    ha=has(i);
    lines{i}=get(ha,'children');
    
    % First make sure the axis system does not change
    axis(ha,'manual')
    hold(ha,'on')
    set(ha,'ClimMode','manual');

    % Copy the lines to a zombie figure
    hf_zombie{i}=figure('windowstyle','normal','units','normalized','outerposition',[0 0 1 1]);
    ha_zombie{i}=axes('parent',hf_zombie{i});
    set(flipud(lines{i}),'parent',ha_zombie{i});
    
    % Have the same axes limit for the zombie to have the correct export to
    % the figure
    set(ha_zombie{i},'Xlim',get(ha,'Xlim'),'Ylim',get(ha,'Ylim'),'Zlim',get(ha,'Zlim'),'Clim',get(ha,'Clim'));
    % Other display properties
    colormap(ha_zombie{i},colormap(ha))
    set(ha_zombie{i},'XDir',get(ha,'XDir'),'YDir',get(ha,'YDir'),'ZDir',get(ha,'ZDir'));
    
    % Make axes the size of the figure    
    set(ha_zombie{i},'Visible','off') % Remove axes
    set(ha_zombie{i},'units','normalized')
    set(ha_zombie{i},'position',[0 0 1 1]) % Make size of figure
    drawnow;
    set(hf_zombie{i},'units','inches') % Cannot be normalized, so take any unit
    pos = get(hf_zombie{i},'Position');
    set(hf_zombie{i},'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    drawnow

    %% Export to pdf
    points=strfind(filename,'.'); 
    if isempty(points); points=length(filename)+1; end
    filename_i=strcat(filename(1:(points-1)),'--',num2str(i),'.',format);
    if any(strcmp(format,{'pdf','eps'}))
        warning('Set to OpenGL instead of painters since vector printing by ML can lead to white lines.')
    end
    print(hf_zombie{i},filename_i,['-d',format],'-r0','-opengl')
    
    drawnow

    %% Make a 'fake' image, to get the export to png function of matlab2tikz
    h_image{i}=imagesc(get(ha,'XLim'),get(ha,'YLim'),NaN(2,2),'parent',ha);
end

%% Matlab2tikz this stuff
matlab2tikz(filename,'figurehandle',hf,'imagesAsPng',true,'height', '\figureheight', 'width', '\figurewidth' );

for i=1:nr_has
    ha=has(i);
    set(flipud(lines{i}),'parent',ha);
    delete(hf_zombie{i})
    delete(h_image{i})
end

disp('*** SUCCESS!! MAKE SURE TO RENAME THE IMPORT FILE IN THE TIKZ FILE TO THE PDF!')
end
