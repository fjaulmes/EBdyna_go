function h=vectarrow(p0,p1,handle,color_line,color_head)
%Arrowline 3-D vector plot.
%   vectarrow(p0,p1) plots a line vector with arrow pointing from point p0
%   to point p1. The function can plot both 2D and 3D vector with arrow
%   depending on the dimension of the input
%
%   Example:
%       3D vector
%       p0 = [1 2 3];   % Coordinate of the first point p0
%       p1 = [4 5 6];   % Coordinate of the second point p1
%       vectarrow(p0,p1)
%
%       2D vector
%       p0 = [1 2];     % Coordinate of the first point p0
%       p1 = [4 5];     % Coordinate of the second point p1
%       vectarrow(p0,p1)
%
%   See also Vectline

%   Rentian Xiong 4-18-05
%   $Revision: 1.0

if exist('handle','var') && ~isempty(handle)
    h_axis=handle;
else
    h_axis=gca;
end
washold=ishold(h_axis);

if max(size(p0))==3
    if max(size(p1))==3
        x0 = p0(1);
        y0 = p0(2);
        z0 = p0(3);
        x1 = p1(1);
        y1 = p1(2);
        z1 = p1(3);
        if exist('color_line','var') && ~isempty(color_line)
            h(1)=plot3(h_axis,[x0;x1],[y0;y1],[z0;z1],'color',color_line);   % Draw a line between p0 and p1
        else
            h(1)=plot3(h_axis,[x0;x1],[y0;y1],[z0;z1],'color','r');   % Draw a line between p0 and p1
        end
        
        p = p1-p0;
        alpha = 0.1;  % Size of arrow head relative to the length of the vector
        beta = 0.1;  % Width of the base of the arrow head relative to the length
        
        hu = [x1-alpha*(p(1)+beta*(p(2)+eps)); x1; x1-alpha*(p(1)-beta*(p(2)+eps))];
        hv = [y1-alpha*(p(2)-beta*(p(1)+eps)); y1; y1-alpha*(p(2)+beta*(p(1)+eps))];
        hw = [z1-alpha*p(3);z1;z1-alpha*p(3)];
        
        
        hold(h_axis,'on')
        if exist('color_head','var') && ~isempty(color_head)
            h(2)=plot3(h_axis,hu(:),hv(:),hw(:),'color',color_head);  % Plot arrow head
        else
            h(2)=plot3(h_axis,hu(:),hv(:),hw(:),'color','r');  % Plot arrow head
        end
        if ~washold
            hold(h_axis,'off');
        end
    else
        error('p0 and p1 must have the same dimension')
    end
elseif max(size(p0))==2
    if max(size(p1))==2
        x0 = p0(1);
        y0 = p0(2);
        x1 = p1(1);
        y1 = p1(2);
        if exist('color_line','var') && ~isempty(color_line)
            h(1)=plot(h_axis,[x0;x1],[y0;y1],'color',color_line,'linewidth',1.5);   % Draw a line between p0 and p1
        else
            h(1)=plot(h_axis,[x0;x1],[y0;y1],'color',[0.5 0.5 0.5],'linewidth',1.5);   % Draw a line between p0 and p1
        end
        
        p = p1-p0;
        alpha = 0.1;  % Size of arrow head relative to the length of the vector
        beta = 0.1;  % Width of the base of the arrow head relative to the length
        
        hu = [x1-alpha*(p(1)+beta*(p(2)+eps)); x1; x1-alpha*(p(1)-beta*(p(2)+eps))];
        hv = [y1-alpha*(p(2)-beta*(p(1)+eps)); y1; y1-alpha*(p(2)+beta*(p(1)+eps))];
        
        hold on
        if exist('color_head','var') && ~isempty(color_head)
            h(2)=plot(h_axis,hu(:),hv(:),'color',color_head,'linewidth',2);  % Plot arrow head
        else
            h(2)=plot(h_axis,hu(:),hv(:),'color','r','linewidth',2);  % Plot arrow head
        end
        if ~washold
            hold(h_axis,'off')
        end
    else
        error('p0 and p1 must have the same dimension')
    end
else
    error('this function only accepts 2D or 3D vector')
end
