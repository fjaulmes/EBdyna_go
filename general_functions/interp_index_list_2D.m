function [indexes,slopes] = interp_index_list_2D (size_map,Xindex_precise,Zindex_precise)
%%interp_index_list_2D Gives a list of the indexes of the square to
%%interpolate and the 'local' coordinates

% Use with interpolation function
size_1=size_map(1);
%% Coordinate in index
% Find X and Z in equidistant grid in index value
% Xindex_precise=((x(:,1,:)-dim.R0)/dim.DX)+dim.mid_Xzero;
% Zindex_precise=(x(:,2,:)/dim.DX)+dim.mid_Z;

% Limit 'position' used for interpolation inside grid
Xindex_precise(Xindex_precise<1)=1;
Zindex_precise(Zindex_precise<1)=1;
Xindex_precise(Xindex_precise>size_1)=size_1;
Zindex_precise(Zindex_precise>size_map(2))=size_map(2);

%% Left corner of square
Xindex=floor(Xindex_precise);
Zindex=floor(Zindex_precise);

% If outside of grid, use square to the left / bottom
expr_X=Xindex==size_1;
expr_Z=Zindex==size_map(2);
Xindex(expr_X)=Xindex(expr_X)-1; 
Zindex(expr_Z)=Zindex(expr_Z)-1;

%% Indexes of square in which to interpolate
% indexes.lb = sub2ind(size_map, Xindex  , Zindex  ); % left  bottom corner
% indexes.lt = sub2ind(size_map, Xindex  , Zindex+1); % left  top    corner
% indexes.rb = sub2ind(size_map, Xindex+1, Zindex  ); % right bottom corner
% indexes.rt = sub2ind(size_map, Xindex+1, Zindex+1); % right top    corner

indexes.lb = Xindex+  (Zindex-1)*size_1; % left  bottom corner
indexes.lt = Xindex+  (Zindex  )*size_1; % left  top    corner
indexes.rb = Xindex+1+(Zindex-1)*size_1; % right bottom corner
indexes.rt = Xindex+1+(Zindex  )*size_1; % right top    corner

%% local coordinate (from 0 to 1)
slopes.x=Xindex_precise-Xindex;
slopes.z=Zindex_precise-Zindex;
end