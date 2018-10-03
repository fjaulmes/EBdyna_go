function  [indexes,slopes] = interp_index_list_3D(size_map,Xindex_precise,Yindex_precise,Zindex_precise)
%%interp_index_list_2D Gives a list of the indexes of the square to
%%interpolate and the 'local' coordinates
% Use with interpolation function

%% Coordinate in index
% Find X and Z in equidistant grid in index value (+1 since value 0 is index 1)
% Xindex_precise=(x*DX_inv)+1;
% Yindex_precise=(y*DY_inv);
% Zindex_precise=(z*DZ_inv)+1;

% Limit 'position' used for interpolation inside grid
Xindex_precise(Xindex_precise<1)=1;
Yindex_precise(Yindex_precise<1)=1;
Zindex_precise(Zindex_precise<1)=1;
Xindex_precise(Xindex_precise>size_map(1))=size_map(1);
Yindex_precise(Yindex_precise>size_map(2))=size_map(2);
Zindex_precise(Zindex_precise>size_map(3))=size_map(3);

%% Left corner of square
Xindex=floor(Xindex_precise);
Yindex=floor(Yindex_precise);
Zindex=floor(Zindex_precise);

% If outside of grid, use square to the left / bottom
expr_X=Xindex==size_map(1);
expr_Y=Yindex==size_map(2);
expr_Z=Zindex==size_map(3);
Xindex(expr_X)=Xindex(expr_X)-1; 
Yindex(expr_Y)=Yindex(expr_Y)-1; 
Zindex(expr_Z)=Zindex(expr_Z)-1;

%% Indexes of square in which to interpolate
% indexes.ldb = sub2ind(size_map, Xindex  ,Yindex,  Zindex  ); %left  down  bottom corner
% indexes.ldt = sub2ind(size_map, Xindex  ,Yindex,  Zindex+1); %left  down  top    corner
% indexes.lub = sub2ind(size_map, Xindex  ,Yindex+1,Zindex  ); %left  upper bottom corner
% indexes.lut = sub2ind(size_map, Xindex  ,Yindex+1,Zindex+1); %left  upper top    corner
% indexes.rdb = sub2ind(size_map, Xindex+1,Yindex,  Zindex  ); %right down  bottom corner
% indexes.rdt = sub2ind(size_map, Xindex+1,Yindex,  Zindex+1); %right down  top    corner
% indexes.rub = sub2ind(size_map, Xindex+1,Yindex+1,Zindex  ); %right upper bottom corner
% indexes.rut = sub2ind(size_map, Xindex+1,Yindex+1,Zindex+1); %right upper top    corner

size_1=size_map(1);
size_21=size_map(2)*size_1;

indexes.ldb = Xindex  +(Yindex-1)*size_1+(Zindex-1)*size_21;  %left  down  bottom corner
indexes.ldt = Xindex  +(Yindex-1)*size_1+(Zindex  )*size_21;  %left  down  top    corner
indexes.lub = Xindex  +(Yindex  )*size_1+(Zindex-1)*size_21; %left  upper bottom corner
indexes.lut = Xindex  +(Yindex  )*size_1+(Zindex  )*size_21; %left  upper top    corner
indexes.rdb = Xindex+1+(Yindex-1)*size_1+(Zindex-1)*size_21; %right down  bottom corner
indexes.rdt = Xindex+1+(Yindex-1)*size_1+(Zindex  )*size_21; %right down  top    corner
indexes.rub = Xindex+1+(Yindex  )*size_1+(Zindex-1)*size_21; %right upper bottom corner
indexes.rut = Xindex+1+(Yindex  )*size_1+(Zindex  )*size_21; %right upper top    corner

%% local coordinate (from 0 to 1)
slopes.x=Xindex_precise-Xindex;
slopes.y=Yindex_precise-Yindex;
slopes.z=Zindex_precise-Zindex;