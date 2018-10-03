function map_output = smooth_small_map( map_input )

filter_cell=[0.25 1 0.25 ; 1 2 1; 0.25 1 0.25 ];
map_output=conv2(map_input,[0.25 1 0.25 ; 1 2 1; 0.25 1 0.25 ],'same')/sum(sum(filter_cell));


end

