function map_output = smooth_small_map( map_input )

% diag_coef=2*(1-sqrt(2)/2);
diag_coef=0.5;
filter_cell=[diag_coef 1 diag_coef ; 1 2 1; diag_coef 1 diag_coef ];
map_output=conv2(map_input,filter_cell,'same')/sum(sum(filter_cell));


end

