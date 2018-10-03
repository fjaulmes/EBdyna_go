function n_rank = eval_dist_ref_axis_pot( Nlength, X_coords, Z_coords ,line_ref_pot_X, line_ref_pot_Z)
%take a point on the contuor and find it distance to the reference
%potential line
clear dist

for(n=1:Nlength)
    
    dist_values=sqrt((line_ref_pot_X-X_coords).^2+(line_ref_pot_Z-Z_coords).^2);
    
    dist(n)=min(dist_values)
    
end

[n_rank eps]=min(dist);

end

