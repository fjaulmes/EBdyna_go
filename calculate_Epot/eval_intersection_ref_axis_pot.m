function n_rank = eval_intersection_ref_axis_pot( Nbegin,Nend, X_coords, Z_coords ,line_ref_pot_X, line_ref_pot_Z)
%take a point on the contuor and find it distance to the reference
%potential line
dist=0;

for(n=Nbegin:Nend)
    
    dist_values=sqrt((line_ref_pot_X-X_coords(n)).^2+(line_ref_pot_Z-Z_coords(n)).^2);
    
    dist(n-Nbegin+1)=min(dist_values);
    
end

[eps n_rank ]=min(dist);

end

