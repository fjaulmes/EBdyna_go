function axis_rank = eval_intersection_ref_axis_pot(  X_coord, Z_coord ,line_ref_pot_X, line_ref_pot_Z)
%take a point on the contuor and find it distance to the reference
%potential line

 
dist_values=sqrt((line_ref_pot_X-X_coord).^2+(line_ref_pot_Z-Z_coord).^2);
    


[eps axis_rank ]=min(dist_values);

end

