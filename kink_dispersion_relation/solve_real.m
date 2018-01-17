function [ value ] = solve_real( RHS,real_values,x0 )
value=interp1(real_values,(real(RHS)-1),x0,'*linear');
end

