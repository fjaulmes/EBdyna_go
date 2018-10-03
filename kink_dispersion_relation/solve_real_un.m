function [ value ] = solve_real( RHS,imag_values,x0 )
value=interp1(imag_values,(real(RHS)-1),x0,'*linear');
end

