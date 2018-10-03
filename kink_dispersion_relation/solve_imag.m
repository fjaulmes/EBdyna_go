function [ value ] = solve_imag( RHS,real_val_norm,x0 )
value=interp1(real_val_norm,(imag(RHS)),x0,'*linear');
end

