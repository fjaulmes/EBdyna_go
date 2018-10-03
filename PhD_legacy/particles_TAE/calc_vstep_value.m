function [ v_step ] = calc_vstep_value( v_prev,v,v_next )
%estimate step value through second order polynomial interpolation
dY1=v-v_prev;
dY2=v_next-v_prev;

A=-dY1+0.5*dY2;
B=2*dY1-0.5*dY2;

dYs=A*(1.5)^2+B*(1.5);

v_step=v_prev+dYs;

end

