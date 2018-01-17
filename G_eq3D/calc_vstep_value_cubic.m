function [ v_step ] = calc_vstep_value_cubic(v_prev_prev, v_prev,v,v_next )
%estimate step value through third order polynomial interpolation
dY1=v_prev_prev-v_prev;
dY2=v-v_prev;
dY3=v_next-v_prev;

%MATPOLY=[[-(1/6) -0.5 +(1/6)][0.5 0.5 0][-(1/3) 1.0 -(1/6)]];
%A=-(1/6)*dY1-0.5*dY2+(1/6)*dY3;
%B=0.5*dY1+0.5*dY2;
%C=-(1/3)*dY1+1.0*dY2-(1/6)*dY3;
%dYs=A*(1.5)^3+B*(1.5)^2+C*1.5;

A=-3.375/6+2.25/2-0.5;
B=-0.5*3.375+0.5*2.25+1.5;
C=3.375/6-1.5/6;

dYs=A*dY1+B*dY2+C*dY3;

v_step=v_prev+dYs;

end

