
disp('Finding the reference middle potential line for phi=');
disp(phi);

%rx_precise=xPsih_zero*rx_value;
if (phi<=pi)
    omega=(NP-1)*(phi)/(2*pi)+1;
    omega_op=(NP-1)*(pi+phi)/(2*pi)+1;
else
    omega=(NP-1)*(phi)/(2*pi)+1;
    omega_op=(NP-1)*(phi-pi)/(2*pi)+1;
end

if (phi<=pi)
    Line_ref_X=[interp2((1:NP),(1:Nradial),X_PR_map',omega_op,(Nradial:-1:1)) ; interp2((1:NP),(2:Nradial),X_PR_map(:,2:Nradial)',omega,(2:Nradial))];
    Line_ref_Z=[interp2((1:NP),(1:Nradial),Z_PR_map',omega_op,(Nradial:-1:1)) ; interp2((1:NP),(2:Nradial),Z_PR_map(:,2:Nradial)',omega,(2:Nradial))];
else
    Line_ref_X=[interp2((1:NP),(1:Nradial),X_PR_map',omega,(Nradial:-1:1)) ; interp2((1:NP),(2:Nradial),X_PR_map(:,2:Nradial)',omega_op,(2:Nradial))];
    Line_ref_Z=[interp2((1:NP),(1:Nradial),Z_PR_map',omega,(Nradial:-1:1)) ; interp2((1:NP),(2:Nradial),Z_PR_map(:,2:Nradial)',omega_op,(2:Nradial))];
end

if DISPLAY_OUTPUTS==1
    hold on
    plot(Line_ref_X,Line_ref_Z,'r');
    pause(0.1);
end