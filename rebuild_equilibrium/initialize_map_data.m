

number_of_data_points=size(finesse_data,1);


%we need to recenter the X data that is on magnetic axis (Raxis)
%to major radius (R0)
% n=1;
% for(p=1:NP)
%     for(r=1:Nradial)
%         X_PR_map(p,r)=finesse_data(n,1);
%         Z_PR_map(p,r)=finesse_data(n,2);
%         n=n+1;
%     end
% end
% Zouter=Z_PR_map(:,end);
% [theta_pos_max_Zouter Zvalue]=max(Zouter);
% R0_pos=interp1(theta_scale,X_PR_map(:,end),theta_pos_max_Zouter);


% BE CAREFUL!
% FINESSE CREATES AN OFFSET ON ALL DATA for the X axis

% disp('FINESSE CREATES AN OFFSET ON ALL DATA for the X axis')
% disp('The Xaxis pos given by the output needs to be repositionned accordingly...')
% disp('The value of a has been rescaled...')
% Xoffset=((finesse_data(1,1))-X_axis)
% 
% a_recalc=0.5*(max(finesse_data(:,1))+abs(min(finesse_data(:,1))))
% 
% X_axis=X_axis+Xoffset
% 3.600308791889E-03  2.600409285518E-01

load finesse_tokamak_parameters.mat

epsilon=3.437000000000E-01;
% a=0.5*(abs(min(finesse_data(:,1)))+max(finesse_data(:,1)))
R0=a/epsilon
alpha=2.900000000000E+00 

psi1_input=a^2*Bphi0/alpha

X_axis=a*1.187240408016E-01 
Z_axis=a*(-6.145134411388E-03) 


plasma_beta_tot=3.600308791889E-03 
plasma_beta_pol=2.600409285518E-01



load input_data.mat

NPSI=size(input_data1,1);
DPSI=1/(NPSI-1);






