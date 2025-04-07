%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%[R,Z,R2,Z2]=plot_COMPASS_vessel(ploti);
%
%Returns the coordinate of the vessel cross section of COMPASS
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUTS %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%ploti : optional : if different than 0, then plot the vessel
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OUTPUTS %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%R,Z : major radius and vertical position (m) of the main vessel cross section
%R2,Z2 : limiter positions


function [R,Z,R2,Z2 ph]=plot_CU_vessel(ploti)

ph=ploti;

if nargin==0
	ploti=1;
end

r_maj=0.894; 

load('wallRZ_CU.mat');

R=[wall_CU.R wall_CU.R(1)];
Z=[wall_CU.Z wall_CU.Z(1)];

% no info on limiters yet
R2=R;
Z2=Z;


if ploti

    hold on
	ph=plot(R,Z,'k','linewidth',3);
    hold off

%     plot(R2,Z2,'color',[0.0 0.0 0.0],'linewidth',3);

end

%     xlim([0.15 0.8]);ylim([-0.5 0.5])
