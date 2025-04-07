load('wallRZ_CU.mat')

openfig('wallRZ_CU_smoothed.fig')
hw = findobj(gca,'Type','line');
xw=get(hw,'Xdata') ;
zw=get(hw,'Ydata') ;

hold on
grid on
plot(wall_CU.R,wall_CU.Z,'r--')

R_upper_vessel=xw(abs(zw-0.58)<1e-4);
Z_upper_vessel=zw(abs(zw-0.58)<1e-4);
plot(R_upper_vessel,Z_upper_vessel,'linewidth',3)
xw(abs(zw-0.58)<1e-4)=[];
zw(abs(zw-0.58)<1e-4)=[];

R_outer_vessel=xw(xw>0.81-1e-4);
Z_outer_vessel=zw(xw>0.81-1e-4);
plot(R_outer_vessel,Z_outer_vessel,'linewidth',3)
zw(xw>0.81-1e-4)=[];
xw(xw>0.81-1e-4)=[];

R_inner_vessel=xw(xw<0.655);
Z_inner_vessel=zw(xw<0.655);
R_inner_vessel(Z_inner_vessel<-0.6+1e-4)=[];
Z_inner_vessel(Z_inner_vessel<-0.6+1e-4)=[];
[Z_inner_vessel SORT_IND]=sort(Z_inner_vessel);
R_inner_vessel=R_inner_vessel(SORT_IND);
plot(R_inner_vessel,Z_inner_vessel,'linewidth',3)
for zz=1:length(Z_inner_vessel)
    INN_WALL=(xw==R_inner_vessel(zz))&(zw==Z_inner_vessel(zz));
    zw(INN_WALL)=[];
    xw(INN_WALL)=[];
end

R_lower_vessel=xw;
Z_lower_vessel=zw;
plot(R_lower_vessel,Z_lower_vessel,'linewidth',3)


save('CU_vessel_limits.mat', 'R_inner_vessel', 'R_lower_vessel', 'R_outer_vessel', 'R_upper_vessel', 'Z_inner_vessel', 'Z_lower_vessel', 'Z_outer_vessel', 'Z_upper_vessel')

save('./CU_vessel_limits_standard_no_rescale.mat', 'R_inner_vessel', 'R_lower_vessel', 'R_outer_vessel', 'R_upper_vessel', 'Z_inner_vessel', 'Z_lower_vessel', 'Z_outer_vessel', 'Z_upper_vessel')

% figure
% hold on
% plot(R_inner_vessel,Z_inner_vessel)
% plot(R_lower_vessel,Z_lower_vessel)
% plot(R_outer_vessel,Z_outer_vessel)
% plot(R_upper_vessel,Z_upper_vessel)
