load('wallRZ_CU_july2021.mat')

figure;
hold on
grid on
plot(wall_CU.R(1:578),wall_CU.Z(1:578))
plot(wall_CU.R(1436:end),wall_CU.Z(1436:end))

REMOVE_BUMP=9;

plot(wall_CU.R(578:753),wall_CU.Z(578:753))  % length is 175
plot(wall_CU.R(1261:1436),wall_CU.Z(1261:1436))  % length is 175

plot(wall_CU.R(578:753-REMOVE_BUMP),wall_CU.Z(578:753-REMOVE_BUMP))  % length is 175
plot(wall_CU.R(1261+REMOVE_BUMP:1436),wall_CU.Z(1261+REMOVE_BUMP:1436))  % length is 175

plot(wall_CU.R(753:1261),wall_CU.Z(753:1261))  
% REMOVE_BUMP=1
% plot(wall_CU.R(753+REMOVE_BUMP:1261-REMOVE_BUMP),wall_CU.Z(753+REMOVE_BUMP:1261-REMOVE_BUMP))  
% 
% wall_CU.R(754:754+REMOVE_BUMP)=[];
% wall_CU.Z(1261:1261+REMOVE_BUMP)=[];
% plot(wall_CU.R(753:1261),wall_CU.Z(753:1261))  

BUMP_IND1=753-REMOVE_BUMP;
BUMP_IND2=1261+REMOVE_BUMP;

R_lower_vessel = wall_CU.R(578:753-REMOVE_BUMP);
R_inner_vessel = wall_CU.R(754:1260); 
R_upper_vessel = wall_CU.R(1261+REMOVE_BUMP:1436);
R_outer_vessel = [ wall_CU.R(1436:end) ; wall_CU.R(1:578)  ];

Z_lower_vessel = wall_CU.Z(578:753-REMOVE_BUMP);
Z_inner_vessel = wall_CU.Z(754:1260);
Z_upper_vessel = wall_CU.Z(1261+REMOVE_BUMP:1436);
Z_outer_vessel = [wall_CU.Z(1436:end) ; wall_CU.Z(1:578)   ];

R_inner_vessel = [wall_CU.R(BUMP_IND1); R_inner_vessel ; wall_CU.R(BUMP_IND2)];
Z_inner_vessel = [wall_CU.Z(BUMP_IND1); Z_inner_vessel ; wall_CU.Z(BUMP_IND2)];

R_outer_vessel(272)=[];
Z_outer_vessel(272)=[];

%%
plot(R_lower_vessel(1:end-66),Z_lower_vessel(1:end-66))
R_inner_vessel=[R_lower_vessel(end-66:end-1) ; R_inner_vessel ; R_upper_vessel(2:68)];
Z_inner_vessel=[Z_lower_vessel(end-66:end-1) ; Z_inner_vessel ; Z_upper_vessel(2:68)];

R_lower_vessel=R_lower_vessel(1:end-66);
Z_lower_vessel=Z_lower_vessel(1:end-66);

R_upper_vessel=R_upper_vessel(68:end);
Z_upper_vessel=Z_upper_vessel(68:end);

plot(R_lower_vessel,Z_lower_vessel)
plot(R_inner_vessel,Z_inner_vessel)
plot(R_upper_vessel,Z_upper_vessel)


%%

% plot(wall_CU.R(1:63),wall_CU.Z(1:63))
% plot(wall_CU.R(63:157),wall_CU.Z(63:157))
% plot(wall_CU.R(157:219),wall_CU.Z(157:219))
% plot(wall_CU.R(219:end),wall_CU.Z(219:end))
% 
% R_lower_vessel = wall_CU.R(1:63);
% R_inner_vessel = wall_CU.R(63:157); 
% R_upper_vessel = wall_CU.R(157:220);
% R_outer_vessel = wall_CU.R(220:end);
% 
% Z_lower_vessel = wall_CU.Z(1:63);
% Z_inner_vessel = wall_CU.Z(63:157);
% Z_upper_vessel = wall_CU.Z(157:220);
% Z_outer_vessel = wall_CU.Z(220:end);
% 
% R_inner_vessel=R_inner_vessel(11:85);
% Z_inner_vessel=Z_inner_vessel(11:85);
% 
% R_lower_vessel=[R_lower_vessel(1:28) ; R_lower_vessel(38:end)];
% Z_lower_vessel=[Z_lower_vessel(1:28) ; Z_lower_vessel(38:end)];
% R_lower_vessel(29)=R_lower_vessel(29)-0.0001;
% 
% R_upper_vessel=[R_upper_vessel(1:26) ; R_upper_vessel(36:end)];
% Z_upper_vessel=[Z_upper_vessel(1:26) ; Z_upper_vessel(36:end)];
% R_upper_vessel(27)=R_upper_vessel(27)+0.0001;

% save('CU_vessel_limits.mat', 'R_inner_vessel', 'R_lower_vessel', 'R_outer_vessel', 'R_upper_vessel', 'Z_inner_vessel', 'Z_lower_vessel', 'Z_outer_vessel', 'Z_upper_vessel')
% 
% save('./CU_vessel_limits_standard_no_rescale.mat', 'R_inner_vessel', 'R_lower_vessel', 'R_outer_vessel', 'R_upper_vessel', 'Z_inner_vessel', 'Z_lower_vessel', 'Z_outer_vessel', 'Z_upper_vessel')

figure
hold on
plot(R_inner_vessel,Z_inner_vessel)
plot(R_lower_vessel,Z_lower_vessel)
plot(R_outer_vessel,Z_outer_vessel)
plot(R_upper_vessel,Z_upper_vessel)
