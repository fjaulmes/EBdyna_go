load('MASTU_WALL.mat')

figure
hold on
plot(MASTU_WALL.R(1:10),MASTU_WALL.Z(1:10))
plot(MASTU_WALL.R(10:11),MASTU_WALL.Z(10:11))
plot(MASTU_WALL.R(11:28),MASTU_WALL.Z(11:28))
plot(MASTU_WALL.R(28:end),MASTU_WALL.Z(28:end))

plot(MASTU_WALL.R(19:20),MASTU_WALL.Z(19:20))
MASTU_WALL.Z(19)=MASTU_WALL.Z(19)+0.01
MASTU_WALL.Z(20)=MASTU_WALL.Z(20)-0.01
plot(MASTU_WALL.R(11:28),MASTU_WALL.Z(11:28))

R_inner_vessel=[MASTU_WALL.R(1:10)];
Z_inner_vessel=[MASTU_WALL.Z(1:10)];

R_lower_vessel=MASTU_WALL.R(10:11);
Z_lower_vessel=MASTU_WALL.Z(10:11);

R_outer_vessel=MASTU_WALL.R(11:28);
Z_outer_vessel=MASTU_WALL.Z(11:28);


R_upper_vessel=MASTU_WALL.R(28:end);
Z_upper_vessel=MASTU_WALL.Z(28:end);




save('MASTU_vessel_limits.mat', 'R_inner_vessel','R_lower_vessel', 'R_outer_vessel', 'R_upper_vessel', 'Z_inner_vessel', 'Z_lower_vessel', 'Z_outer_vessel', 'Z_upper_vessel')
