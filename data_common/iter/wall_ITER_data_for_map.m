load('wallRZ_ITER.mat')

figure
hold on
plot(wall_ITER.R(1:10),wall_ITER.Z(1:10))
plot(wall_ITER.R(10:12),wall_ITER.Z(10:12))
plot(wall_ITER.R(12:33),wall_ITER.Z(12:33))
plot(wall_ITER.R(33:43),wall_ITER.Z(33:43))
plot(wall_ITER.R(43:end-1),wall_ITER.Z(43:end-1))


R_outer_vessel=[wall_ITER.R(43:end-1) wall_ITER.R(1:10)];
Z_outer_vessel=[wall_ITER.Z(43:end-1) wall_ITER.Z(1:10)];

R_upper_vessel=wall_ITER.R(10:12);
Z_upper_vessel=wall_ITER.Z(10:12);

R_inner_vessel=[wall_ITER.R(12:20) wall_ITER.R(22:33)];
Z_inner_vessel=[wall_ITER.Z(12:20) wall_ITER.Z(22:33)];

R_lower_vessel=wall_ITER.R(33:43);
Z_lower_vessel=wall_ITER.Z(33:43);


save('ITER_vessel_limits.mat', 'R_inner_vessel','R_lower_vessel', 'R_outer_vessel', 'R_upper_vessel', 'Z_inner_vessel', 'Z_lower_vessel', 'Z_outer_vessel', 'Z_upper_vessel')
