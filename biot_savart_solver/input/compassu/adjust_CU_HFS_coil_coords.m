%
% assumes you have imported the 2 colum2s as X_coils and Z_coils


for ind_point=1:length(X_coils)
    if X_coils(ind_point)<0.43
        X_coils(ind_point)=X_coils(ind_point)-0.16;
    end
end

figure
grid on 
hold on
plot(X_coils,Z_coils)
save CU_coil_XZ_coords X_coils Z_coils