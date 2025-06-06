% Script: validate_psi_and_B_fields.m
% Purpose: Check correctness of psi interpolation and ∇⋅B = 0

clc;  close all;

%% === Load interpolated maps ===
% load('./data_tokamak/XZsmall_fields_tokamak_pre_collapse.mat');
% load('./data_tokamak/motions_map_dimensions.mat');
% 
% mu0 = 4*pi*1e-7;

%% === Compute Gradients ===
[dPsi_dZ, dPsi_dR] = gradient(psi_map, Z_grid, R_grid);
% dPsi_dR = gradient_spline(psi_map, R_grid);
% dPsi_dZ = gradient_spline(psi_map, Z_grid);

% Compute B_r and B_z from psi (note minus signs due to ∇ x A)
Br_check = -dpsidZ_map ./ Rpos_XZ_map' / (2*pi);
Bz_check =  dpsidR_map ./ Rpos_XZ_map' / (2*pi);

%% === Divergence of B ===
[~ , dBr_dR] = gradient(Br_check, Z_grid, R_grid);
[dBz_dZ , ~] = gradient(Bz_check, Z_grid, R_grid);

divB = dBr_dR + dBz_dZ;

% For FIESTA data
[~ , dBr_dR] = gradient(FIESTA.map2D.Br, FIESTA.map2D.scale_Z, FIESTA.map2D.scale_R);
[dBz_dZ , ~] = gradient(FIESTA.map2D.Bz, FIESTA.map2D.scale_Z, FIESTA.map2D.scale_R);

divB_FIESTA = dBr_dR + dBz_dZ;

%% === Plot Comparison of B fields ===
figure;
subplot(2,2,1);
imagesc(R_grid, Z_grid, Br_map' - Br_check',[-1e-3,1e-3]);
axis xy; colorbar; title('Br error: interpolated - derived'); xlabel('R'); ylabel('Z');

subplot(2,2,2);
imagesc(R_grid, Z_grid, Bz_map' - Bz_check',[-1e-3,1e-3]);
axis xy; colorbar; title('Bz error: interpolated - derived'); xlabel('R'); ylabel('Z');

subplot(2,2,3);
imagesc(R_grid, Z_grid, divB',[-1,1]);
axis xy; colorbar; title('Divergence of B (gradient psi)'); xlabel('R'); ylabel('Z');

subplot(2,2,4);
imagesc(FIESTA.map2D.scale_R, FIESTA.map2D.scale_Z, divB_FIESTA',[-1,1]);
axis xy; colorbar; title('Divergence of B FIESTA'); xlabel('R'); ylabel('Z');

%% === Statistics ===
fprintf('Max |Br_error| = %g\n', max(abs(Br_map(:) - Br_check(:))));
fprintf('Max |Bz_error| = %g\n', max(abs(Bz_map(:) - Bz_check(:))));
fprintf('Max |divB|     = %g\n', max(abs(divB(:))));