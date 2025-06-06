% Script: map_psi_FIESTA_streamlined.m
% Purpose: Generate required inputs for EBdyna from CDB and FIESTA
% Outputs: pressure_profile.mat, motions_map_dimensions.mat, XZsmall_fields_tokamak_pre_collapse.mat

clc; clearvars; close all;

%% === Configuration ===
SHOT_NUMBER = 5400;
TIME_EQUIL = 1.15;         % Equilibrium time in seconds
sim_folder_name='CU_5400_1150_RT42_01'
DATA_PLASMA_FOLDER=['../',sim_folder_name,'/data_plasma/']
DATA_PHYS_FOLDER=['../data_common/physics_data/']

SAVEFILE = true;
Rmin = 0.59; Rmax = 1.24;   % Grid boundaries for interpolation
R0 = 0.894
NX = 400; NZ = 800;       % Grid resolution
Delta_Z=NZ*(Rmax-Rmin)/NX;
Zmin = -0.5*Delta_Z; Zmax = 0.5*Delta_Z;
interp_method = 'spline'; % Options: 'gridfit', 'scattered', 'spline'

% extract FIESTA object and save it
extract_FIESTA_from_cudb;

%% === Load constants ===
load([DATA_PHYS_FOLDER,'physics_constants.mat']);

%% === Load FIESTA structure ===
load([DATA_PLASMA_FOLDER,'FIESTA_equil.mat']);

%% === CDB access to 1D profiles from OMP_PROFILES (with n0/n0_D2 fallback) ===
cdb = cdb_client();
fields = {'Te', 'Ti', 'ne', 'ni', 'vtor', 'R','psi'};
for i = 1:length(fields)
    signal = cdb.get_signal([fields{i} '/OMP_PROFILES:' num2str(SHOT_NUMBER)]);
    data_struct.(fields{i}) = interp1(signal.time_axis.data, signal.data, TIME_EQUIL);
    if strcmp(fields{i}, 'R')
        psi_norm = signal.axis1.data;
    end
end

% Try to get optional neutral density profiles
try
    signal = cdb.get_signal(['n0/OMP_PROFILES:' num2str(SHOT_NUMBER)]);
    data_struct.n0 = interp1(signal.time_axis.data, signal.data, TIME_EQUIL);
catch
    warning('n0 profile not found in OMP_PROFILES.');
end

try
    signal = cdb.get_signal(['n0_D2/OMP_PROFILES:' num2str(SHOT_NUMBER)]);
    data_struct.n0_D2 = interp1(signal.time_axis.data, signal.data, TIME_EQUIL);
catch
    warning('n0_D2 profile not found in OMP_PROFILES.');
end

Te_prof = data_struct.Te;
Ti_prof = data_struct.Ti;
ne_prof = data_struct.ne;
ni_prof = data_struct.ni;
vtor_prof = data_struct.vtor;
psi_scale = data_struct.psi;

%% === Interpolated 2D Maps ===
map2D = FIESTA.map2D;
[Zg, Rg] = meshgrid(map2D.scale_Z, map2D.scale_R);

xnodes = linspace(Rmin, Rmax, NX);
ynodes = linspace(Zmin, Zmax, NZ);
%[X_grid, Z_grid] = meshgrid(xnodes, ynodes);
[Z_grid, X_grid] = meshgrid(ynodes, xnodes);
% Z_grid=Zg;
% X_grid=Rg;

switch lower(interp_method)
    case 'scattered'
        % F_pressure = scatteredInterpolant(Rg(:), Zg(:), map2D.pressure(:), 'natural', 'none');
        F_psi      = scatteredInterpolant(Rg(:), Zg(:), map2D.psi(:), 'natural', 'none');
        % F_psi_n    = scatteredInterpolant(Rg(:), Zg(:), map2D.psi_n(:), 'natural', 'none');
        % F_Bphi     = scatteredInterpolant(Rg(:), Zg(:), map2D.Bphi(:), 'natural', 'none');
        % F_Br       = scatteredInterpolant(Rg(:), Zg(:), map2D.Br(:), 'natural', 'none');
        % F_Bz       = scatteredInterpolant(Rg(:), Zg(:), map2D.Bz(:), 'natural', 'none');

        % pressure_map = F_pressure(X_grid, Z_grid);
        psi_map      = F_psi(X_grid, Z_grid);
        [dpsidZ_map, dpsidR_map] = gradient(psi_map, Z_grid, R_grid);

        % psi_n_map    = F_psi_n(X_grid, Z_grid);
        % Bphi_map     = F_Bphi(X_grid, Z_grid);
        % Br_map       = F_Br(X_grid, Z_grid);
        % Bz_map       = F_Bz(X_grid, Z_grid);
    case 'spline'
        X_vec = X_grid(:);
        Y_vec = Z_grid(:);

        sp_psi = csapi({FIESTA.map2D.scale_R, FIESTA.map2D.scale_Z}, FIESTA.map2D.psi);
        psi_vec = fnval(sp_psi, [X_vec'; Y_vec']);  % Must be a matrix with 2 lines
        psi_map = reshape(psi_vec, size(X_grid));        
        %sp_br = csapi({FIESTA.map2D.scale_R, FIESTA.map2D.scale_Z}, FIESTA.map2D.Br);
        %Br_map = reshape(fnval(sp_br, [X_vec'; Y_vec']), size(X_grid));        
        %sp_bz = csapi({FIESTA.map2D.scale_R, FIESTA.map2D.scale_Z}, FIESTA.map2D.Bz);
        %Bz_map = reshape(fnval(sp_bz, [X_vec'; Y_vec']), size(X_grid));        
        % Compute smooth partial derivatives using fnder
        sp_dpsi_dR = fnder(sp_psi, [1 0]);  % ∂ψ/∂R
        dpsidR_map = reshape(fnval(sp_dpsi_dR, [X_vec'; Y_vec']), size(X_grid));        
        sp_dpsi_dZ = fnder(sp_psi, [0 1]);  % ∂ψ/∂Z
        dpsidZ_map = reshape(fnval(sp_dpsi_dZ, [X_vec'; Y_vec']), size(X_grid));        

    case 'gridfit'
        gf_data=reshape(FIESTA.map2D.psi,1,length(FIESTA.map2D.scale_R)*length(FIESTA.map2D.scale_Z));
        FIESTA_data_X_gf=repmat(FIESTA.map2D.scale_R,1,length(FIESTA.map2D.scale_Z))-R0;
        FIESTA_data_Z_gf = reshape(repmat(FIESTA.map2D.scale_Z, length(FIESTA.map2D.scale_R), 1), [], 1);
        psi_map      = gridfit(FIESTA_data_X_gf,FIESTA_data_Z_gf,gf_data,xnodes,ynodes,'smoothness',GRIDFIT_SMOOTHNESS);
        [dpsidZ_map, dpsidR_map] = gradient(psi_map, Z_grid, R_grid);

        % pressure_map = gridfit(Rg(:), Zg(:), map2D.pressure(:), xnodes, ynodes);
        % psi_map      = gridfit(Rg(:), Zg(:), map2D.psi(:), xnodes, ynodes);
        % psi_n_map    = gridfit(Rg(:), Zg(:), map2D.psi_n(:), xnodes, ynodes);
        % Bphi_map     = gridfit(Rg(:), Zg(:), map2D.Bphi(:), xnodes, ynodes);
        % Br_map       = gridfit(Rg(:), Zg(:), map2D.Br(:), xnodes, ynodes);
        % Bz_map       = gridfit(Rg(:), Zg(:), map2D.Bz(:), xnodes, ynodes);

    otherwise
        error('Unknown interpolation method: %s', interp_method);
end

%% === Coordinates for EBdyna grid ===
R_grid = xnodes;
X_grid = R_grid-R0;
Z_grid = ynodes;
NR = NX;
NZ = NZ;
DX = mean(diff(R_grid));
DZ = mean(diff(Z_grid));
mid_X = round(NR / 2);
mid_Z = round(NZ / 2);
Rpos = R_grid;
Rpos_XZ_map = repmat(Rpos(:)', NZ, 1)';
Br_map = -dpsidZ_map ./ Rpos_XZ_map / (2*pi);
Bz_map =  dpsidR_map ./ Rpos_XZ_map / (2*pi);
R_axis=FIESTA.r_mag;
Z_axis=FIESTA.z_mag;

%% Secondary 2D maps derived from psi
psi_axis = FIESTA.psi_axis;           % in Wb
psi_boundary = FIESTA.psi_boundary;   % in Wb
% Compute normalized poloidal flux
psi_n_map = (psi_map - psi_axis) ./ (psi_boundary - psi_axis);
% Mask of points below the X-point (private flux region)
private_flux_mask = Z_grid < FIESTA.Z_xpoint;
% Force psi_n = 1 in the private flux region
psi_n_map=private_flux_mask+psi_n_map.*(1-private_flux_mask);
psi_n_map = min(max(psi_n_map, 0), 1);

% Interpolate F onto the 2D psi_n grid
F_XZ_map = interp1(FIESTA.prof.psiN, FIESTA.prof.f, psi_n_map, 'makima', 'extrap');
% Compute Bphi using Bphi = F / R
Bphi_map = F_XZ_map ./ Rpos_XZ_map;

%% === Save pressure_profile.mat ===
if SAVEFILE
    if isfield(data_struct, 'n0')
        n0_prof = data_struct.n0;
    else
        n0_prof = [];
    end
    if isfield(data_struct, 'n0_D2')
        n0_D2_prof = data_struct.n0_D2;
    else
        n0_D2_prof = [];
    end
    save([DATA_PLASMA_FOLDER,'pressure_profile.mat'], 'psi_norm','psi_scale', 'Te_prof', 'Ti_prof', 'ne_prof', 'ni_prof', 'vtor_prof', 'n0_prof', 'n0_D2_prof');
end

%% === Save motions_map_dimensions.mat ===
if SAVEFILE
    save([DATA_PLASMA_FOLDER,'motions_map_dimensions.mat'], 'R0','R_axis','Z_axis','NR', 'NZ', 'X_grid', 'Z_grid', 'Rpos_XZ_map', 'DX', 'DZ', 'mid_X', 'mid_Z');
end

%% === Save XZsmall_fields_tokamak_pre_collapse.mat ===
if SAVEFILE
    save([DATA_PLASMA_FOLDER,'XZsmall_fields_tokamak_pre_collapse.mat'],  'psi_map', 'psi_n_map', 'Bphi_map', 'Br_map', 'Bz_map', 'R_grid', 'Z_grid');
end

disp('EBdyna input files generated.')
