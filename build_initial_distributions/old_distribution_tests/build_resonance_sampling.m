function build_resonance_sampling
%build_phase_space_sampling Builds evenly spaced sampling in 3D phase space
%   Loads in plasma parameters and creates particles which start on
%   Z=Z_axis of equilibrium.

seed=1;
% seed='shuffle';
rng(seed);
%% Load default parameters
par.paths=initialize_folder_names_struct;   % Folder names
[const,maps,dim]=load_distr_maps(par);
%% Parameters to tweak distribution
N_vpll     =2e6;    % Number of radial positions
input.N_total=N_vpll;  % Total number of particles

input.Ekin=50e3*ones(input.N_total,1);%+(rand(input.N_total,1)-0.5)*1e3;    % Total kinetic energy [eV]
input.Z=1;                      % Mass number
input.m=const.mD;               % Mass [kg]

%% Make mesh
%  Fill in R positions and E_perp
% R_ax=dim.R0+dim.X_axis;
% R=linspace(R_ax+dim.a*0,dim.R0+dim.a*1,N_R);
psi=rand(input.N_total,1)*0.5*maps.psi_global;
psi_norm=interp1(dim.psi_scale,1:dim.NB_PSI,psi);
theta_norm=1+rand(input.N_total,1)*512;

R=dim.R0+interp2(1:size(maps.X_PR,1),1:size(maps.X_PR,2),maps.X_PR,psi_norm,theta_norm);
Z=interp2(maps.Z_PR,psi_norm,theta_norm);
phi=rand(input.N_total,1);

vpll = (rand(N_vpll,1)-0.5).*sqrt(2*input.Ekin*(const.eV/input.m));   

x=[R,Z,phi];    % Position matrix

save('../particles_equilibrium/input/resonance_sampling.mat','-v7.3','input','x','vpll');
end