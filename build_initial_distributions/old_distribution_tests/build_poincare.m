function build_poincare
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
N_R     =2^10;    % Number of radial positions
N_Eperp =1;   % Number of different pitches
N_phi   =1;
input.N_total=N_R*N_Eperp*N_phi;  % Total number of particles

Ekin=0.04*ones(input.N_total,1);%+(rand(input.N_total,1)-0.5)*1e3;    % Total kinetic energy [eV]
input.Z=1;                      % Mass number
input.m=const.mD;               % Mass [kg]

%% Make mesh
%  Fill in R positions and E_perp
% R_ax=dim.R0+dim.X_axis;
% R=linspace(R_ax+dim.a*0,dim.R0+dim.a*1,N_R);
% rho=linspace(0.7,1,N_R)';
% psi_overline=rho.^2;
psi_overline=linspace(0.5,1,N_R/2)';
psi=(1-psi_overline)*maps.psi_global;
psi_norm=interp1(dim.psi_scale,1:dim.NB_PSI,psi);
q=interp1(dim.psi_scale,dim.q_initial_profile,psi);
theta=2*pi*rand(input.N_total/2,1);
phi=q.*theta;
psi_norm=[psi_norm;psi_norm];
theta=mod([theta;theta],2*pi);
phi=mod([phi;phi+pi],2*pi); % Replicate the positions on the other side of the tokamak to fill all islands (i.e. those created with an n=2 RMP)
theta_norm=interp1(linspace(0,2*pi,size(maps.X_PR,1)),1:size(maps.X_PR,1),theta);
R=dim.R0+interp2(1:size(maps.X_PR,1),1:size(maps.X_PR,2),maps.X_PR,psi_norm,theta_norm);

% Z=zeros(input.N_total,1)+dim.Z_axis;
Z=interp2(maps.Z_PR,psi_norm,theta_norm);

E_perp=0;
%phi=-0.05;

% sign_vpll=[-1,1];
sign_vpll=1;

% Order in which particles are sorted as vector
order=randperm(input.N_total)'; % Random order
% order=(1:input.N_total)'; % Sequential order

[R,Eperp,sign_vpll]=ndgrid(R,E_perp,sign_vpll); % Make mesh (Ekin,R,Eperp,phi)

% Sort them in order
R=R(order);
Z=Z(order);
input.Ekin=Ekin(order);
Eperp=Eperp(order).*input.Ekin;
phi=phi(order);
sign_vpll=sign_vpll(order);

% Parallel velocity as input variable. Gyro phase in start of simulation
Epll=input.Ekin-Eperp;
vpll = sign_vpll.*sqrt(2*Epll *(const.eV/input.m));   

x=[R,Z,phi];    % Position matrix

save('../particles_equilibrium/input/poincare_distribution_0.04eV_rand_edge.mat','-v7.3','input','x','vpll');
end