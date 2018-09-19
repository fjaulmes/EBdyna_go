clear all;
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
N_R     =2^7;    % Number of radial positions
input.N_total=N_R;  % Total number of particles

Ekin=(100)*ones(input.N_total,1);%+(rand(input.N_total,1)-0.5)*1e3;    % Total kinetic energy [eV]
input.Z=1;                      % Charge number
input.m=const.mD;               % Mass [kg]

%% Make mesh
%  Fill in R positions and E_perp
% R_ax=dim.R0+dim.X_axis;
% R=linspace(R_ax+dim.a*0,dim.R0+dim.a*1,N_R);
% rho=linspace(0.7,1,N_R)';
% psi_overline=rho.^2;
SIGN_PSI_AXIS=sign(dim.psi_scale(1));
psi_overline=linspace(0.6,0.98,N_R/2)';
psi=SIGN_PSI_AXIS*(1-psi_overline)*maps.psi_global;
psi_norm=interp1(dim.psi_scale,1:dim.NB_PSI,psi);
q=interp1(dim.psi_scale,dim.q_initial_profile,psi);
theta=2*pi*rand(input.N_total/2,1);
phi=q.*theta;
psi_norm=[psi_norm;psi_norm];
theta=mod([theta;theta],2*pi);
phi=mod([phi;phi+pi],2*pi); % Replicate the positions on the other side of the tokamak to fill all islands (i.e. those created with an n=2 RMP)
theta_norm=interp1(linspace(0,2*pi,size(maps.X_PR,1)),1:size(maps.X_PR,1),theta);

R=dim.R0+interp2(1:size(maps.X_PR,1),1:size(maps.X_PR,2),maps.X_PR,psi_norm,theta_norm);
Z=interp2(maps.Z_PR,psi_norm,theta_norm);


%% Velocity

E_perp=0.2; % fraction in perp energy
sign_vpll=1;

% Order in which particles are sorted as vector
order=randperm(input.N_total)'; % Random order

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

%% save distribution file

x=[R,Z,phi];    % Position matrix
save('./pphi_test_distribution_100eV_edge.mat','-v7.3','input','x','vpll');





