function build_phase_space_sampling
%build_phase_space_sampling Builds evenly spaced sampling in 3D phase space
%   Loads in plasma parameters and creates particles which start on
%   Z=Z_axis of equilibrium.

seed=1;
% seed='shuffle';
rng(seed);
%% Load default parameters
par.paths=initialize_folder_names_struct;   % Folder names
[const,maps,dim]=load_distr_maps(par);
dim2=load([par.paths.DATA_FOLDER,'q_profile.mat']...
    ,'q_initial_profile');
dim=combine_structs(dim,dim2);
%% Parameters to tweak distribution
N_R     =1000;    % Number of radial positions
N_Eperp =200;   % Number of different pitches
N_phi   =1;
input.N_total=N_R*N_Eperp*N_phi;  % Total number of particles

Ekin=35e3*ones(input.N_total,1);%+(rand(input.N_total,1)-0.5)*1e3;    % Total kinetic energy [eV]
input.Z=1;                      % Mass number
input.m=const.mD;               % Mass [kg]

%% Make mesh
%  Fill in R positions and E_perp
% R_ax=dim.R0+dim.X_axis;
% R=linspace(R_ax+dim.a*0,dim.R0+dim.a*1,N_R);
q=linspace(1.24,1.8,N_R)';
psi_norm=interp1(dim.q_initial_profile,1:dim.NB_PSI,q);
R=dim.R0+interp2(1:dim.NB_PSI,1:dim.NB_PSI,maps.X_PR,psi_norm,1,'*cubic');

E_perp=linspace(0,1,N_Eperp);
phi=2*pi/49*33;
sign_vpll=[1];

% Order in which particles are sorted as vector
% order=randperm(input.N_total)'; % Random order
order=(1:input.N_total)'; % Sequential order

[R,Eperp,phi,sign_vpll]=ndgrid(R,E_perp,phi,sign_vpll); % Make mesh (Ekin,R,Eperp,phi)

% Sort them in order
R=R(order);
input.Ekin=Ekin(order);
Eperp=Eperp(order).*input.Ekin;
phi=phi(order);
sign_vpll=sign_vpll(order);

% Parallel velocity as input variable. Gyro phase in start of simulation
Epll=input.Ekin-Eperp;
vpll = sign_vpll.*sqrt(2*Epll *(const.eV/input.m));    % Co-field
% vpll =-sqrt(2*Epll *(const.eV/input.m));  % Counter-field

% Other coordinates
Z=zeros(input.N_total,1)+dim.Z_axis;

x=[R,Z,phi];    % Position matrix

save('../particles_equilibrium/input/initial_distr_phase_space4.mat','-v7.3','input','x','vpll');
end