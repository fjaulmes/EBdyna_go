function build_thermal_helium_ST
%build_phase_space_sampling Builds evenly spaced sampling in 3D phase space
%   Loads in plasma parameters and creates particles which start on
%   Z=Z_axis of equilibrium.

if verLessThan('matlab','8.4')
    intp_method='cubic';
else
    intp_method='phchip';
end


seed=1;
% seed='shuffle';
rng(seed);
%% Load default parameters
par.paths=initialize_folder_names_struct;   % Folder names
[const,maps,dim]=load_distr_maps(par);
dim2=load([par.paths.DATA_FOLDER,'q_profile.mat'],'q_initial_profile'); dim=combine_structs(dim,dim2);
dim2=load(strcat(par.paths.DATA_FOLDER,'psi_profiles.mat'),'psi_pol_initial_profile'); dim=combine_structs(dim,dim2);

%% ST-coordinates
psi_0=max(dim.psi_pol_initial_profile)-dim.psi_pol_initial_profile;  % Define psi as 0 zero in center of plasma and finite at edge
chi=cumtrapz(psi_0,dim.q_initial_profile);                           % toroidal flux (psi) = int_0^psi  q d(psi) with trapezium integration
psi_star=psi_0-chi;               

r  =   sqrt(2*chi);  % Definition of r-coordinate

% r at q=1 (with top of psi_star)
ind_r0      = interp1(dim.q_initial_profile,1:length(dim.q_initial_profile),1,intp_method);     % Index q=1
h_ind_r0_int= floor(ind_r0/2);                                                              % Integer index between center and q=1, for finding rmix
r0          = interp1(dim.q_initial_profile,r,1,intp_method);                            % r-value at q=1

% Mixing radius (psi_star = 0). Interpoldation from h_ind_r0_int since we do not want the center point
ind_rmix = interp1(psi_star(h_ind_r0_int:end), (h_ind_r0_int:length(dim.psi_pol_initial_profile)	),0,intp_method);     % index mixing radius
rmix     = interp1(psi_star(h_ind_r0_int:end),r(h_ind_r0_int:end                                    ),0,intp_method);     % r-valua mixing radius

%% Parameters to tweak distribution
N_R     =1000;  % Number of radial positions
N_Eperp =200;   % Number of different pitches
input.N_total=N_R*N_Eperp;  % Total number of particles

input.Ekin=2.8e3*ones(input.N_total,1);%+(rand(input.N_total,1)-0.5)*1e3;    % Total kinetic energy [eV]
input.Z=2;                      % Mass number
input.m=const.mHe;               % Mass [kg]

%% Make mesh
%  Fill in R positions and E_perp
r_max=rmix*1.5;
r_pos=rand(input.N_total,1)*r_max;
theta_ind=rand(input.N_total,1)*size(maps.X_PR,1);

psi_norm=interp1(r,1:dim.NB_PSI,r_pos);
R=dim.R0+interp2(maps.X_PR,psi_norm,theta_ind,'*cubic');
Z=interp2(maps.Z_PR,psi_norm,theta_ind,'*cubic');

vtot=sqrt(input.Ekin*2*const.eV/input.m);
prll_chance=1-2*rand(input.N_total,1);
vpll=vtot.*prll_chance;

x=[R,Z,rand(input.N_total,1)];    % Position matrix

save(['../particles_equilibrium/input/',datestr(now,'yyyy-mm-dd'),'_thermal_helium_ST.mat'],'-v7.3','input','x','vpll');
end