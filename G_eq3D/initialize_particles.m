function [x,v,input,output,ejected]=initialize_particles
%%initialize_particles Loads in distribution and takes
%%the part this processor should take care of.
%
global const dim par qom
% Functions to convert energy [eV] to velocity [m/s] and vica verca
E_to_v=@(m,E) sqrt(2*E*(const.eV/m));
v_to_E=@(m,v) 0.5*(m/const.eV)*v.^2;

if ~par.TEST_DIST
    dist=load(par.LOADNAME);
end
%if isfield(dist,'input')
%    if ~isfield(dist.input,'x')
%        % put back particles at the right place
%        if (min(dist.x(:,1)-dim.R0)<min(dim.scale_X))
%            dist.x(:,1)=dist.x(:,1)+dim.R0;
%            disp('WARNING : Shifting initial data of -R0');
%        end
%        dist.input.x=dist.x;
%    end
%end
%% 1. Load the distribution (file) with all particles
if par.TEST_DIST
    %% 1.0. TEST DISTRIBUTION
    warning('debug particles')
    input.m=const.mD;            % Mass
    input.Z=1;                   % Charge
    qom=input.Z*const.eV/input.m;       %Kinectic constant
    
    [x,vpll,vperp]=make_test_dist(const,dim,1,input.m);
    % Input and energy
    input.N_job=size(x,1);
    input.N_total=input.N_job;
    input.particle_nr=(1:input.N_job)';
    
    Eperp=v_to_E(input.m,vperp);
    input.Ekin=v_to_E(input.m,vpll)+Eperp;
    
    % Determine velocity vector
    [N,b]=make_normal_vector(input,x,input.N_job); % Vector of perpendicular velocity
    v=bsxfun(@times,vpll,b)+bsxfun(@times,vperp,N);
    %% 1.1. Old distribution file
elseif isfield(dist,'Nalphas_simulated')
    load_from_old_dist_description;
    %% 1.2. New distribution file (with input struct)
	%% or   PRECESSION LOADED FILE
elseif isfield(dist,'input')
    %% 1.2.1. Split for this job
    input=dist.input;
    input=dist.input;
	if isfield(dist,'x')
	   input.x=dist.x;
    end
    if isfield(input,'x')
	   x=input.x;
	end
	if isfield(input,'v')
	   v=input.v;
    end
    if isfield(dist,'vpll')
	   input.vpll=dist.vpll;
    end
    if isfield(dist,'rare')
        expr_job=false(input.N_total,1);
        expr_job(dist.rare)=true;
        input.rare=dist.rare;
    else
        expr_job=get_expr_job([input.N_total, 1],par.NB_PROCESS,par.PROCESS_NUMBER);  % Get expression which particles this job will simulate
    end
    input.particle_nr=find(expr_job);
    input.N_job=sum(expr_job(:)); % Number of particles in simulation
	
	%	%% Velocity
   if ~isfield(input,'v') && isfield(input,'vpll')
        %Find perpendicular velocity
        vpll=input.vpll(expr_job);
        input.vpll_ini=vpll;
        input=rmfield(input,'vpll');
        vtot=E_to_v(input.m,input.Ekin);
		vtot=vtot(expr_job);
		x=x(expr_job,:);
        vperp=sqrt(vtot.^2-vpll.^2);
        % 1.2.4.2.2. Determine v by direction of b
        [N,b]=make_normal_vector(input,x,input.N_job); % Vector of perpendicular velocity
        v=bsxfun(@times,vpll,b)+bsxfun(@times,vperp,N);
		input.v=v;
    elseif ~isfield(input,'v')
        error('Distribution file has improper velocity data')
    end
	if ~isfield(input,'N_total'); input.N_total=size(dist.x,1); end

    
    %% 1.2.2. Load general parameters
    input=reduce_struct(input,input.N_total,expr_job);
    qom=input.Z*const.eV/input.m;       %Kinectic constant
    
    % Any precession information
    if par.GET_PREC_INFO && isfield(dist,'prec')
        dist.prec=remove_fields(dist.prec,{'ind'});
        input.prec=reduce_struct(dist.prec,input.N_total,expr_job);
    end
    
	par.APPLY_FC=0;
	if isfield(input,'ang_rot')
		par.APPLY_FC=1;
		input.Fc_field=input.ang_rot.^2/qom;
	else
		input.Fc_field=x(:,1)*0;
    end
	
    %% 1.2.3. Location
    if isfield(input,'x') && isfield(input,'v')
        warning('Using position and velocity from prec distribution file')
        [x,v,input,output,ejected]=make_start_arrays(input);
	else
		x=dist.x(expr_job,:);
    end    

	%% 1.3. New distribution file (without input struct)
	% not implemented properly yet
%elseif isfield(dist,'x') 
%	x=dist.x;
%	if isfield(dist,'v')
%	   v=dist.v;
%	end
%	%% Velocity
%   if ~isfield(dist,'v') && isfield(dist,'vpll')
%        %1.2.4.2.1 Find perpendicular velocity
%        vpll=dist.vpll(expr_job);
%        input.vpll_ini=vpll;
%        dist=rmfield(dist,'vpll');
%        vtot=E_to_v(input.m,input.Ekin);
%        vperp=sqrt(vtot.^2-vpll.^2);
%        % 1.2.4.2.2. Determine v by direction of b
%        [N,b]=make_normal_vector(input,x,input.N_job); % Vector of perpendicular velocity
%        v=bsxfun(@times,vpll,b)+bsxfun(@times,vperp,N);
%    else
%        error('Distribution file has improper velocity data')
%    end
%    obsolete : %% 1.3. PRECESSION LOADED FILE
%elseif  par.GET_PREC_INFO && isfield(dist,'prec') 
%    if isfield(dist,'input')
%		% Possible rare value added for specific study
%		% Save the rare value in the input struct and use only those
%		input=dist.input;
%		if isfield(dist,'rare')
%			input.rare=dist.rare;   % Save which particles are in this simulation by storing rare inside input
%			% Determine the number of particles in this simulation
%			expr_job=false(size(input.x,1),1);
%			expr_job(input.rare)=true;
%			if ~isfield(input,'particle_nr'); input.particle_nr=find(expr_job); end
%			input.N_job=sum(expr_job(:));
%			input=reduce_struct(input,input.N_total,expr_job);
%		else 
%			expr_job=true(size(input.x,1),1);
%			if input.N_job~=size(input.x,1); error('Size of N_job and x-array do not match'); end
%		end
%		qom=input.Z*const.eV/input.m;       %Kinectic constant
%		
%		% Load prec-data as an `input' variable.
%		input.prec=remove_fields(dist.prec,{'ind'}); % Remove ind, since information is already stored in pop (if combined). This to avoid confusing with particle_nrs and such
%		disp('PREC-data loaded from PREC-file')
%		
%		if isfield(dist,'par')
%			if dist.par.NB_PROCESS~=par.NB_PROCESS || dist.par.PROCESS_NUMBER ~= par.PROCESS_NUMBER
%				error('`Boxes'' not equally devided in PREC-file');
%			elseif dist.par.DISTNAME ~= par.DISTNAME
%				error('Distribution file not used consistently')
%			end
%		end
%	else
%		load_from_old_dist_description;
%	end

else
    error('PREC/DISTRIBUTION file does not contain proper fields for identification')
end
if isfield(dist,'x')
    % remove useless field
    dist=rmfield(dist,'x');
end
%% 2. Shift initial position so input x is in gyrocenter
if ~isfield(input,'x_gc')
    input.x_gc=input.x; % Store gc input
    [~,B]=B_interpolation(input.x_gc,'2D'); %equilibrium magnetic field
    Bfield=sqrt(dot(B,B,2));
    input.x=input.x_gc-bsxfun(@times,(1/qom)./Bfield,cross(v,b,2));
end
x_gc=input.x_gc;
x=input.x;

% Ejected particles
if exist('dist','var') && isfield(dist,'ejected')
    ejected=dist.ejected(expr_job);
else
    ejected=false(input.N_job,1);
end

%% 3. Find physical quantities / store input

% in case the toroidal rotation is part of the input data
par.APPLY_FC=0;
% Fc_field=x(:,1)*0;
if isfield(input,'ang_rot')
    par.APPLY_FC=1;
    disp('Simulation including centrifugal effects') 
end
if exist('output')~=1
	[x,v,input,output,ejected]=make_start_arrays(input,ejected,qom);
end
% keeping track of the particles initial gc positions and velocities
% (the input.x_gc and input.v are replace at the end of a simulation so you
%  that can run several simulations consecutively)
input.x_ini=input.x_gc;
input.v_ini=input.v;
[~,B]=B_interpolation(x);
Bfield_sq=dot(B,B,2);
Bfield=sqrt(Bfield_sq);
vpll=dot(B,v,2)./Bfield;
input.vpll_ini=vpll;
return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A NESTED FUNCTION FOR OUTPUT AND INPUT ARRAY
function [x,v,input,output,ejected]=make_start_arrays(input,ejected,qom)
global par dim maps const
v_to_E=@(m,v) 0.5*(m/const.eV)*v.^2;
if nargin==1
    ejected=false(input.N_job,1);
end
% Start the x and v arrays
x=input.x;
v=input.v;

% Pre-allocate output
output.x=NaN([size(x),par.NB_TIME_STAMPS]);
output.v=NaN([size(v),par.NB_TIME_STAMPS]);

output.x_ej=NaN(size(x));
output.vpll_ej=NaN(input.N_job,1);
output.mm_ej=NaN(input.N_job,1);

psi=interp2(dim.scale_Z,dim.scale_X,maps(1).psi_XZ,x(:,2),x(:,1)-dim.R0,'*linear');% find psi-value

% condition >=0 valid for JET 
% for ASDEX, use psi<=0 instead
if mean(dim.psi_scale_correct)<0
    ejected=ejected | isnan(psi) | psi>=0;
else
    ejected=ejected | isnan(psi) | psi<=0;
end

x(ejected,:)=NaN;
v(ejected,:)=NaN;

output.time_step_loss=NaN(input.N_job,1);
output.time_step_loss(ejected)=0;

output.nr_midplane_crossing=zeros(input.N_job,1);
output.nr_vpll_crossing    =zeros(input.N_job,1);

% Magnetic moment
[~,B]=B_interpolation(x); %magnetic field in 3D field
Bfield=sqrt(dot(B,B,2));

vperp=sqrt(dot(v,v,2)-(dot(v,B,2)./Bfield).^2);
Eperp=v_to_E(input.m,vperp);

input.mm=Eperp./Bfield; % Magnetic moment


% pphi_kin (also used to determining pphi_an)
dt=par.dt;
par.dt=0.5*dt;

[~,v_plus_1]=time_step_integration_GT_eq_struct(x,v,input.Fc_field);
par.dt=dt;
input.pphi_kin= get_pphi_kin(input,x,v_plus_1,psi)/(input.Z*maps(1).psi_global);

% extra field to record distance covered by particles(mostly for Poincare plot)
output.Delta_l=zeros(input.N_job,1);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [N,b,gyro_angle]=make_normal_vector(input,x,N_simulated)
%% Find normal vector b
% Determine direction B with inclusion of RMP-field
if size(x,1)<2; error('Code errors likely for less than 2 particles per job. Make sure there are more particles than jobs'); end;
[~,B]=B_interpolation(x,'2D');
Bfield_inv=dot(B,B,2).^-0.5; % inverse magnitude
b=bsxfun(@times,Bfield_inv,B);

%% Make 2 perpendicular vectors
% Perpendicular unit vector with perpendicular to b and eZ (unit in Z-direction)
u=zeros(N_simulated,3);
u(:,1)=sqrt(1./(1+(b(:,1)./b(:,3)).^2)); % x-component
u(:,3)=-(b(:,1)./b(:,3)).*u(:,1);        % phi-component

% Other perpendicular vector by cross product
w=cross(u,b,2);

%% Make vector normal vector of perpendicular velocity
% normal vector N = cos(phase)* u + sin(phase)* w
seed=1;
rng(seed);
% gyro_angle=rand(size(x,1),1)*2*pi;
gyro_angle=rand(2000,1)*2*pi;
sub_ind=mod(input.particle_nr-1,2000)+1;
gyro_angle=gyro_angle(sub_ind);

N=   bsxfun(@times,cos(gyro_angle),u)...
    +bsxfun(@times,sin(gyro_angle),w);
end

%%
function [struct] = reduce_struct(struct,N_total,expr_job)
fnames=fieldnames(struct);
for i=1:length(fnames)
    if isstruct(struct.(fnames{i}))
        struct.(fnames{i})=reduce_struct(struct.(fnames{i}),N_total,expr_job);
    elseif size(struct.(fnames{i}),1)==N_total
        struct.(fnames{i})=struct.(fnames{i})(expr_job,:,end);
    end
end
end