function [x,v,input,output,ejected]=initialize_particles
%%initialize_particles Loads in distribution and takes
%%the part this processor should take care of.
%
global const dim par qom maps
% Functions to convert energy [eV] to velocity [m/s] and vica verca
E_to_v=@(m,E) sqrt(2*E*(const.eV/m));
v_to_E=@(m,v) 0.5*(m/const.eV)*sum(v.^2,2);

if ~par.TEST_DIST
    dist=load(par.LOADNAME);
end
if par.CONTINUE_PREV_SIM
    dist_prev=load(par.LOADNAME_PREV);    
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
%     input.v(1:2:end)=sqrt(3)*input.v(1:2:end);
%     input.v(2:2:end)=0.1*input.v(2:2:end);
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
else
    error('PREC/DISTRIBUTION file does not contain proper fields for identification')
end
if ~isfield(input,'Ekin')
    input.Ekin=v_to_E(input.m,input.v);
end
if ~isfield(input,'N_total')
    input.N_total=length(input.Ekin);
end
    
par.MARKER_WEIGHT        = par.INPUT_POWER * par.time_scale(end) / (sum(input.Ekin)*const.eV);       % representative weight of each marker (used in NDD) = (power * time) / sum(Ekin)

if isfield(dist,'rare')
    expr_job=false(input.N_total,1);
    expr_job(dist.rare)=true;
    input.rare=dist.rare;
else
    if par.SD_MARKERS_END==0
        expr_job=get_expr_job([input.N_total, 1],par.NB_PROCESS,par.PROCESS_NUMBER);  % Get expression which particles this job will simulate
    else
        expr_job_1=get_expr_job([par.SD_MARKERS_END, 1],par.NB_PROCESS,par.PROCESS_NUMBER);  % Get expression which particles this job will simulate
        expr_job_2=get_expr_job([input.N_total-par.SD_MARKERS_END, 1],par.NB_PROCESS,par.PROCESS_NUMBER);  % Get expression which particles this job will simulate
        disp('splitting in 2 parts the initial data: sd and ionization')
        nb1=length(find(expr_job_1));
        expr_job=[ expr_job_1 ;expr_job_2 ];
        %storing the number of sd markers to adjust the birth matrix!
        par.nb_sd=nb1;
    end
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
    [N,gyro_angle]=make_normal_vector_simple(input,x,input.N_job); % Vector of perpendicular velocity
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
    %         warning('Using position and velocity from prec distribution file')
    [x,v,input,output,ejected]=make_start_arrays(input);
else
    x=dist.x(expr_job,:);
end




%%  for constant birth rate throughout simulation (for slowing down study mostly)
if par.PROGRESSIVE_BIRTH  & par.NB_BIRTH_CHUNKS>0
    if input.N_job < par.NB_BIRTH_CHUNKS
        error('input.N_job < par.NB_BIRTH_CHUNKS !!!!')
    end
    
    birth_time_scale=linspace(0,par.time_scale(end),par.NB_BIRTH_CHUNKS);
    birth_matrix=zeros(par.NB_BIRTH_CHUNKS,input.N_job);
    if par.SD_MARKERS_END==0
        N_CHUNK=floor(input.N_job/par.NB_BIRTH_CHUNKS)
        NB_FLOOR_P1=input.N_job-N_CHUNK*par.NB_BIRTH_CHUNKS; % number of times when we will add a bit more
        NB_FLOOR=par.NB_BIRTH_CHUNKS-NB_FLOOR_P1;      % number of times when we add the reference
        for ind_bm=1:NB_FLOOR
            birth_matrix(ind_bm,1:ind_bm*N_CHUNK)=1;
        end
        for ind_bm=NB_FLOOR+1:par.NB_BIRTH_CHUNKS
            birth_matrix(ind_bm,1:NB_FLOOR*N_CHUNK+(ind_bm-NB_FLOOR)*(N_CHUNK+1))=1;
        end
    else
        par.nb_ionizations=input.N_job-par.nb_sd;
        birth_matrix_sd=ones(par.NB_BIRTH_CHUNKS,par.nb_sd);
        birth_matrix_ionization=zeros(par.NB_BIRTH_CHUNKS,par.nb_ionizations);
        N_CHUNK=floor(par.nb_ionizations/par.NB_BIRTH_CHUNKS)
        NB_FLOOR_P1=par.nb_ionizations-N_CHUNK*par.NB_BIRTH_CHUNKS; % number of times when we will add a bit more
        NB_FLOOR=par.NB_BIRTH_CHUNKS-NB_FLOOR_P1;      % number of times when we add the reference
        for ind_bm=1:NB_FLOOR
            birth_matrix_ionization(ind_bm,1:ind_bm*N_CHUNK)=1;
        end
        for ind_bm=NB_FLOOR+1:par.NB_BIRTH_CHUNKS
            birth_matrix_ionization(ind_bm,1:NB_FLOOR*N_CHUNK+(ind_bm-NB_FLOOR)*(N_CHUNK+1))=1;
        end
        birth_matrix=[birth_matrix_sd  birth_matrix_ionization];
    end
    par.birth_matrix = birth_matrix;
    par.birth_time_scale = birth_time_scale;
else
    par.PROGRESSIVE_BIRTH = 0;
    par.NB_BIRTH_CHUNKS = 1; 
    par.birth_matrix = ones(par.NB_BIRTH_CHUNKS,input.N_job);
    par.birth_time_scale = 0;
end

try 
    inputvalues=load(par.DISTNAME);
catch
    error('input file not found !')
end



%% simplistic relativistic correction
if par.RELATIVISTIC_CORR==1
    gamma=1./sqrt(1-(input.v./const.C0).^2);
    input.v=gamma.*input.v;
end
%%
if isfield(dist,'x')
    % remove useless field
    dist=rmfield(dist,'x');
end
%% 2. Shift initial position so input x is in gyrocenter
if ~isfield(input,'x_gc') | (par.INPUT_TRUEPOS)
    [~,B]=B_interpolation(input.x,'2D'); %equilibrium magnetic field
    if (~par.INPUT_TRUEPOS)
        input.x_gc=input.x; % Store gc input
        Bfield=sqrt(dot(B,B,2));
        Bfield_inv=dot(B,B,2).^-0.5; % inverse magnitude
        b=bsxfun(@times,Bfield_inv,B);
        input.x=input.x_gc-bsxfun(@times,(1/qom)./Bfield,cross(v,b,2));
    else
        Bfield_sq(:)=dot(B(:,:),B(:,:),2);
        input.x_gc(:,:)=input.x(:,:)+bsxfun(@times,1./(qom*Bfield_sq(:)),cross(input.v(:,:),B(:,:),2));
    end
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

% finally take care of pphi initial value (normalized by Z*psi_global)
[~,v_plus_1]=time_step_integration_GT_eq_struct(input.x,input.v,input.Fc_field);
psi=interp2(dim.scale_Z,dim.scale_X,maps(1).psi_XZ,x(:,2),x(:,1)-dim.R0,'*cubic');% find psi-value
input.pphi_kin= get_pphi_kin(input,x,v_plus_1,psi)/(input.Z*maps(1).psi_global);

%%
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

output.x_ej_prev4=NaN(size(x));
output.x_ej_prev3=NaN(size(x));
output.x_ej_prev2=NaN(size(x));
output.x_ej_prev=NaN(size(x));
output.x_ej=NaN(size(x));
output.x_ej_next=NaN([size(x),par.RECORD_AFTER_LOSS]);
output.flag_loss_counter=zeros(input.N_job,1);
output.x_CXn=NaN(size(x));
output.x_CXi=NaN(size(x));
output.Ekin_ej=NaN(input.N_job,1);
output.vpll_ej=NaN(input.N_job,1);
output.mm_ej=NaN(input.N_job,1);
output.ejected_wall=NaN(input.N_job,1);
if par.CALCULATE_CX
    output.CX_rate=zeros(input.N_job,par.NB_TIME_STAMPS); 
    output.CX_cum=zeros(input.N_job,par.NB_TIME_STAMPS);
    output.CX_flag=zeros(input.N_job,par.NB_TIME_STAMPS);
    % output.ejected_CX=zeros(input.N_job,par.NB_TIME_STAMPS);
    output.ejected_CX=NaN(input.N_job,1);
end
if par.CALCULATE_NDD
    output.ndd_val=zeros(input.N_job,par.NB_TIME_STAMPS);
    output.ndd_cum=zeros(input.N_job,par.NB_TIME_STAMPS);
end
if par.CALCULATE_PDEP
    output.Pdep_e=zeros(input.N_job,par.NB_TIME_STAMPS);
    output.Pdep_i=zeros(input.N_job,par.NB_TIME_STAMPS);
    output.Etot_e=zeros(input.N_job,par.NB_TIME_STAMPS);
    output.Etot_i=zeros(input.N_job,par.NB_TIME_STAMPS);
    output.delta_Ekin=zeros(input.N_job,par.NB_TIME_STAMPS);
    output.cum_tor_ang_mom=zeros(input.N_job,par.NB_TIME_STAMPS);
end
if par.COULOMB_COLL
    % thermal part of deposited energy should be
    % used to evaluate additional power during steady-state
    output.Edep_th=zeros(input.N_job,1);
end

psi=interp2(dim.scale_Z,dim.scale_X,maps(1).psi_XZ,x(:,2),x(:,1)-dim.R0,'*linear');% find psi-value

% condition >=0 valid for JET 
% for ASDEX, use psi<=0 instead
if ~par.USE_VESSEL_LIMIT 
    if mean(dim.psi_scale_correct)<0
        ejected=ejected | isnan(psi) | psi>=par.PSI_LIMIT_EJECTED;
    else
        ejected=ejected | isnan(psi) | psi<=par.PSI_LIMIT_EJECTED;
    end
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
if par.calculate_length_trajectory
    output.Delta_l=zeros(input.N_job,1);
end
if par.COULOMB_COLL
    if par.CALCULATE_DEVIATION
		output.deviation=zeros(input.N_job,par.NB_TIME_STAMPS);
	end
	output.ejected_sd=zeros(input.N_job,par.NB_TIME_STAMPS);
end
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
