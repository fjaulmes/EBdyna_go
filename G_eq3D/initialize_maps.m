function [m,d] = initialize_maps(maker,load_worker)
global const maps dim map_data par
%initialize_sim_data Loads in data files for simulation
%   Returns structs with constants, B-maps and dimension parameters
if isempty(const)
    filename=[par.paths.DATA_FOLDER,'physics_constants.mat'];
    const=load(filename);
end

switch nargin
    case 1
        load_worker=false;
end

%% Make the maps and store in shared memory
if load_worker
    % Broadcast by making a file that the process is busy making shared RAM
    busy=1;  %#ok<NASGU>
    save('busy.txt','busy','-ascii');
    
    % Define keys
    if ispc
        map_data.keys.key_maps=[par.ID_NAME,'_key_maps'];
    elseif isunix    
        map_data.keys.key_maps=54321;
    end
    
    [m,d] = initialize_maps(true,false);  % Load maps and dimensions conventionally
    map_data.dim=d;
    
    % COPY IN SHARED RAM
    % Put the matrix with largest size (nr dims) last
    mnames=fieldnames(m); p=zeros(size(mnames));
    for i=1:length(mnames)
        p(i)=ndims(m.(mnames{i}));
    end
    [~,I]=sort(p);
    m=orderfields(m,I);
    
    % Make a matrix to be compatible with sharedmatrix function (UNIX)
    m(2).n_maps=2; m(1).n_maps=2;     % Put an extra (empty) array in the maps-struct, since sharedmatrix needs at least 2 elements;
	m=remove_fields(m,{'n_maps'});
    
    if ispc
        SharedMemory('clone',map_data.keys.key_maps  ,m);
    elseif isunix    
        sharedmatrix('clone',map_data.keys.key_maps  ,m);
    end
    
    % Store keys in file
    save('shared_memory_keys.mat','-struct','map_data'); % Save the keys
    pause(1)
    % Broadcast that the file is ready to be read (or remove the broadcasted busy file)
    delete('busy.txt');
    disp('Shared RAM made succesfully!');
    return
end
%% Load maps from memory
if ~maker
    if nargout~=0
        error('Cannot return attached variables (Shared Memory Linked). These can only be loaded by globalizing parameters')
    end
    if exist('shared_memory_keys.mat','file')   % Load from shared RAM
	    disp('Loading from shared RAM...');
        map_data=load('shared_memory_keys.mat');
        % Attach shared RAM
        if ispc
            maps  = SharedMemory('attach',map_data.keys.key_maps);
        elseif isunix    
            maps  = sharedmatrix('attach',map_data.keys.key_maps);
        end
        dim=map_data.dim;
    else
        warning('No shared memory keys found!')
        % Made by a process instead of in advance!
        rng('shuffle');
        pause(mod(par.PROCESS_NUMBER,16)*10+rand*5)  % Make sure the processes do not make the file precisely at the same time
        
        % Check if anyone is making the shared RAM / has made the shared RAM
        if exist('busy.txt','file')
            busy=load('busy.txt');
            if busy % wait and try again
                pause(320)
				busy=load('busy.txt');
				if busy
					error('SOME PROGRESS IS TAKING TOO LONG TO MAKE THE SHARED RAM (OR FAILED WITH BUSY-FILE ON)')
				end
				initialize_maps(false,false);
				return
			end
            error('Busy-file indicates shared RAM is made, but no keys are found')
        else % Make shared RAM yourself
            busy=1; %#ok<NASGU>
            save('busy.txt','busy','-ascii');
            warning(['MAKING SHARED DATA BY: PROCESS ',num2str(par.PROCESS_NUMBER)])
            
            initialize_maps(true,true);  % Make the maps
            
            % Load memory as shared memory
            maps=[];
            initialize_maps(false,false); % Load them as shared memory
            return
        end
    end
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                   - ONLY FOR A 'MAKER'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Load in map files....')
%% Load in map files
filename=[par.paths.DATA_FOLDER,'XZsmall_fields_tokamak_pre_collapse.mat'];
m=load(filename,'Bphi_XZsmall_map','BpolX_initial_XZsmall_map','BpolZ_initial_XZsmall_map','psi_global','psi_XZsmall_map','psi_norm_XZsmall_map','theta_XZsmall_map');

d2=load(filename,'size_X','size_Z');
filename=[par.paths.DATA_FOLDER,'motions_map_dimensions.mat'];
d=load(filename,'mid_X','mid_Z','mid_Zzero','DX','scale_X','scale_Z','NB_PSI','R0','Z_axis','mid_Xzero','X_axis','psi_scale');
d=combine_structs(d,d2);
d2=load(strcat(par.paths.DATA_FOLDER,'psi_profiles.mat'),'psi_pol_initial_profile'); 
d.psi_scale_correct=d2.psi_pol_initial_profile;

if par.mode==5 || par.APPLY_SAWTOOTH
    filename=[par.paths.DATA_FOLDER,'q_profile.mat'];
    d2=load(filename,'q_initial_profile'); d=combine_structs(d,d2);
end

% Remove map name in fields, which I find unnessecary
fnames=fieldnames(m)';
for i=1:length(fnames)
    if strendswith(fnames{i},'small_map')
        m.(fnames{i}(1:end-9))=m.(fnames{i});
        m=rmfield(m,fnames{i});
    elseif strendswith(fnames{i},'_map')
        m.(fnames{i}(1:end-4))=m.(fnames{i});
        m=rmfield(m,fnames{i});
    end
end
fnames

% Check direction of toroidal field 
if sign(mean(m.Bphi_XZ(:)))<=0
    disp('Toroidal field in negative direction');
else
	disp('Toroidal field in positive direction');
end

%% Alter the 2D maps

% Truncate the 2D map
% not necessary and dangerous : removed
original_size=[d.size_X,d.size_Z];
% d.size_X=1100;
% d.size_Z=1800;
% if ~isequal([d.size_X,d.size_Z],original_size); warning('Truncating 2D map sizes'); end

% Change size of arrays
d.scale_X=d.scale_X(1:par.step.R:d.size_X);
d.scale_Z=d.scale_Z(1:par.step.Z:d.size_Z);

% All field maps (all except )
fnames=fieldnames(m)';
for fname=fnames
    if ~isempty(strfind(fname{1},'_XZ')) && isequal(size(m.(fname{1})),original_size)
        m.(fname{1})=m.(fname{1})(1:par.step.R:d.size_X,1:par.step.Z:d.size_Z);
    end
end

%  Change the middle values of the grid (size)
d.mid_X         =(d.mid_X    -1)/par.step.R+1;  % Center of plasma
d.mid_Xzero     =(d.mid_Xzero-1)/par.step.R+1;  % Center of tokamak (R=R0+mid_Xzero)
d.mid_Z         =(d.mid_Z    -1)/par.step.Z+1;  % Center of plasma
d.mid_Zzero     =(d.mid_Zzero-1)/par.step.Z+1;  % Center of plasma

% Change the grid size parameter
d.DX=d.DX*par.step.R;
if isfield(d,'DZ')
    d.DZ=d.DZ*par.step.Z;
else
    d.DZ=d.DX*par.step.Z/par.step.R;
end

% Split two theta maps, to prevent interpolation error at 2*pi border
[m]=split_theta_map(m);

%% Determine inverse of size DX and DZ, for quicker / easier interpolation
d.DX_inv=1/d.DX;
d.DZ_inv=1/d.DZ;   %NOTE: Z in 2D maps has equal distance to DX

% Replace values outside LCFS with NaN and store in m.B_2D-struct
% expr=m.psi_norm_XZ>=d.NB_PSI;
m.B_2D.BR       =m.BpolX_initial_XZ;
m.B_2D.BZ       =m.BpolZ_initial_XZ;
m.B_2D.Bphi     =m.Bphi_XZ;
% m.B_2D.BR(expr)=NaN; m.B_2D.BZ(expr)=NaN; m.B_2D.Bphi(expr)=NaN;
m=remove_fields(m,{'BpolX_initial_XZ','BpolZ_initial_XZ','Bphi_XZ'});



%% Add RMP-field to data
if par.APPLY_RMP || par.APPLY_TFR
    if par.interp_scheme==0 && strcmp(par.coord_syst,'flux')
        error('Interpolation scheme not compatible with flux based 3D grid')
    end
    
    % Load in the 3D field
    [m,d]=load_3D_maps(m,d);
    disp('3D fields loaded....')
	
    % Superimpose the 2D and 3D fields
    if par.interp_scheme~=0 && strcmp(par.coord_syst,'flux')
        % For interpolation we need to couple the psi_norm / index of a 2D map
        % to the 3D map. Therefore we need the ratio of these dimensions map
        % e.g. 513 in 2D-map and 257 in 3D-map -> ratio of 2
        
        % Ratio of psi values used in 3D grid vs. 2D grid
        ratio_2D_3D=(d.NB_PSI-1)/(d.n3D.size_3D(2)-1);
        
        % Function which gets psi index of 3D map by index 2D map
        d.n3D.ind_2D_to_3D=@(ind) 1+(ind-1)/ratio_2D_3D;
        d.n3D.ind_3D_to_2D=@(ind) 1+(ind-1)*ratio_2D_3D;
        
        % For theta and phi we know the boundaries are at 0 and 2*pi (so the
        % last coincide for proper interpolation). Therefore the scale is:
        theta_scale =linspace(0,2*pi,d.n3D.size_3D(1));
        
        % Define incrediments
        d.Dpsi  =1;
        d.Dtheta=theta_scale(2)  -theta_scale(1);
        
        % Inverse incrediments (to multiply instead of divide)
        d.Dtheta_inv=1/d.Dtheta;
        
        m.psi_norm_3D_XZ=d.n3D.ind_2D_to_3D(m(1).psi_norm_XZ); % make a new map with the index of the 3D field
    end
end

%% Correct for the size in 2D map size, for interpolation purposes
% Note that if not all steps (in par.step) are 1, then d.size_X needs to be reducted as:
size_2D=size(m.psi_XZ);
d.size_X=size_2D(1);
d.size_Z=size_2D(2);

%% Re-order B-field for the proper interpolation
% Interpolation scheme's
% 0     -   self-written linear, with vector potential
% 1     -   self-written linear         [individual B-components]
% 2     -   interp linear               [individual B-components]
% 3     -   ba_interp linear            [combined   B-components]
% 4     -   griddedInterpolant linear   [individual B-components]
% 5     -   interp cubic                [individual B-components]
% 6     -   ba_interp cubic             [combined   B-components]
% 7     -   griddedInterpolant cubic    [individual B-components]
% 8     -   interp2 spline              [individual B-components]
% 9     -   griddedInterpolant spline   [individual B-components]

switch par.interp_scheme
    case 0
        % We need to make the 2D maps for this interpolation
        m=make_potential_interpolation_maps_2D(d,m);
        return
    case {3,6}
        m.B_2D=cat(3,m.B_2D.BR,m.B_2D.BZ,m.B_2D.Bphi);
        if par.APPLY_RMP || par.APPLY_TFR
            m.B_3D=cat(4,m.n3D.BR,m.n3D.BZ,m.n3D.Bphi);
			m.n3D=remove_fields(m.n3D,{'BR','BZ','Bphi'}); % Remove B-field from the 3D-struct (to save RAM)
        end
end

%% Add sawtooth parameters (multiple profiles / functions in time)
if par.APPLY_SAWTOOTH
    [d,m]=sawtooth_collapse(d,m);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% - FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function to load files and their symmetry
function [load_data,symm]=load_3D_file(file_name,fields)
load_data=load(file_name,fields{:});

if ~all(isfield(load_data,fields))
    error('Fields missing in requested data')
end

symm=load(file_name,'symm');
if ~isfield(symm,'symm')
    symm=1;
else
    symm=symm.symm;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [m,d]=load_3D_maps(m,d)
global par

% Define fields originally needed for calculations
if par.interp_scheme==0
    BA_fields={'AR','AZ','Aphi'};
    load_fields={'Beta_r','Beta_z','Beta_phi','Beta_phi_r','Beta_phi_z','Beta_r_z','Beta_z_r','Beta_z_phi','Beta_r_phi','Beta_phi_r_z','Gamma_r','Gamma_z','Gamma_phi','Alpha_r','Alpha_z','Aphi'};
	if ~par.superimpose_2D_3D; error('Superposition of 2D and 3D required for interpolation scheme 0');	end
else
    BA_fields={'BR','BZ','Bphi','Aphi'}; % B-components and Aphi for determining p_phi in the end
    load_fields=BA_fields;
end
if par.CALCULATE_PPHI_3D
    BA_fields=cat(2,BA_fields,{'dAR_dphi','dAZ_dphi','dAphi_dphi'});
    load_fields=cat(2,load_fields,{'dAR_dphi','dAZ_dphi','dAphi_dphi'});
end

%% Which file names
if      par.APPLY_RMP && ~par.APPLY_TFR
    load_name=par.RMP_file;
elseif ~par.APPLY_RMP &&  par.APPLY_TFR
    load_name=par.TFR_file;
elseif  par.APPLY_RMP &&  par.APPLY_TFR
    % find the name of a combined file
    sl=strfind(par.TFR_file,'/');
    if ~isempty(sl)
        RMP_TFR_file=strcat(par.RMP_file(1:end-4),'_',par.TFR_file(sl(end)+1:end));
    else
        RMP_TFR_file=strcat(par.RMP_file(1:end-4),'_',par.TFR_file(1:end));
    end
    load_name=RMP_TFR_file;
else
    error('Switch error, not clear if RMP / TFR apply')
end
load_name_intp0=get_intp0_filename(load_name);

%% Try loading files

% Delete inp0 if not of the proper step size
if par.interp_scheme==0 && exist(load_name_intp0,'file')
    % Check the step size
    step2=load(load_name_intp0,'step'); step2=step2.step;
    if ~isequal(par.step,step2)
        delete(load_name_intp0)
    end
end
% 1. intp0-file
if par.interp_scheme==0 && exist(load_name_intp0,'file')
    % LOAD FILE
    [m.n3D,d.n3D.symm]=load_3D_file(load_name_intp0,load_fields);
    
    % Dimensions
    d.n3D.size_3D=size(m.n3D.(load_fields{1})); % size 3D-map
    d.n3D.size_3D=1+d.n3D.size_3D;
    phi_scale   =linspace(0,2*pi/d.n3D.symm,d.n3D.size_3D(3));
    d.Dphi  =phi_scale(2)    -phi_scale(1);
    d.Dphi_inv=1/d.Dphi;
    return
    % 2. make intp0-file
elseif par.interp_scheme==0 && ~exist(load_name_intp0,'file')
    if par.APPLY_RMP && par.APPLY_TFR && ~exist(RMP_TFR_file,'file')
        [m    ,d.n3D.symm]=make_RMP_TFR_file(m,RMP_TFR_file);
    else
        [m.n3D,d.n3D.symm]=load_3D_file(load_name,BA_fields);
    end
    
    % Make coarser 3D grid NOTE DO NOT SUPERIMPOSE 2D FIELD
    if strcmp(par.coord_syst,'toroidal')
        for BA=BA_fields
            m.n3D.(BA{1})=m.n3D.(BA{1})(1:par.step.R:d.size_X,1:par.step.Z:d.size_Z,1:par.step.phi:end);
        end
    end
    
    % Get dimensions
    d.n3D.size_3D=size(m.n3D.(BA_fields{1})); % size 3D-map
    phi_scale   =linspace(0,2*pi/d.n3D.symm,d.n3D.size_3D(3));
    d.Dphi  =phi_scale(2)    -phi_scale(1);
    d.Dphi_inv=1/d.Dphi;
    
    % Make intp0-maps
    m=make_potential_interpolation_maps_3D(d,m);
    
    % Remove obsolete fields
    m.n3D=remove_fields(m.n3D,[],load_fields);
    
    % Save intp0-file
    field_to_save=m.n3D;
    field_to_save.symm=d.n3D.symm;
    field_to_save.step=par.step;
    % save(load_name_intp0,'-v7.3','-struct','field_to_save'); (Saving is not nessecary and thus commented out to save disk space)
    return
    % 3. Load file with normal B-fields
else
    if par.APPLY_RMP && par.APPLY_TFR && ~exist(RMP_TFR_file,'file')
        [m    ,d.n3D.symm]=make_RMP_TFR_file(m,RMP_TFR_file);
        m.n3D=remove_fields(m.n3D,[],load_fields);
    else
        [m.n3D,d.n3D.symm]=load_3D_file(load_name,load_fields);
    end
    if ~all(isfield(m.n3D,load_fields)); error('Required fields not present in RMP_TFR_file'); end;
       
    for BA=BA_fields
        if par.superimpose_2D_3D && strcmp(par.coord_syst,'toroidal') && any(strcmp(BA{1},{'BR','BZ','Bphi'})) && ~par.APPLY_SAWTOOTH
            % Combine 2D and 3D field in toroidal case if no sawtooth applies
            m.n3D.(BA{1})=bsxfun(@plus,m.n3D.(BA{1})(1:par.step.R:d.size_X,1:par.step.Z:d.size_Z,1:par.step.phi:end),m.B_2D.(BA{1}));
        elseif strcmp(par.coord_syst,'toroidal')
            if ~par.superimpose_2D_3D; warning('2D and 3D FIELDS SUPERPOSITION DISABLED'); end
            m.n3D.(BA{1})=m.n3D.(BA{1})(1:par.step.R:d.size_X,1:par.step.Z:d.size_Z,1:par.step.phi:end);
        end
    end
    % Dimensions
    d.n3D.size_3D=size(m.n3D.(BA_fields{1})); % size 3D-map
    phi_scale   =linspace(0,2*pi/d.n3D.symm,d.n3D.size_3D(3));
    d.Dphi  =phi_scale(2)    -phi_scale(1);
    d.Dphi_inv=1/d.Dphi;
end


%% NESTED FUNCTIONS
    function intp0name=get_intp0_filename(filename)
        % Extend the load name with -intp0 to identify the maps with fields for interpolation method 0
        loc_extension=strfind(filename,'.mat');
        if isempty(loc_extension)
            intp0name=strcat(filename,'-intp0');
        else
            intp0name=strcat(filename(1:loc_extension-1),'-intp0','.mat');
        end
    end
end

%% Function to make RMP_TFR file
function [m,symm]=make_RMP_TFR_file(m,RMP_TFR_file) %#ok<INUSD>
global par
warning('combined RMP / TFR file absent, trying superposition of both files')

% Check the present fields in both files, only those can be combined
field_to_combine=extract_similar_fields(par.RMP_file,par.TFR_file);
disp('RMP and TFR wil combine the following fields: ')
disp(field_to_combine);

% Load RMP
[n3D_RMP,symm_RMP]=load_3D_file(par.RMP_file,field_to_combine);
disp('RMP file loaded')

% Load TFR
[n3D_TFR,symm_TFR]=load_3D_file(par.TFR_file,field_to_combine);
disp('TFR file loaded')

% Find the symmetry in both files (tor mode of RMP<= 16)
if symm_TFR~=symm_RMP
    ratio_symm=symm_TFR/symm_RMP;
    if mod(ratio_symm,1); error('Symmetries/points in 3D fields do not match!'); end;
    TFR_indexes=1:size(n3D_TFR.(field_to_combine{1}),3);
    TFR_indexes=[repmat(TFR_indexes(1:end-1),1,ratio_symm) TFR_indexes(end)];
else
    TFR_indexes=1:size(n3D_TFR.(field_to_combine{1}),3);
end
symm=symm_RMP;

% Superposition both 3D fields
for BA=field_to_combine
    m.n3D.(BA{1})=n3D_RMP.(BA{1})+n3D_TFR.(BA{1})(:,:,TFR_indexes);
    
    % Save RAM by deleting the obsolete individual fields
    n3D_RMP.(BA{1})=[];
    n3D_TFR.(BA{1})=[];
end
disp('RMP and TFR files superpositioned')

% Save the resulting combined file
field_for_save=m.n3D;
field_for_save.symm=symm;
disp('Saving combined file')
% save(RMP_TFR_file,'-v7.3','-struct','field_for_save') (SAVING is not nessecary and thus commented out)
clear field_for_save
% disp('Combined file saved')

%% INCORPORATED FUNCTIONS
    function tot_fnames=extract_similar_fields(file_A,file_B)
        %extrace_similar_fields Find fieldnames that present in both file_A
        %and file_B
        vars_A = whos('-file',file_A); vars_A={vars_A.name};
        vars_B = whos('-file',file_B); vars_B={vars_B.name};
        tot_fnames=vars_A;
        for i=1:length(vars_A)
            if ~any(strcmp(vars_A{i},vars_B))
                tot_fnames{i}='Missing field';
            end
        end
        tot_fnames(strcmp('Missing field',tot_fnames))=[];
        tot_fnames(strcmp('symm',tot_fnames))=[];
    end
end


function m=make_potential_interpolation_maps_2D(dim,m)
%make_potential_interpolation_maps Map creation for divergence free interpolation
% Creates the variables for the algorithm of Hugo.
% Please note that the R and Z-dimensions are 1 smaller (on top / LFS),
% since the with this cell cannot be calculated.
% So if the first 3D-grid was [Nr xNz x Nphi] the
% new 3D-grid is [Nr-1 x Nz-1 x Nphi-1]

%% BETAS first
%                                 R       Z
m.n2D.Beta_r  =	-diff(m.psi_XZ  (1:end-1 ,:      ),1,2)*dim.DZ_inv;

m.n2D.Beta_z  = +diff(m.psi_XZ  (:       ,1:end-1),1,1)*dim.DX_inv;

m.n2D.Beta_phi= m.B_2D.Bphi     (1:end-1 ,1:end-1);

%% BETAS second

m.n2D.Beta_phi_r= diff(m.B_2D.Bphi(:,1:end-1),1,1);

m.n2D.Beta_phi_z= diff(m.B_2D.Bphi(1:end-1,:),1,2);

%% BETAS third
m.n2D.Beta_phi_r_z =+diff(m.B_2D.Bphi,1,2);
m.n2D.Beta_phi_r_z = diff(m.n2D.Beta_phi_r_z,1,1);

%% ALPHAS

m.n2D.Alpha_r= (-diff(m.psi_XZ(2:end,:),1,2)*dim.DZ_inv - m.n2D.Beta_r);
m.n2D.Alpha_z= (+diff(m.psi_XZ(:,2:end),1,1)*dim.DX_inv - m.n2D.Beta_z);

end

function m=make_potential_interpolation_maps_3D(dim,m)
%make_potential_interpolation_maps Map creation for divergence free interpolation
% Creates the variables for the algorithm of Hugo.
% Please note that the R and Z-dimensions are 1 smaller (on top / LFS),
% since the with this cell cannot be calculated.
% So if the first 3D-grid was [Nr xNz x Nphi] the
% new 3D-grid is [Nr-1 x Nz-1 x Nphi-1]

psi_3D=bsxfun(@minus,m.psi_XZ,bsxfun(@times,dim.scale_X'+dim.R0,m.n3D.Aphi));

%% BETAS first
%                                 R       Z       phi
m.n3D.Beta_r  =	-diff(m.n3D.AZ(1:end-1 ,1:end-1,:      ),1,3)/dim.Dphi ...
				-diff(psi_3D  (1:end-1 ,:      ,1:end-1),1,2)/dim.DZ;

m.n3D.Beta_z  = +diff(psi_3D  (:       ,1:end-1,1:end-1),1,1)/dim.DX ...
				+diff(m.n3D.AR(1:end-1 ,1:end-1,:      ),1,3)/dim.Dphi;

m.n3D.Beta_phi= -diff(m.n3D.AR(1:end-1 ,:      ,1:end-1),1,2)/dim.DZ ...
				+diff(m.n3D.AZ(:       ,1:end-1,1:end-1),1,1)/dim.DX;
m.n3D.Beta_phi=bsxfun(@plus,m.n3D.Beta_phi,m.B_2D.Bphi(1:end-1,1:end-1));

%% BETAS second
m.n3D.Beta_phi_r=-diff(m.n3D.AR(:,:,1:end-1),1,2)/dim.DZ;
m.n3D.Beta_phi_r= diff(m.n3D.Beta_phi_r     ,1,1)/dim.DX;
m.n3D.Beta_phi_r= bsxfun(@plus, m.n3D.Beta_phi_r,diff(m.B_2D.Bphi(:,1:end-1),1,1)./dim.DX);

m.n3D.Beta_phi_z= diff(m.n3D.AZ(:,:,1:end-1),1,2)/dim.DZ;
m.n3D.Beta_phi_z= diff(m.n3D.Beta_phi_z     ,1,1)/dim.DX;
m.n3D.Beta_phi_z= bsxfun(@plus, m.n3D.Beta_phi_z,diff(m.B_2D.Bphi(1:end-1,:),1,2)./dim.DZ);

m.n3D.Beta_r_z  =-diff(m.n3D.AZ(1:end-1,:,:),1,3)/dim.Dphi;
m.n3D.Beta_r_z  = diff(m.n3D.Beta_r_z       ,1,2)/dim.DZ;

m.n3D.Beta_z_r  =+diff(m.n3D.AR(:,1:end-1,:),1,3)/dim.Dphi;
m.n3D.Beta_z_r  = diff(m.n3D.Beta_z_r       ,1,1)/dim.DX;

m.n3D.Beta_z_phi=-diff(psi_3D  (:,1:end-1,:),1,3)/dim.Dphi;
m.n3D.Beta_z_phi= diff(m.n3D.Beta_z_phi     ,1,1)/dim.DX;

m.n3D.Beta_r_phi=+diff(psi_3D  (1:end-1,:,:),1,3)/dim.Dphi;
m.n3D.Beta_r_phi= diff(m.n3D.Beta_r_phi     ,1,2)/dim.DZ;

%% BETAS third
m.n3D.Beta_phi_r_z =+diff(m.B_2D.Bphi       ,1,2)/dim.DZ;
m.n3D.Beta_phi_r_z =+diff(m.n3D.Beta_phi_r_z,1,1)/dim.DX;
m.n3D.Beta_phi_r_z = repmat(m.n3D.Beta_phi_r_z,[1 1 size(m.n3D.Beta_r,3)]);

%% GAMMA
m.n3D.Gamma_r  =-diff(m.n3D.AR(:,2:end,:),1,3)/dim.Dphi;
m.n3D.Gamma_r  = diff(m.n3D.Gamma_r      ,1,1)/dim.DX;
m.n3D.Gamma_r  = (m.n3D.Gamma_r  +m.n3D.Beta_z_r  )/dim.DZ;

m.n3D.Gamma_z  =-diff(m.n3D.AZ(2:end,:,:),1,3)/dim.Dphi;
m.n3D.Gamma_z  = diff(m.n3D.Gamma_z      ,1,2)/dim.DZ;
m.n3D.Gamma_z  = (m.n3D.Gamma_z  -m.n3D.Beta_r_z  )/dim.DX;

m.n3D.Gamma_phi=-diff(psi_3D(:,2:end,:),1,3)/dim.Dphi;
m.n3D.Gamma_phi= diff(m.n3D.Gamma_phi      ,1,1)/dim.DX;
m.n3D.Gamma_phi= (m.n3D.Gamma_phi-m.n3D.Beta_z_phi)/dim.DZ;

%% ALPHA
m.n3D.Alpha_r=-(diff(m.n3D.AZ(2:end,1:end-1,:),1,3)/dim.Dphi + diff(psi_3D(2:end,:,1:end-1),1,2)/dim.DZ + m.n3D.Beta_r)/dim.DX;
m.n3D.Alpha_z=+(diff(m.n3D.AR(1:end-1,2:end,:),1,3)/dim.Dphi + diff(psi_3D(:,2:end,1:end-1),1,1)/dim.DX - m.n3D.Beta_z)/dim.DZ;

end

function [d,m]=sawtooth_collapse(d,m)
global par
if verLessThan('matlab','8.4')
    intp_method='cubic';
else
    intp_method='phchip';
end

% (Complete) time array
d.st.time=get_time(par.st.n_reconnection,par.st.n_relaxation,par.st.t_reconnection,par.st.t_relaxation,'linear');

p_recon=par.st.n_reconnection+1;    % Index where the reconnection ends
p_relax=par.st.n_reconnection+2;    % Index where the relaxation starts

%% SFL-coordinates (initial)
if length(d.psi_scale_correct)~=length(d.q_initial_profile'); error('Length of psi and q profiles not the same. Please correct input in DATA TOKAMAK folder'); end;

psi_0=d.psi_scale_correct;
d.st.sign_psi_pol=sign(psi_0(1));


if d.st.sign_psi_pol==-1
    psi_0=-d.psi_scale_correct+min(psi_0);
else
    psi_0=max(d.psi_scale_correct)-d.psi_scale_correct;  % Define psi as 0 zero in center of plasma and finite at edge
end
psi_0=-psi_0;

chi=cumtrapz(psi_0,d.q_initial_profile);                           % toroidal flux (psi) = int_0^psi  q d(psi) with trapezium integration

m.chi_XZ=interp1(1:d.NB_PSI,chi,m.psi_norm_XZ,intp_method);      % Make a 2D map of chi, to determine toroidal flux
% psi_0=d.psi_scale_correct;
d.st.psi_star=psi_0-chi;                               % for m=n=1, psi_star = psi-phi_H, with helical poloidal flux p_H = n/m chi

% d.st.psi_star=d.st.sign_psi_pol*d.st.psi_star;

%% RADIAL POSITIONS
d.st.r  =   sqrt(2*abs(chi));  % Definition of r-coordinate

% r at q=1 (with top of psi_star)
ind_r0      = interp1(d.q_initial_profile,1:length(d.q_initial_profile),1,intp_method);     % Index q=1
ind_r0_int  = floor(ind_r0);                                                                % Integer index q=1
h_ind_r0_int= floor(ind_r0/2);                                                              % Integer index between center and q=1, for finding rmix
r0          = interp1(d.q_initial_profile,d.st.r,1,intp_method);                            % r-value at q=1

% Mixing radius (psi_star = 0). Interpoldation from h_ind_r0_int since we do not want the center point
d.st.ind_rmix = interp1(d.st.psi_star(h_ind_r0_int:end), (h_ind_r0_int:length(d.psi_scale_correct)	),0,intp_method);     % index mixing radius
rmix          = interp1(d.st.psi_star(h_ind_r0_int:end),d.st.r(h_ind_r0_int:end                             ),0,intp_method);     % r-valua mixing radius


%% RECONNECTION r-parameters (based on interpolation / finding equal flux)
% Assumption: r1 (position from original axis to center `old' core) moves with constant speed
% this is where to modify to change reconneciton dynamics
d.st.r1_t=@(t) r0*(1-t/par.st.t_reconnection);
d.st.r1=zeros(size(d.st.time));
d.st.r1(1:p_recon)=d.st.r1_t(d.st.time(1:p_recon));

% Find psi reconnection (go up in psi,r-plot)
d.st.psi_plus=interp1(d.st.r(1:ind_r0_int+1),d.st.psi_star(1:ind_r0_int+1),d.st.r1,intp_method);

% Find corresponding r2 (go right in psi,r-plot
d.st.r2=interp1(d.st.psi_star(ind_r0_int+1:end),d.st.r(ind_r0_int+1:end),d.st.psi_plus,intp_method);
d.st.r2(1)=d.st.r1(1);  % Make sure the reconnections starts at r0 for both parameters

% Check if r2 has reasonable values
if any(d.st.r2<d.st.r1); error('radial positions badly conditioned'); end
if (any(d.st.r2>rmix) || d.st.r2(end)~=rmix); error('Check r2 vector for conditioning'); end
d.st.r2(p_relax:end)=rmix;
% Position of new flux contour:
d.st.r3=sqrt(d.st.r2.^2-d.st.r1.^2);


%% FLUX VALUES
% q-values to determine dpsi / dr (which is dpsi/dchi * dchi/dr, with chi=0.5*r^2)
q1 =    interp1(d.st.r,d.q_initial_profile,d.st.r1,intp_method);
q2 =    interp1(d.st.r,d.q_initial_profile,d.st.r2,intp_method);

% dpsi/dr based on r1 and q-values:
psi_1_acc= q2.*(1-q1) .*d.st.r1;
psi_2_acc= q1.*(1-q2) .*d.st.r2;

%% SAWTOOTH PARAMETERS  (K and KC)
% RECONNECTION with KC
% KC: reconnection parameter equal to deformation with continuous E-field / no signular potential.
% Note this parameter is defined by the commented line, yet calculated by
% the next to remove singularity at r1=0 (end reconnection).
%kc =2./(r2+r1) .* (psi_1_acc + psi_2_acc)./(psi_2_acc./r2 - psi_1_acc./r1);
d.st.kc =2./(d.st.r2+d.st.r1) .* (psi_1_acc + psi_2_acc)./(q1 - q2);
d.st.kc(1)=0;

% RELAXATION WITH K -> 0
% Let k decrease to 0
d.st.k_t=@(t) d.st.kc(p_recon)*(1-(t-par.st.t_reconnection)/(par.st.t_relaxation));
d.st.k	=cat(1,d.st.kc(1:p_recon)    ,d.st.k_t(d.st.time(p_relax:end)));
d.st.k(end)=0;  % Make sure the final value of k is 0

% Determine kr based on geometry
d.st.kr=(d.st.r2-d.st.r1)./(d.st.r2+d.st.r1);
d.st.kr_inv=1./d.st.kr;

%% Time derivatime values
% if dynamics is modified above, this part will need to be modified as well
d.st.r1_dot=zeros(size(d.st.r1));
d.st.r1_dot(1:p_recon)=-r0/par.st.t_reconnection;       % Constant during relaxation and zero during reconnection (and beyond).

d.st.r2_dot     =d.st.r1_dot .* psi_1_acc./psi_2_acc;
d.st.r2_dot(1)	= -d.st.r1_dot(1);

d.st.r3_dot =(d.st.r2_dot.*d.st.r2-d.st.r1_dot.*d.st.r1)./d.st.r3;

d.st.kr_dot = 2*(d.st.r1.*d.st.r2_dot-d.st.r2.*d.st.r1_dot)./(d.st.r2+d.st.r1).^2;
d.st.kr_dot(p_relax:end)=0;     % Hard core krdot =0 after reconnection

d.st.k_dot  =gradient(d.st.k,d.st.time);
d.st.k_dot (p_relax:end)=-d.st.kc(p_recon)/par.st.t_relaxation;

d.st.k_m_kr_dot = d.st.k_dot-d.st.kr_dot;

%% Stable time (numerical errors)
disp(['Linear interpolating field till timestamp ',num2str(par.st.stable_time)])

%% Function to determine time spacing
    function time=get_time(n_1,n_2,trec,trel,option)
        switch option
            case 'linear'
                %% Linear time growth
                time=[linspace(0,trec,n_1+1),linspace(trec,trec+trel,n_2+1)]';
            case 'normal'
                %% Normal time growth (i.e. more points at start and transition
                norm_dist=@(mu,sigma,x) 1/sqrt(2*pi) / sigma * exp(-(x-mu).^2/(2*sigma^2)); % Normal distribution function
                
                % Distribution function with 0.5 constant and 0.5 normally distributed functions
                f_1=@(t) 0.5/trec+0.5*norm_dist(0,0.1*trec,t)+0.5*norm_dist(trec,0.1*trec,t);
                f_2=@(t) 0.5/trel+norm_dist(trec,0.1*trel,t);
                
                % Time array reconnection
                t_1=linspace(0,trec,1e7);   % Start with linear time distribution
                for i=1:10
                    c=cumtrapz(t_1,f_1(t_1));   % CDF (Cumulative Distribution function
                    c(end)=1;                   % Hardcode such that all time points are included
                    
                    t_1=interp1(c,t_1,linspace(0,1,n_1+1));   % Linearly in CDF distributed time points
                end
                if any(c>1); error('Not enough original points for proper CDF approximation'); end
                
                % Time array relaxation
                t_2=linspace(trec,trec+trel,1e7);  % Start with linear time distribution
                for i=1:10
                    c=cumtrapz(t_2,f_2(t_2));   % CDF (Cumulative Distribution function
                    c(end)=1;                   % Hardcode such that all time points are included
                    
                    t_2=interp1(c,t_2,linspace(0,1,n_2+1));   % Linearly in CDF distributed time points
                end
                if any(c>1); error('Not enough original points for proper CDF approximation'); end
                
                time=[t_1 t_2]';
                
        end
        time(n_1+1)=time(n_1+1)-0.5*par.dt;     % Set the last time point in reconnection phase just before the relaxation, so integration will not notice the jump. 
    end
end