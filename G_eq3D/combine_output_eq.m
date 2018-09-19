function [] =combine_output_eq(ID,type_sim)
global par maps dim qom const
%combine_input_files Function that files from processors in equilibrium.
%   Combines files written by

format compact 

%% Looking for files
if nargin==0
    ID='1915611';
    type_sim='prec';
elseif nargin==1
    type_sim='full';
end

evaluate_matprec_files=true;
% if run locally used folder name with ID parameter
save_folder=['./output_',ID,'/']
% save_folder=['./output/']
type_sim

files=dir(save_folder)
fnames=[];
fnames_matprec=[];

disp('length(files)');
disp(length(files)-2)

for i=3:length(files)
    filename=files(i).name
    if strcmp(type_sim,'full') && isempty(strfind(filename,'matprec')) && isempty(strfind(filename,'filepart')) && isempty(strfind(filename,[ID,'_raw']))
		%disp(['scanning file #' num2str(i)]);
        fnames{end+1}=[save_folder,files(i).name]; %#ok<AGROW>
    elseif strcmp(type_sim,'prec') && isempty(strfind(filename,'matprec')) && isempty(strfind(filename,'filepart')) 
		%disp(['scanning file #' num2str(i)]);
        fnames{end+1}=[save_folder,files(i).name]; %#ok<AGROW>
    elseif (~isempty(strfind(filename,ID)) || ~isempty(strfind(filename,[ID,'_raw'])))&& ~isempty(strfind(filename,'matprec')) && isempty(strfind(filename,'filepart'))
		%disp(['scanning file #' num2str(i)]);
        fnames_matprec{end+1}=[save_folder,files(i).name]; %#ok<AGROW>
    end
end

nr_files=length(fnames)
%% Load first file
base=load(fnames{1}); % Skip 1 and 2 (. and ..)
par=base.par;
base_fields=fieldnames(base);

%% Determine output savenameswitch type_sim
begin_name=strfind(base.par.SAVENAME,'full');
save_name_full=[base.par.SAVENAME(1:(begin_name+length('full')-1)),'.mat'];
switch type_sim
    case 'full'
        save_name=save_name_full;
    case 'prec'
        begin_name=strfind(base.par.SAVENAME_STATS,type_sim);
        save_name=[base.par.SAVENAME_STATS(1:(begin_name+length(type_sim)-1)),'.mat'];
    case 'raw'
        begin_name=strfind(base.par.SAVENAME_RAW,type_sim);
        save_name=[base.par.SAVENAME_RAW(1:(begin_name+length(type_sim)-1)),'.mat'];
        
        % [maps,dim] = initialize_maps(true,false);
        % v_to_E=@(m,v) 0.5*(m/const.eV).*v.^2;
    otherwise
        begin_name=strfind(base.par.SAVENAME,'full');
        save_name=[base.par.SAVENAME(1:(begin_name-1)),'unknowntype.mat'];
end
disp('********************************************************************')
disp('Combined output will be saved in:')
disp(save_name)
disp('********************************************************************')

%% Extend fnames for simulations that have all particles ejected in prec-simulation
PROCESS_NUMBERS_IGNORED=[];
% this does not seem very useful and is pretty confusing
if evaluate_matprec_files
   rerun_script=false;
   for i=1:length(fnames_matprec)
       copyfile(fnames_matprec{i},[save_folder,'TEMP.mat']);
       new=load([save_folder,'TEMP.mat']);
       if all(new.ejected)
           warning(['No data off process ',num2str(new.par.PROCESS_NUMBER),' since every particle is ejected in equilibrium'])
           PROCESS_NUMBERS_IGNORED(end+1)=new.par.PROCESS_NUMBER; %#ok<AGROW>
	   else
%        else
            warning('UNEXPECTED MATPREC FILE; NOT ALL PARTICLES EJECTED')
%            if isempty(maps)
%                [maps,dim]            =initialize_maps(true,false);
%                qom = input.Z*const.eV/input.m;
%            end
%            % Extract the precession data from the matprec-file
%            new.output=evaluate_output(new.input,new.output,new.ejected);
%            new.prec = extract_precession_information_struct(new.input,new.output,new.ejected);
%
%            new.output=remove_fields(new.output,[],{'x','v','x_gc','nr_vpll_crossing'});
%            save(new.par.SAVENAME_STATS,'-v7.3','-struct','new');
%            movefile(new.par.SAVENAME_STATS',save_folder)
%            delete(fnames_matprec{i})
%            rerun_script=true;
        end
    end
%    if rerun_script
%        combine_output_eq(ID,type_sim);
%        return
%    end
    PROCESS_NUMBERS_IGNORED=unique(PROCESS_NUMBERS_IGNORED);
end


%% Begin combining files
%if (nr_files+length(PROCESS_NUMBERS_IGNORED))~=isempty(strfind(filename,type_sim))
%    error('Failed to match number of processes with files')
%end
total=[];
N_total=base.input.N_total;
total.process_time=[];

% WAITBAR
waitbar_on=false;
if waitbar_on
    delete(findall(0,'type','figure','tag',mfilename)); %#ok<UNRCH>
    set(0,'defaulttextinterpreter','tex')
    wb=waitbar(0,'Combining files','name','Combining output files','CreateCancelBtn','setappdata(gcbf,''canceling'',1)','tag',mfilename);
end

%% Loop over all the rest and add them in the same field
for index=1:nr_files
    % Loading
    new=load(fnames{index});        % Load the next file
    
    % Checking validity
    if any(new.par.PROCESS_NUMBER==PROCESS_NUMBERS_IGNORED)
        error(['Simulation box ',num2str(new.par.PROCESS_NUMBER),' is both present as output and as MATPREC'])
    end
    
    new_fields=fieldnames(new); % Fieldnames
    if ~isequal(new_fields,base_fields)
        error('Inconsistency in stored data fields')
    end
    % Add output files
    for j=1:length(new_fields)
        total=combine_fields(total,new,new_fields{j},new.input.particle_nr,new.input.N_job,N_total);
    end
    
    % Waitbar data
    if waitbar_on
        waitbar(index/nr_files,wb,...
            ['Combining files, @ ',num2str(index/nr_files*100,'%10.1f'),'%']) %#ok<UNRCH>
        if getappdata(wb,'canceling')
            delete(wb); set(0,'defaulttextinterpreter','latex')
            error('pressed cancel button')
        end
    end
    disp(['Combined file ',num2str(new.par.PROCESS_NUMBER),' of ',num2str(nr_files),'. At ',num2str(100*index/nr_files,'%10.1f'),'%'])
end

%% After combining operations / saving
switch type_sim
    case 'full'
        if isfield(total.output,'time_step_loss') && ~isfield(total.output,'loss')
            total.output.loss=NaN(base.par.NB_TIME_STEPS,1);
            for i=1:par.NB_TIME_STEPS
                total.output.loss(i)=sum(i==total.output.time_step_loss);
            end
        end
        save(save_name,'-v7.3','-struct','total');
    case 'raw'
        %% Extract E-profile
        E=v_to_E(total.input.m,sqrt(dot(total.output.v,total.output.v,2)));
        E_mean=NaN(1,dim.NB_PSI,par.NB_TIME_STEPS);
        
        for t=1:par.NB_TIME_STAMPS
            X_ind=((total.output.x(:,1,t)-dim.R0)*dim.DX_inv)+dim.mid_Xzero;
            Z_ind=( total.output.x(:,2,t)        *dim.DZ_inv)+dim.mid_Z;
            
            psi_norm=ba_interp2(maps(1).psi_norm_XZ       ,Z_ind,X_ind,'linear');
            
            E_matrix=bsxfun(@eq,floor(psi_norm(:,:)),0:dim.NB_PSI-1);
            E_matrix=bsxfun(@times,E_matrix,E(:,:,t));
            E_matrix(E_matrix==0)=NaN;
            for j=1:dim.NB_PSI
                E_to_mean=E_matrix(:,j);
                if all(isnan(E_to_mean))
                    continue
                end
                E_mean(1,j,t)=mean(E_to_mean(isfinite(E_to_mean)));
            end
            disp(['Profile of time ',num2str(t),' made. Now at ',num2str(t/par.NB_TIME_STEPS*100,'%10.1f'),'%'])
        end
        %% SAVE COMBINED FILE
        output=load(save_name_full,'output'); output=output.output;
        output.E_mean=squeeze(E_mean);
        save(save_name_full,'-append','output')
        % save(save_name,'-struct','total')
    case 'prec'
        prec=total.prec;
		prec.input=total.input;
        output=total.output;
		%output.vphi=output.v(:,3);
		%output=remove_fields(output,{'v'})
		%save('./output/total_file.mat','total','-v7.3')
		save('./output/output_file.mat','output','-v7.3')
		save('./output/prec_file.mat','-struct','prec','-v7.3')
        if exist(base.par.DISTNAME,'file')
            save(base.par.DISTNAME,'-append','prec','-v7.3');
            disp('SAVED PRECESSION WITHIN DISTRIBUTION FILE!')            
        else
            save(save_name)
            warning('Distribution file to save precession information does not exist')
        end
end

disp('Combined output file created!')

%% Clean up
if waitbar_on; delete(wb); set(0,'defaulttextinterpreter','latex'); end %#ok<UNRCH>
end

%% FUNCTIONS

function total=combine_fields(total,new,fname,part_nr,N_job,N_total)
if strcmp(fname,'ind')
    %% Couple indexes of all particle types
    fnames=fieldnames(new.ind);
    for j=1:length(fnames)
        if ~isfield(total,'ind')
            total.ind=[];
        end
        if ~isfield(total.ind,fnames{j})
            total.ind.(fnames{j})=[];
        end
        total.ind.(fnames{j})=cat(1,total.ind.(fnames{j}),new.ind.(fnames{j}));
    end
elseif isstruct(new.(fname))
    %% STRUCTS
    % New (doesn't exist already)
    if ~isfield(total,fname)
        total.(fname)=[];
    end
    % Old struct (not first file to combine)
    fnames=fieldnames(new.(fname));
    for j=1:length(fnames)
        total.(fname)=combine_fields(total.(fname),new.(fname),fnames{j},part_nr,N_job,N_total);
    end
elseif strcmp(fname,'process_time')       % These parameters are combined / cat
    %% PROCESS TIME
    total.process_time(new.par.PROCESS_NUMBER)=new.process_time;
elseif any(strcmp(fname,{'loss','nr_profile','profile_2D'}))    % These parameters are added
    %% VALUES WHICH SHOULD BE ADDED
    if ~isfield(total,fname)
        total.(fname)=new.(fname);
    else
        total.(fname)=total.(fname)+new.(fname);
    end
else
    %% VECTORS AND MATRICES
    dim=find(size(new.(fname))==N_job); % Dimension of particles
    % Copy single values once (if not process_time a pop value or)
    if isempty(dim)
        if ~isfield(total,fname)
            total.(fname)=new.(fname);
        end
        return
    end
    % New vector/matrix
    if ~isfield(total,fname)
        if dim==1
            make_size=[N_total,size(new.(fname),2),size(new.(fname),3)];
        else
            error('particles not ordered in first dimension')
        end
        if islogical(new.(fname))
            total.(fname)=true(make_size);
        else
            total.(fname)=NaN(make_size);
        end
    end
    % Copy data in matrix / vector
    total.(fname)(part_nr,:,:)=new.(fname);
end
end
