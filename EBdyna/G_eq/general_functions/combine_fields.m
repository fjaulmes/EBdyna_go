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
	    new_length=length(new.(fname));
		total_f=total.(fname);
		if new_length<length(total.(fname))
        total_f(1:new_length)=total_f(1:new_length)+new.(fname);
		else
		total_f=total_f+new.(fname);
		end
		total.(fname)=total_f;
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
		    new.(fname)=new.(fname)';
			dim=find(size(new.(fname))==N_job);
			if dim==1
				make_size=[N_total,size(new.(fname),2),size(new.(fname),3)];
			else
				error('particles not ordered in first dimension');
			end
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
