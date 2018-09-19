function [  ] = BS_RMP_combine_output(type_field,types)
global par BA_fields
%combine_output Combines files ending with p1 to a total file
%   After parallel processing of BS3D we would like to make maps of the
%   total fields. Therefore we combine the files
%#ok<*NASGU>


%% PLEASE PROVIDE THE REQUESTED MODE / CONFIGURATION OF THE RMP COILS (defined in BS_COMPASS_RMP_add_sings)
tor_modes=[2];
parity_names={'odd'};


%% OUTPUT PARAMETERS / load combined file of single coil config
save_folder='./output_RMP/';
save_name_type=@(type) strcat('./output/RMP_toroidal_single_coil_',type,'_',datestr(now,'yyyy-mm-dd'),'.mat'); % Native support for single coil mode with previous combined file (1 file to rule them all).
save_name_mode=@(tor_mode,p,par,type) strcat('./output/RMP_n=',num2str(tor_mode),'_',p,'_',par.coord_sys,'_',datestr(now,'yyyy-mm-dd'),'_',type,'.mat');

% Load coil cofiguration
[ RMP_coils] = BS_COMPASS_RMP_load_coils;

if nargin<2 || isempty(types)
    types={'Bu','Bl','Au','Al'};
end
% Combine the two types seperately
for type=types
    
if nargin>0 && ~isempty(type_field) || exist(save_name_type(type{1}),'file')
    warning('Combined file found of single coil vacuum field. Periodic repetition with requested modes will produce output file')
    mode='single_coil_mode';
    
    if nargin==0 || isempty(type_field)
        type_field=load(save_name_type(type{1}));
    end
        
    %% Find all files to combine for parameters
    if isempty(par)
        files=dir(save_folder);
        fnames=[];
        for i=3:length(files)
            if ~isempty(strfind(files(i).name,'individual')) && ~isempty(strfind(files(i).name,'RMP'))
                fnames{end+1}=strcat(save_folder,files(i).name); %#ok<*AGROW>
            end
        end
        
        % Load first file
        base=load(fnames{1}); % Skip 1 and 2 (. and ..)
        par=base.par;
        clear fnames
    end
    %% Find calculated field
    if isempty(BA_fields)
        BA_fields=fieldnames(base.field.(type{1}))';     % Find which fields have been determined (BR,BZ,Bphi)
        %BA_fields={'BR','BZ','Bphi','AR','AZ','Aphi','dAR_dphi','dAZ_dphi','dAphi_dphi'};
    end
else
    %% Find all files to combine
    files=dir(save_folder);
    fnames=[];
    for i=3:length(files)
        if ~isempty(strfind(files(i).name,'individual')) && ~isempty(strfind(files(i).name,'RMP'))
            fnames{end+1}=strcat(save_folder,files(i).name); %#ok<*AGROW>
        end
    end
    
    % Load first file
    base=load(fnames{1}); % Skip 1 and 2 (. and ..)
    par=base.par;
    par.N_PROCESS=200
    
    %if length(fnames)~=par.N_PROCESS
	%    fnames
    %    error('Failed to match number of processes with files')
    %end
    
    %% Find calculated field
    BA_fields=fieldnames(base.field.(type{1}))';     % Find which fields have been determined (BR,BZ,Bphi)
    %BA_fields={'BR','BZ','Bphi','AR','AZ','Aphi','dAR_dphi','dAZ_dphi','dAphi_dphi'};
    % Check if the loaded coils have (at least one) magnetic field determined
    for c=RMP_coils.coil_name
        if ~isfield(base.field,c) && ~isempty(RMP_coils.(c{1}))
            error('Use BS-solver first to find the field of the requested field')
        end
    end
    
    % Check how many coils of each type have been determined and exclude something other than 1 or 16
    
    nr_coils=length(base.field.(type{1}));
        
    switch nr_coils(1)
        case 1
            mode='single_coil_combination';
            if mod(par.size_total(3)-1,8)~=0;    error('Positions RMP coils / phi-positions don''t match');   end;
        case 4
            mode='all_coils';
        case 8
            mode='all_coils';
        otherwise
            error('The number of coils is neither 1 nor 8. Determine all coils or just 1 of each to use periodic repetition of vacuum field (RMP mode still acts')
    end
    disp(['EXECUTING IN MODE: ',mode])
end

%% Combine the output files (if all coils, apply mode right away)
switch mode
    case 'single_coil_combination'
        %% SIGNLE COIL COMBINATION
        % Pre-allocate struct with total field for each type
        % Add all files for each type seperately. Use RMP configuration and periodicy later on.
        save_name_single_coil=save_name_type(type{1});
        disp(['Making ',save_name_single_coil,' file'])
        for BA=BA_fields
            type_field.(BA{1})=zeros(par.size_total);    % Pre-allocate
        end
        for i=1:length(fnames)
            new=load(fnames{i}); % Load the field and parameters
            for BA=BA_fields
                type_field.(BA{1})(new.par.indexes) = new.field.(type{1})(1).(BA{1})(:);
            end
			disp(['Succesfully added ',num2str(i),' of ',num2str(par.N_PROCESS),' at: ',num2str(i*100/length(fnames),'%2.0f'),'%'])
        end
        
        %% SAVE field of each type of coil
        save(save_name_single_coil,'-v7.3','-struct','type_field');
        disp(['Saved the combined file with a single coil for type: ',type{1}])
        
        BS_RMP_combine_output(type_field,type);
        continue
    case 'single_coil_mode'
        %% APPLY THE MODE
        
        % The first and last element are duplicates, so therefore do
        % are copied twice.
        % EXAMPLE: 4 coils and 9 phi positions of grid (8 distinct)
        % index_grid:       [1      2       3       4       5       6       7       8       9   ]
        % 1st coil:         [1/9    2       3       4       5       6       7       8       1/9 ]
        % 2nd coil:         [7      8       1/9     2       3       4       5       6       7   ]
        % 3rd coil:         [5      6       7       8       1/9     2       3       4       5   ]
        % 4th coil:         [3      4       5       6       7       8       1/9     2       3   ]
        % i.e. formula if one has only coil 1, depends on nr that is doubled (call it y), with:     y=9-(coil-1)*2
        %       Then we copy values from y:8 and 1:y.       (e.g. for coil 2: y=7, for coil 3: y=5, for coil 4: y=3)
        
        total_phi_elements=par.size_total(3);                  % # total elements                      [9]
        translate_phi=(total_phi_elements-1)/8;                     % # phi-positions between RMP coils     [8/4=2]
        
        for tor_mode=tor_modes
            % Reduce size combined grid, to occompany mode
            if mod(total_phi_elements-1,tor_mode)==0
                symm=tor_mode;
            else
                symm=1;
            end
            size_total=par.size_total;
            size_total(3)=(total_phi_elements-1)/symm+1;     % Cut the total combined field to this value [5 with n=2] Other half is symmetric
            total_field.symm=symm;
            for p=parity_names
                if tor_mode==4 && any(strcmp(p{:},{'+90','-90'}))
                    continue
                end
                save_name=save_name_mode(tor_mode,p{1},par,type{1});
                % Load RMP configuration
                RMP_coils=BS_COMPASS_RMP_add_signs(RMP_coils,tor_mode,p{1});
                                
                % Pre-allocate the total field
                for BA=BA_fields
                     total_field.(BA{1})=zeros(size_total);
                     for c=1:8  
                        % Find coil parameters, i.e. sign and which indexes to copy
                        sign=RMP_coils.(type{1})(c).sign;
                        
                        phi_translation= (c-1) * translate_phi;        % Phi-positions that are translated in respect to the first
                        indexes_translated = mod( (1:size_total(3)) - phi_translation , total_phi_elements-1);  % Start with [1:5], reverse by 2, with correct indexes
                        indexes_translated(indexes_translated==0)=total_phi_elements-1;
                        
                        % Superposition of each field (with proper sign)
                        total_field.(BA{1})=total_field.(BA{1})...
                            +sign*type_field.(BA{1})(:,:,indexes_translated);
                     end
                end
                
                %% SAVE
                save(save_name,'-v7.3','-struct','total_field');
                disp(['Combined output file created in mode: ',mode,' tor_mode: ',num2str(tor_mode),' parity: ',p{1}])
            end
        end
        continue
    case 'all_coils'
        %% all coil mode
        for tor_mode=tor_modes
            for p=parity_names
                if tor_mode==4 && any(strcmp(p{:},{'+90','-90'}))
                    continue
                end
                save_name=save_name_mode(tor_mode,p{1},par,type{1});
                
                % Pre-allocate the total field
                for BA=BA_fields
                    total_field.(BA{1})=zeros(par.size_total);
                end
                % Load RMP configuration
                RMP_coils=BS_COMPASS_RMP_add_signs(RMP_coils,tor_mode,p{1});
                
                % Combine the files with this RMP configuration
                for i=1:length(fnames)
                    new=load(fnames{i}); % Load the field and parameters
                    for BA=BA_fields
                        for c=1:length(new.field.(type{1}))
                            % Add them with the proper sign
                            sign_coils=RMP_coils.(type{1})(c).sign;
							data_coils=new.field.(type{1})(c).(BA{1});
							data_coils=data_coils';
                            total_field.(BA{1})(new.par.indexes)=total_field.(BA{1})(new.par.indexes)+...
                                sign_coils*data_coils;
                        end
                    end
                    
                    disp(['Succesfully added ',num2str(i),' of ',num2str(par.N_PROCESS),' at: ',num2str(i*100/length(fnames),'%2.0f'),'%'])
                end
                
                %% SAVE all-COIL
                save(save_name,'-v7.3','-struct','total_field')
                disp(['Combined output file created in mode: ',mode,' tor_mode: ',num2str(tor_mode),' parity: ',p{1}])
            end
        end
        continue
end
end
%% Combine the B-fields of all coils
if nargin==0
    clearvars -except tor_modes parity_names BA_fields save_name_mode types par
    for tor_mode=tor_modes
        for p=parity_names
            if tor_mode==4 && any(strcmp(p{:},{'+90','-90'}))
                continue
            end
            save_name_1=save_name_mode(tor_mode,p{1},par,types{1});
            symm=load(save_name_1,'symm');
			if isempty(fieldnames(symm));
				symm=1;
			else
				symm=symm.symm;
			end
            final_save_name=[save_name_1(1:end-(length(types{1})+5)),'.mat'];
            for type=types
                save_name=save_name_mode(tor_mode,p{1},par,type{1});
                field.(type{1})=load(save_name,BA_fields{:});
            end
            for BA=BA_fields
                for i=2:length(types)
                    field.(types{1}).(BA{1})=field.(types{1}).(BA{1})...
                        +field.(types{i}).(BA{1});
					field.(types{i}).(BA{1})=[];
                end
            end
            field=field.(types{1});
			field.symm=symm;
            save(final_save_name,'-v7.3','-struct','field')
            disp('Combined complete mode')
        end
    end
end
end