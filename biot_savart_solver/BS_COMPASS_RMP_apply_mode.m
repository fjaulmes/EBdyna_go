function [ RMP_coils ,field] = BS_COMPASS_RMP_apply_mode(RMP_map,RMP_coils )
%BS_AUG_RMP_apply_mode applies toroidal mode
% Asks for the sign from the add signs function. Then adds superpositions
% the B and A fields of the RMP coils. 

%% USE AS STANDALONE, FROM RAW DATA
if nargin==0
    tor_mode=2;
    parity='odd';
    %% Get file names individual files
    files=dir('./output/');
    fnames=[];
    for i=3:length(files)
        if ~isempty(strfind(files(i).name,'individual'))
            fnames{end+1}=strcat('./output/',files(i).name);
        end
    end
    %% Apply mode and save
    remove={'BX','BY','AX','AY'};
    for i=1:length(fnames)
        loaded_data=load(fnames{i});
        [ ~,field] = BS_COMPASS_RMP_apply_mode(loaded_data.RMP_map,loaded_data.RMP_coils,tor_mode,parity);
        par=loaded_data.par;
        field=remove_fields(field,remove);
        save_name_even2=strcat('./output/BS_COMPASS_RMP_n=',num2str(tor_mode),'_even_',par.coord_sys,'field_',datestr(now,'yyyy-mm-dd'),'_processmode',num2str(par.PROCESS_NUMBER),'.mat');
        save(save_name_even2,'-v7.3','field','par');
        disp(strcat('Saved n=',num2str(tor_mode),' file: ',num2str(loaded_data.par.PROCESS_NUMBER),' of ',num2str(loaded_data.par.N_PROCESS),'at : ',num2str(i*100/length(fnames),'%2.0f'),'%'))
        clear loaded_data field par save_name_even2
    end
    return
%% USE AS FUNCTION
else
    narginchk(2,2);
end

%% Add them with proper sign
if ~isfield(RMP_map,'Bu')
    error('Size script in this function uses Bu. Please note other fields if used differently')
end
fields_determined=fieldnames(RMP_map.Bu); % type of fields that have been calculated

for field_type=1:length(fields_determined) % for each field
    % Pre-allocate
    field.(fields_determined{field_type})=zeros(size(RMP_map.Bu(1).(fields_determined{field_type})));
    
    for type=1:3 % For each type of coil
        for i=1:length(RMP_coils.(RMP_coils.coil_name{type}))
            % Add to field
            field.(fields_determined{field_type})=field.(fields_determined{field_type})...
                +RMP_coils.(RMP_coils.coil_name{type})(i).sign*RMP_map.(RMP_coils.coil_name{type})(i).(fields_determined{field_type});
        end
    end
end

%% Determine phi and R-components from original Cartesian
% field.Aphi=-RMP_map.X./RMP_map.R.*field.AY+RMP_map.Y./RMP_map.R.*field.AX;
% field.Bphi=-RMP_map.X./RMP_map.R.*field.BY+RMP_map.Y./RMP_map.R.*field.BX;
% 
% field.AR=RMP_map.X./RMP_map.R.*field.AX+RMP_map.Y./RMP_map.R.*field.AY;
% field.BR=RMP_map.X./RMP_map.R.*field.BX+RMP_map.Y./RMP_map.R.*field.BY;

% field.AZ=field.AZ;
% field.BZ=field.BZ;
%% Remove obsolete fields
% remove={'BX','BY','AX','AY'};
% field=remove_fields(field,remove); 
return
end