function change_ID(ID_old,ID_new,save_folder)
%change_ID Function to change the ID of a simulation
%   Function to alter the ID of a simulation. To be used incidently. 
% 	Finds all files with the old ID in their name in a certain folder and renames them accordingly.
% 	This function will also update the par-struct.

switch nargin
    case 0
        ID_old='1424'; 				% Enter old ID here
        ID_new='1432'; 				% Enter new ID here
        save_folder='./output_3/'; 	% Specify output folder
    case 2
        save_folder=['./output_1/'];
    case 3
    otherwise
        error('Not enough input arguments');
end


files=dir(save_folder);
fnames=[];
for i=3:length(files)
    if ~isempty(strfind(files(i).name,ID_old)) && isempty(strfind(files(i).name,'filepart'))
        fnames{end+1}=strcat(save_folder,files(i).name);
    end
end
if isempty(fnames)
    error('No files for ID change found')
end
fnames_new=fnames;
fields_ID={'ID','ID_NAME','SAVENAME_STATS','SAVENAME_RAW','SAVENAME'};
for i=1:length(fnames)
    new=load(fnames{i});
    
    for j=1:length(fields_ID)
        strt_symb=strfind(new.par.(fields_ID{j}),ID_old);
        if isempty(strt_symb) || length(strt_symb)>1
            error('ID not found')
        end
        new.par.(fields_ID{j})=[new.par.(fields_ID{j})(1:(strt_symb-1)),ID_new,new.par.(fields_ID{j})((strt_symb+length(ID_old)):end)];
    end
    strt_symb=strfind(fnames_new{i},ID_old);
    fnames_new{i}=[fnames_new{i}(1:(strt_symb-1)),ID_new,fnames_new{i}((strt_symb+length(ID_old)):end)];
    
    save(fnames_new{i},'-v7.3','-struct','new');
    disp(['SAVED ',num2str(i)])
    delete(fnames{i});
end
disp('FINISHED')