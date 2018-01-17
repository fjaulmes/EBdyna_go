function load_RAM(remove)
%%load_RAM
% Will load in maps from I/O-file system and load this in the shared RAM.
% Multiple processes can use this. Make sure that each process detaches the
% pointer.

%% Machine identifier
global mach
if isempty(mach)
	if 	~isempty(strfind(computer(),'WIN')) || ~isempty(strfind(computer(),'win'))
		disp('Detected WINDOWS machine')
		mach='WINDOWS';
	else
		disp('Detected LINUX machine')
		mach='LINUX';
	end
end
narginchk(0,1);

%% Remove variable
switch nargin
    case 0
        remove=false;
end

switch remove
    case '0'
        remove=false;
    case '1'
        remove=true;
end

%%% Make ba_interp and sharedmatrix functions available
%old_path=pwd;
%cd ..
%% Compiling functions if neccesary
%make_sharedmatrix_function;
%make_ba_interp3_function;
%make_ba_interp2_function;
%cd(old_path)

%% Make the maps
if ~remove 
	clearvars -global -except mach
    if exist('shared_memory_keys.mat','file')
        delete('shared_memory_keys.mat')
    end
    G_eq([],[],true);
    disp('Making maps in shared memory succeeded! Keys stored in: ');
    which('shared_memory_keys.mat')
    return
else
    %% Remove the maps
    disp('Trying to free shared RAM')
    map_data=load('shared_memory_keys.mat');
    delete('shared_memory_keys.mat')
    key_names=fieldnames(map_data.keys);
    for i=1:length(key_names)
        switch mach
            case 'WINDOWS'
                SharedMemory('free',map_data.keys.(key_names{i}));        
            case 'LINUX'
                sharedmatrix('free',map_data.keys.(key_names{i}));
        end
    end
    disp('Shared RAM freed')
    
end


end
