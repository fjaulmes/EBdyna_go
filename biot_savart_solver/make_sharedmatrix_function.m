function make_sharedmatrix_function
%%make_sharedmatrix_function Installs the functions called sharedmatrix
% Identify the windows or linux machines and use the install file 
% accordingly. 
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

%% Switch to compilation path
path_shared_matrix=which('sharedmatrix.m');
old_wd=pwd;
cd(path_shared_matrix(1:end-14))

%% Determine the system and install accordingly
switch mach
		case 'WINDOWS';
		try
			c=rand(3);
			SharedMemory('clone','key_test',c);
			clear c
			SharedMemory('free','key_test');
		catch
			warning('(re)installing SharedMemory function')
			cd ./sharedmatrix_windows_by_andrew_smith
			SharedMemory_install;
		end
		
		case 'LINUX';
		try
			c=rand(3);
			sharedmatrix('clone',11111,c);
			clear c
			sharedmatrix('free',11111);
		catch
			warning('(re)installing sharedmatrix function')
			sharedmatrix_install;
		end
end
cd(old_wd);
end