function make_sharedmatrix_function
%%make_sharedmatrix_function Installs the functions called sharedmatrix
% Identify the windows or linux machines and use the install file 
% accordingly. 

%% Determine compilation paths
compile_path=fileparts(mfilename('fullpath'));  	% Use directory of current function
old_wd 		=pwd;

%% Test and otherwise install function
try
	c=rand(3);
	
	if ispc
		SharedMemory('clone','key_test',c);
		clear c
		SharedMemory('free','key_test');
	elseif isunix
		sharedmatrix('clone',11111,c);
		clear c
		sharedmatrix('free',11111);
	end
	return
catch
	if ispc
		warning('(re)installing SharedMemory function')
		cd(fullfile(compile_path,sharedmatrix_windows_by_andrew_smith))
		SharedMemory_install;
	elseif isunix
		warning('(re)installing sharedmatrix function')
		cd(compile_path)
		sharedmatrix_install;
	end	
end

%% Revert to old directory
cd(old_wd);
end