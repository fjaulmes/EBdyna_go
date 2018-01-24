function make_ba_interp2_function
%%make_sharedmatrix_function Installs the functions called sharedmatrix
% Identify the windows or linux machines and use the install file 
% accordingly. 

%% Switch to compilation path
compile_path=which('ba_interp2.cpp');
old_wd=pwd;
begin_str=strfind(compile_path,'ba_interp2.cpp');
cd(compile_path(1:begin_str-1))

%% Test and otherwise install function
try
	x=linspace(0,1,10);
	y=linspace(0,1,10);
    [x,y]=meshgrid(x,y);
    
	g=x.^2+y.^2;
	
	X=linspace(0,1,50);
	Y=linspace(0,1,50);
    
	ba_interp2(g,X*10,Y*10,'cubic');
catch
	warning('(re)installing ba_interp2 function')
    if ispc
		movefile('./ba_interp2.cpp','./ba_interp2_LINUX.cpp')
		movefile('./ba_interp2_windows.cpp','./ba_interp2.cpp')
		mex -O ba_interp2.cpp
		movefile('./ba_interp2.cpp','./ba_interp2_windows.cpp')
		movefile('./ba_interp2_LINUX.cpp','./ba_interp2.cpp')
	elseif isunix
		mex -O ba_interp2.cpp
    end
end
cd(old_wd);
end