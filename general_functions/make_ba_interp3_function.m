function make_ba_interp3_function
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
compile_path=which('ba_interp3.cpp');
old_wd=pwd;
begin_str=strfind(compile_path,'ba_interp3.cpp');
cd(compile_path(1:begin_str-1))

%% Test and otherwise install function
try
	x=linspace(0,1,10);
	y=linspace(0,1,10);
    z=linspace(0,1,10);
    [x,y,z]=meshgrid(x,y,z);
    
	g=x.^2+y.^2+z.^2;
	
	X=linspace(0,1,50);
	Y=linspace(0,1,50);
	Z=linspace(0,0,50);
    
	ba_interp3(x,y,z,g,X,Y,Z,'cubic');
catch
	warning('(re)installing ba_interp3 function')
    switch mach
        case 'WINDOWS'
            movefile('./ba_interp3.cpp','./ba_interp3_LINUX.cpp')
            movefile('./ba_interp3_windows.cpp','./ba_interp3.cpp')
            mex -O ba_interp3.cpp
            movefile('./ba_interp3.cpp','./ba_interp3_windows.cpp')
            movefile('./ba_interp3_LINUX.cpp','./ba_interp3.cpp')
        case 'LINUX'
            mex -O ba_interp3.cpp
    end
end
cd(old_wd);
end