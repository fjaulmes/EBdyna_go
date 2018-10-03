function make_ba_interp3_function
%%make_ba_interp3_function Installs the functions called ba_interp3


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
	return
catch
	warning('(re)installing ba_interp3 function')
end

%% Determine compilation paths
compile_path=fileparts(mfilename('fullpath'));  	% Use directory of current function
old_wd 		=pwd;

%% Determine compilation file and move to cppFile
if ispc
	cppOriginFile = fullfile( compile_path , 'ba_interp3_windows.cpp' );
elseif isunix
	cppOriginFile = fullfile( compile_path , 'ba_interp3_linux.cpp' );
end

cppFile 	= fullfile( compile_path , 'ba_interp3.cpp');
copyfile( cppOriginFile , cppFile );

%% Compile
cd(compile_path);
mex -O ba_interp3.cpp
cd(old_wd)

%% Clean-up

delete( cppFile );
end