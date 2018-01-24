function make_ba_interp2_function
%%make_ba_interp2_function Installs the functions called ba_interp2

%% Test and otherwise install function
try
	x=linspace(0,1,10);
	y=linspace(0,1,10);
    [x,y]=meshgrid(x,y);
    
	g=x.^2+y.^2;
	
	X=linspace(0,1,50);
	Y=linspace(0,1,50);
    
	ba_interp2(g,X*10,Y*10,'cubic');
	return
catch
	warning('(re)installing ba_interp2 function')
end

%% Determine compilation paths
compile_path=fileparts(mfilename('fullpath'));  	% Use directory of current function
old_wd 		=pwd;

%% Determine compilation file and move to cppFile
if ispc
	cppOriginFile = fullfile( compile_path , 'ba_interp2_windows.cpp' );
elseif isunix
	cppOriginFile = fullfile( compile_path , 'ba_interp2_linux.cpp' );
end

cppFile 	= fullfile( compile_path , 'ba_interp2.cpp');
copyfile( cppOriginFile , cppFile );

%% Compile
cd(compile_path);
mex -O ba_interp2.cpp
cd(old_wd)

%% Clean-up

delete( cppFile );
end