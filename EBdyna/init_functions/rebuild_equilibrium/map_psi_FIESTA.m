
% Evaluating the magnetic equilibrium configuration
% according to output from the EFIT equilibrium code

clc

clear all
close all
format compact

TIME_EQUIL=1.55 %s

ADD_FIESTA_TO_FILENAMES=0;
FIESTA_FILENAME='FIESTA_equil';
% FIESTA_FILENAME_PARAMS=[FIESTA_FILENAME '_params.mat'];
FIESTA_FILENAME=[FIESTA_FILENAME '.mat']

METIS_FILENAME='./CU_32430s_vITER2.mat'


disp('... Initializing paramters...');
initialize_folder_names;
run('../create_physics_constants.m');


Nradial=129;
NP=129;
% DP=(2*pi)/(NP-1);
% for(p=1:NP)
%     theta_scale(p)=(p-1)*DP;    
% end

create_tokamak_parameters_FIESTA;
initialize_XZ_maps_FIESTA;

for (x=1:NR)
    for (z=1:NZ)
        Rpos_map(x,z)=Rpos(x);
    end
end

disp('Building PR maps .... ');
build_PR_maps_FIESTA;




%%

disp('Interpolate XZ maps from FIESTA RZ data .... ');
interpolate_iso_XYZ_maps_FIESTA;

%needs an approximate theta map for orbit classification at least
calculate_ki_angle_geometry;
theta_XZ_map=ki_XZ_map;

%%


inf_X=X1_out-1;
sup_X=X2_out+1;
inf_Z=1;
sup_Z=NZ;

if ADD_FIESTA_TO_FILENAMES==1
    FILENAME=strcat(DATA_FOLDER,'tokamak_map_dimensions_FIESTA.mat')
else
    FILENAME=strcat(DATA_FOLDER,'tokamak_map_dimensions.mat')
end
save (FILENAME,'NR','NZ','sup_X','sup_Z','inf_X','inf_Z','mid_X','mid_Z','DX','Rpos_XZ_map','Rpos','X_scale','Z_scale','Nradial','NP');




%%


save_maps_FIESTA
