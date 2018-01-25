initialize_folder_names;

filename=strcat(DATA_FOLDER,'physics_constants.mat');
load(filename);
filename=strcat(DATA_FOLDER,'motions_map_dimensions.mat');
load(filename,'scale_X','scale_Z','R0');
filename=strcat(DATA_FOLDER,'XZsmall_fields_tokamak_pre_collapse.mat');
load(filename,'psi_XZsmall_map','q_initial_XZsmall_map');

FILENUMBER=200
SUBNUMBER=3
SUBTIME=2

ID='1915577'

figure;
hold on
contour(scale_X+R0,scale_Z,psi_XZsmall_map',[0 0],'b','linewidth',2);
% contour(scale_X,scale_Z,q_initial_XZsmall_map',[1.5 1.5],'r','linewidth',1);
% contour(scale_X,scale_Z,q_initial_XZsmall_map',[2.0 2.0],'r','linewidth',1);
% contour(scale_X,scale_Z,q_initial_XZsmall_map',[2.5 2.5],'r','linewidth',1);
% contour(scale_X,scale_Z,q_initial_XZsmall_map',[3.0 3.0],'r','linewidth',1);
axis xy equal square

for file_index=1:SUBNUMBER:FILENUMBER
    filename=['output_' num2str(ID) '/G_eq_19xx_full_process' num2str(file_index) '.mat']
    load(filename,'output')
    
    for time_index=2:size(output.x_gc,3)
        for part_index=1:size(output.x_gc,1)
            if ~isnan(output.x_gc(part_index,1,time_index))
                plot(output.x_gc(part_index,1,time_index),output.x_gc(part_index,2,time_index),'k.');
            end
        end
    end
end

contour(scale_X+R0,scale_Z,q_initial_XZsmall_map',[1.5 1.5],'r--','linewidth',1.3);
contour(scale_X+R0,scale_Z,q_initial_XZsmall_map',[2.0 2.0],'r--','linewidth',1.3);
contour(scale_X+R0,scale_Z,q_initial_XZsmall_map',[2.5 2.5],'r--','linewidth',1.3);
contour(scale_X+R0,scale_Z,q_initial_XZsmall_map',[3.0 3.0],'r--','linewidth',1.3);

