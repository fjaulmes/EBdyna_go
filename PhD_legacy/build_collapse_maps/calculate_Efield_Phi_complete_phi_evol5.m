
close all;

REINIT_ALL_TOKAMAK_DATA=1;


if REINIT_ALL_TOKAMAK_DATA==1
    clear all;
    initialize_folder_names;
    initialize_collapse_map_calculation_context;
    
    rescaling_to_XZsmall_maps;
    
end

PROCESS_NUMBER=5



%simulation options
CALCULATE_EXB_DATA_FILE=0;
SAVE_DATA_FILE=1;

DT_INTERPOLATION_METHOD='quadratic'     % by default


if SAVE_DATA_FILE==0
    disp('no data file for this simulation')
end



disp('**************************************************************');


% initialize the entire phi maps for ExB speeds
vExB_X_map_phi=zeros(NB_PHI_DATA_HALF+1,NP,size_r);
vExB_Z_map_phi=zeros(NB_PHI_DATA_HALF+1,NP,size_r);
vExB_phi_map_phi=zeros(NB_PHI_DATA_HALF+1,NP,size_r);
Efield_X_map_phi=zeros(NB_PHI_DATA_HALF+1,NP,size_r);
Efield_Z_map_phi=zeros(NB_PHI_DATA_HALF+1,NP,size_r);
Efield_phi_map_phi=zeros(NB_PHI_DATA_HALF+1,NP,size_r);
grad_Phi_tor_map_phi=zeros(NB_PHI_DATA_HALF+1,NP,size_r);


%for frame_rank=1:101
for frame_rank=FLIST(PROCESS_NUMBER):FLIST(PROCESS_NUMBER+1)-1
    
    calculate_Efield_frame_rank;
    
    if (f<10)
        frame_name='00';
    elseif (f<100)
        frame_name='0';
    else
        frame_name='';
    end
    frame_name=strcat(frame_name,num2str(f));
    filename='../E_maps/E0';
    filename=strcat(filename,frame_name,'.mat');
    
    if (SAVE_DATA_FILE==1)
        if CALCULATE_EXB_DATA_FILE==1
            save(filename,'vExB_X_map_phi','vExB_Z_map_phi','Efield_X_map_phi','Efield_Z_map_phi','grad_Phi_tor_map_phi');
            disp('****** SAVING FILE in ../E_maps/ ****');
        else
            save(filename,'Efield_X_map_phi','Efield_Z_map_phi','grad_Phi_tor_map_phi');
        end
    end
    disp('**************************************************************');
    
    
end


if SAVE_DATA_FILE==1
    disp('****** data is available in ../E_maps/ ****');
else
    disp('****** size of data ****');
    size(E_potential_PR_map_phi)
end



