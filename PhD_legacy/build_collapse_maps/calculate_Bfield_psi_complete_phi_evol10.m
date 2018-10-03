
close all;

REINIT_ALL_TOKAMAK_DATA=1;


if REINIT_ALL_TOKAMAK_DATA==1
    clear all;
    initialize_folder_names;
    
    initialize_collapse_map_calculation_context;
    
    rescaling_to_XZsmall_maps
    
end

PROCESS_NUMBER=10



%simulation options
SAVE_DATA_FILE=1;
CALCULATE_VD_DATA_FILE=0;

DT_INTERPOLATION_METHOD='quadratic'     % by default


if SAVE_DATA_FILE==0
    disp('no data file for this simulation')
end


disp('**************************************************************');


% initialize the entire phi maps for ExB speeds

vD_X_map_phi=zeros(NB_PHI,NP,size_r);
vD_Z_map_phi=zeros(NB_PHI,NP,size_r);
vD_phi_map_phi=zeros(NB_PHI,NP,size_r);
Btot_map_phi=zeros(NB_PHI,NP,size_r);
bX_map_phi=zeros(NB_PHI,NP,size_r);
bZ_map_phi=zeros(NB_PHI,NP,size_r);
bphi_map_phi=zeros(NB_PHI,NP,size_r);

B2_tot_evol=zeros(101,1);
energy_phi_evol=zeros(101,NB_PHI);

for frame_rank=FLIST(PROCESS_NUMBER):FLIST(PROCESS_NUMBER+1)-1
    
    calculate_Bfield_frame_rank;
    
    if (f<10)
        frame_name='00';
    elseif (f<100)
        frame_name='0';
    else
        frame_name='';
    end
    frame_name=strcat(frame_name,num2str(f));
    filename='../B_maps/B0';
    filename=strcat(filename,frame_name,'.mat');
    
    if (SAVE_DATA_FILE==1)
        if CALCULATE_VD_DATA_FILE==1
            save(filename,'vD_X_map_phi','vD_Z_map_phi','E_potential_PR_map_phi','Btot_map_phi','bX_map_phi','bZ_map_phi','grad_psi_star_map_phi');
        else
            save(filename,'Btot_map_phi','bX_map_phi','bZ_map_phi','grad_psi_star_map_phi');
        end
        
        disp('****** SAVING FILE in ../B_maps/  ****');
    end
    %disp('**************************************************************');
    
    
end

save(strcat('energy_mag_evol_',num2str(PROCESS_NUMBER)),'B2_tot_evol','deltaB2_tot_evol','energy_phi_evol')


if SAVE_DATA_FILE==1
    disp('****** data is available in ../B_maps/  ****');
else
    disp('****** size of data ****');
    size(Btot_map_phi)
end



