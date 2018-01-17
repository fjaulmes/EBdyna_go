close all
clear all
initialize_folder_names
initialize_TAE_map_calculation_context;
TAE_parameters;
initialize_XZ_maps_dimensions;
rescaling_to_XZsmall_maps;


size_r_TAE=pTAE_sup-pTAE_inf+1
finesse_data_X=reshape((Rpos_PR_map(:,pTAE_inf:pTAE_sup)-R0),NP*size_r_TAE,1);
finesse_data_Z=reshape(Z_PR_map(:,pTAE_inf:pTAE_sup),NP*size_r_TAE,1);


%%

bphi_XZsmall_map=sqrt(1-bX_XZ_small_map.^2+bZ_XZ_small_map.^2);


clc
NB_FRAME=100

EX_PR_map=zeros(NP,Nradial);
EZ_PR_map=zeros(NP,Nradial);

W_TAE_evol=zeros(NB_FRAME,1);

for f=1:NB_FRAME
    f
    
    %assuming propagation at average Alfven velocity
    OMEGA_VA=omega_TAE;
    
    
    if (f<10)
        frame_name='00';
    elseif (f<100)
        frame_name='0';
    else
        frame_name='';
    end
    frame_name=strcat(frame_name,num2str(f));

    filenameB='../B_maps/B0';
    filenameB=strcat(filenameB,frame_name,'.mat');   
    filenameE='../E_maps/E0';
    filenameE=strcat(filenameE,frame_name,'.mat');   
    map_WTAE_frame_rank;
    
end



W_TAE_oscill_evol=W_TAE_evol;
% W_TAE_oscill_evol=W_TAE_evol/W_TAE_oscill_evol(1)


WTAE_AVG=mean(W_TAE_oscill_evol)


filename='./W_TAE_oscill_evol_global.mat';
save(filename,'W_TAE_oscill_evol');

filename='../data_tokamak/TAE_data.mat';
save(filename,'-append','WTAE_AVG')




