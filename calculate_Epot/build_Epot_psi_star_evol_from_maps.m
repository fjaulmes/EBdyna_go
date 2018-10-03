% rewrap the content of /reconnection_maps
% in a condensed form in Epot_psi_star_dot_evol.mat

clear all;
initialize_folder_names;
filename=strcat(DATA_FOLDER,'tokamak_map_dimensions.mat');
load(filename);
filename=strcat(DATA_FOLDER,'tokamak_PR_map.mat');
load(filename);
filename=strcat(DATA_FOLDER,'psi_star_evol.mat');
load(filename);
filename=strcat(DATA_FOLDER,'psi_profiles_kadomtsev.mat');
load(filename);
initialize_Epot_calculation_parameters;

if ISKINK==1
	NB_FRAMES=100
else
	NB_FRAMES=101
end
NP_half=round((NP+1)/2);
Nomega=NP_half;
NB_PHI=129;
NB_THETA=NP;
THETA_PHI_RATIO=round((NP-1)/(NB_PHI-1))
DELTA_PHI=THETA_PHI_RATIO*PHI_STEP_SIZE;

Epot_omega_psi=zeros(NP,size_r);
psi_star_dot_omega_psi=zeros(NP,size_r);
Epot_Romega=zeros(size_r,NP);
psi_star_dot_Romega=zeros(size_r,NP);

Epot_evol=zeros(NB_FRAMES,size_r,NP);
psi_star_dot_evol=zeros(NB_FRAMES,size_r,NP);


for frame_rank=1:NB_FRAMES
    f=(frame_rank-1)+1
    if (f<10)
        frame_name='00';
    elseif (f<100)
        frame_name='0';
    else
        frame_name='';
    end
    frame_name=strcat(frame_name,num2str(f));
    filename=strcat(RECONNECTION_MAPS_FOLDER,'t0',frame_name,'.mat');
    load(filename);
    phi_rank=1;
    for omega_rank=1:DELTA_PHI:NB_THETA
        Epot_omega_psi=squeeze(E_potential_PR_map_phi(phi_rank,:,:));
        psi_star_dot_omega_psi=squeeze(psi_star_dot_PR_map_phi(phi_rank,:,:));
        Epot_omega_psi=[Epot_omega_psi(omega_rank:end,:) ; Epot_omega_psi(2:omega_rank,:) ];
        psi_star_dot_omega_psi=[psi_star_dot_omega_psi(omega_rank:end,:) ; psi_star_dot_omega_psi(2:omega_rank,:) ];
        E_potential_PR_map_phi(phi_rank,:,:)=Epot_omega_psi;
        psi_star_dot_PR_map_phi(phi_rank,:,:)=psi_star_dot_omega_psi;
        phi_rank=phi_rank+1;
    end
    
    Epot_Romega(:,:)=squeeze(mean(E_potential_PR_map_phi(1:end-1,1:NP,:),1))';
    psi_star_dot_Romega(:,:)=squeeze(mean(psi_star_dot_PR_map_phi(1:end-1,1:NP,:),1))';
%     Epot_Romega=Epot_Romega';
%     psi_star_dot_Romega=psi_star_dot_Romega';
    Epot_evol(frame_rank,:,:)=Epot_Romega(:,:);
    psi_star_dot_evol(frame_rank,:,:)=psi_star_dot_Romega(:,:);
end

FILENAME=strcat(DATA_FOLDER,'Epot_psi_star_dot_evol.mat')
save(FILENAME,'Epot_evol','psi_star_dot_evol');

run('../particles_equilibrium/build_pre_collapse_XZ_maps.m');
