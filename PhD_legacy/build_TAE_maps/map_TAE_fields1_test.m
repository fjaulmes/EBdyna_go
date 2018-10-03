close all
clear all
initialize_folder_names
initialize_TAE_map_calculation_context;
TAE_parameters;
initialize_XZ_maps_dimensions;
rescaling_to_XZsmall_maps;

%%
ion_density_XZ_map=interp1(psi_scale,Ne_profile,psi_XZsmall_map);
ion_density_XZ_map(isnan(ion_density_XZ_map))=0;
P0_PR_map=zeros(NP,Nradial);

for p=1:NP
    for r=1:Nradial
        P0_PR_map(p,r)=P_initial_profile(r);
    end
end

if EVALUATE_PRESSURE==1
    
    P_data=reshape((P0_PR_map(:,1:pTAE_sup)),NP*pTAE_sup,1);
    P0_XZ_map=griddata(finesse_data_X,finesse_data_Z,P_data,XX_small,ZZ_small);
    P0_XZ_map(isnan(P0_XZ_map)) = 0;
    P0_XZ_map=P0_XZ_map';
    
    
    gP0_X=zeros(sizeX,sizeZ);
    gP0_Z=zeros(sizeX,sizeZ);
    
    for (x=3:sizeX-2)
        for (z=3:sizeZ-2)
            gP0_X(x,z)=(1/12)*(-P0_XZ_map(x+2,z)+P0_XZ_map(x-2,z))+(2/3)*(P0_XZ_map(x+1,z)-P0_XZ_map(x-1,z));
            gP0_Z(x,z)=(1/12)*(-P0_XZ_map(x,z+2)+P0_XZ_map(x,z-2))+(2/3)*(P0_XZ_map(x,z+1)-P0_XZ_map(x,z-1));
            %         gpsi_Z(x,z)=0.5*(psi2D(x,z+1)-psi2D(x,z-1));
        end
    end
    gP0_X=gP0_X/DX;
    gP0_Z=gP0_Z/DX;
    
end

bphi_XZsmall_map=sqrt(1-bX_XZ_small_map.^2+bZ_XZ_small_map.^2);

%%
size_r_TAE=pTAE_sup-pTAE_inf+1

clc
NB_FRAME=100
Efield_X_map_phi=zeros(NB_PHI,NP,size_r_TAE);
Efield_Z_map_phi=zeros(NB_PHI,NP,size_r_TAE);
iEpot_map_phi=zeros(NB_PHI,NP,size_r_TAE);
Epot_map_phi=zeros(NB_PHI,NP,size_r_TAE);
% Bstar_map_phi=zeros(NB_PHI,NP,size_r_TAE);
BsX_map_phi=zeros(NB_PHI,NP,size_r_TAE);
BsZ_map_phi=zeros(NB_PHI,NP,size_r_TAE);
Bsphi_map_phi=zeros(NB_PHI,NP,size_r_TAE);
%psi_map_phi=zeros(NB_PHI,NP,size_r_TAE);
iA_map_phi=zeros(NB_PHI,NP,size_r_TAE);
A_map_phi=zeros(NB_PHI,NP,size_r_TAE);

%
% build n=18 ; m=19 and m=20 TAE
Phi_PR_map=zeros(NP,Nradial);
A_PR_map=zeros(NP,Nradial);
EX_PR_map=zeros(NP,Nradial);
EZ_PR_map=zeros(NP,Nradial);
omega_PR_map=zeros(NP,Nradial);
psi_star_PR_map=zeros(NP,Nradial);

theta=2*pi*((1:NP)-1)/(NP-1);


vA3_TAE=vA_TAE/3


filename='../data_tokamak/TAE_data.mat';

load(filename);
% disp('****** SAVING TAE DATA in ../data_tokamak/  ****');
% save(filename,'kTAE','omega_TAE','nTAE','qTAE','pTAE_inf','pTAE_sup','size_r_TAE','TAE_WIDTH','MODE_WIDTH','NB_FRAME','NB_PHI','TAE_angle',...
%     'mtheta1','mtheta2','r1','r2','rTAE','coupling_m_mp1',...
%     'flc_s','Phi0','vA_TAE','vA3_TAE','psiTAE','pphi_TAE','OPP_W');
% nTAE


%%
W_TAE_evol=zeros(NB_FRAME,1);
UNIFORM_TAE_PROPAGATION=1

for f=1:1
    f
    % time should loop to right the begining of another TAE period
    time=((f-1)/(NB_FRAME))*2*pi/omega_TAE
    
    %assuming propagation at average Alfven velocity
    OMEGA_VA=omega_TAE;
    
    map_TAE_2modes_EBfields_frame_rank;
    
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
    %save(filename,'BsX_map_phi','BsZ_map_phi','Bsphi_map_phi','psi_map_phi','iA_map_phi');
%     save(filename,'BsX_map_phi','BsZ_map_phi','Bsphi_map_phi','A_map_phi','iA_map_phi');
    save(filename,'BsX_map_phi','BsZ_map_phi','Bsphi_map_phi');
    disp('****** SAVING FILE in ../B_maps/  ****');
    
    filename='../E_maps/E0';
    filename=strcat(filename,frame_name,'.mat');   
%     save(filename,'Efield_X_map_phi','Efield_Z_map_phi','Epot_map_phi','iEpot_map_phi');
    save(filename,'Efield_X_map_phi','Efield_Z_map_phi');
    
    disp('****** SAVING FILE in ../E_maps/  ****');

    
end

% load('Epot_map_data.mat')

Bavg=mean(mean(Btot_PR_map(:,pTAE_inf:pTAE_sup)))
Bs_amplitude=sqrt(BsX_map_phi.^2+BsZ_map_phi.^2+Bsphi_map_phi.^2);
deltaB=max(max(max(Bs_amplitude)))
deltaB_B0=deltaB/Bavg

W_TAE_oscill_evol=W_TAE_evol;

%W_TAE_oscill_evol=W_TAE_evol/W_TAE_oscill_evol(1)
WTAE_AVG=mean(W_TAE_oscill_evol(1:1))


% correction needed for radial profiles ....

kpllm1_profile(1:pTAE_inf)=kpllm1_profile(pTAE_inf);
kpllm1_profile(pTAE_sup:Nradial)=kpllm1_profile(pTAE_sup);

kpllm2_profile(1:pTAE_inf)=kpllm2_profile(pTAE_inf);
kpllm2_profile(pTAE_sup:Nradial)=kpllm2_profile(pTAE_sup);

m_adjust_profile1(1:pTAE_inf)=m_adjust_profile1(pTAE_inf);
m_adjust_profile1(pTAE_sup:Nradial)=m_adjust_profile1(pTAE_sup);

m_adjust_profile2(1:pTAE_inf)=m_adjust_profile2(pTAE_inf);
m_adjust_profile2(pTAE_sup:Nradial)=m_adjust_profile2(pTAE_sup);

% filename='../data_tokamak/TAE_data.mat';
% save(filename,'-append','WTAE_AVG','deltaB','deltaB_B0','kpllm1_profile','kpllm2_profile','m_adjust_profile1','m_adjust_profile2');

% filename='./W_TAE_oscill_evol1.mat';
% save(filename,'W_TAE_oscill_evol');

% plot(W_TAE_oscill_evol)


exit

