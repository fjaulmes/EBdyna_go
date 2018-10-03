
for f=1:100
time=((f-1)/(NB_FRAME))*2*pi/omega_TAE

W_TAE=0;
W_TAE_phi=zeros(NB_PHI,1);
options = optimset('TolFun',1e-6);

for phi_rank=1:1
    phi_nTAE=2*pi*(phi_rank-1)/(NB_PHI-1)
    map_TAE_2modes_phi_rank;
    A_map_phi(phi_rank,:,:)=A_PR_map(:,pTAE_inf:pTAE_sup);
    iA_map_phi(phi_rank,:,:)=iA_PR_map(:,pTAE_inf:pTAE_sup);
    
    Epot_map_phi(phi_rank,:,:)=Phi_PR_map(:,pTAE_inf:pTAE_sup);
    iEpot_map_phi(phi_rank,:,:)=iPhi_PR_map(:,pTAE_inf:pTAE_sup);
    
    dA_map_phi(phi_rank,:,:)=dA_PR_map(:,pTAE_inf:pTAE_sup);
    dEpot_map_phi(phi_rank,:,:)=dPhi_PR_map(:,pTAE_inf:pTAE_sup);
    


   
end
    imagesc(Phi_PR_map');
    pause(0.2)
end
Phi_PR_map_100=Phi_PR_map;
A_PR_map_100=A_PR_map;

% W_TAE=sum(W_TAE_phi(1:end-1))
% W_TAE_evol(f)=W_TAE;
    
% Bstar_map_phi=sqrt(1/W_TAE)*Bstar_map_phi;
% BsX_map_phi=sqrt(1/W_TAE)*BsX_map_phi;
% BsZ_map_phi=sqrt(1/W_TAE)*BsZ_map_phi;
% Bsphi_map_phi=sqrt(1/W_TAE)*Bsphi_map_phi;
% psi_map_phi=sqrt(1/W_TAE)*psi_map_phi;
% iA_map_phi=sqrt(1/W_TAE)*iA_map_phi;
% 
% Epot_map_phi=sqrt(1/W_TAE)*Epot_map_phi;
% iEpot_map_phi=sqrt(1/W_TAE)*iEpot_map_phi;
% Efield_X_map_phi=sqrt(1/W_TAE)*Efield_X_map_phi;
% Efield_Z_map_phi=sqrt(1/W_TAE)*Efield_Z_map_phi;
