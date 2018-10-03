
% f=80
% time=2*((f-1)/(NB_FRAME))*2*pi/omega_TAE

W_TAE=0;
W_TAE_phi=zeros(NB_PHI,1);
WK_TAE_phi=zeros(NB_PHI,1);
WP_TAE_phi=zeros(NB_PHI,1);
Wpot_TAE_phi=zeros(NB_PHI,1);
options = optimset('TolFun',1e-6);

for phi_rank=1:NB_PHI
    phi_nTAE=2*pi*(phi_rank-1)/(NB_PHI-1)
    map_TAE_2modes_EBfields_phi_rank;
    Wmag_TAE_phi(phi_rank)=nTAE*0.5*DPHI*DX*DX*sum(sum(Rpos_XZsmall_map.*...
        (Bstar_XZ_map_frame_rank.^2)/mu0));%+...
        %epsilon0*(EX_XZ_map.^2+EZ_XZ_map.^2+Ephi_XZ_map.^2)-...
        %(ksi_X.*J0_B1_X+ksi_Z.*J0_B1_Z+ksi_phi.*J0_B1_phi))));
%     WP_TAE_phi(phi_rank)=nTAE*0.5*DPHI*DX*DX*sum(sum(Rpos_XZsmall_map.*...
%         (-P1.*div_ksi)));
    WK_TAE_phi(phi_rank)=nTAE*0.5*DPHI*DX*DX*sum(sum(Rpos_XZsmall_map.*...
        (mH*Ni0).*(ksi_dot_XZ_map_frame_rank).^2));
    WP_TAE_phi(phi_rank)=(WK_TAE_phi(phi_rank)-Wmag_TAE_phi(phi_rank));
%     ksi_amplitude=sqrt((WK_TAE_phi(phi_rank)-Wmag_TAE_phi(phi_rank))/WP_TAE_phi(phi_rank));

    W_TAE_phi(phi_rank)=Wmag_TAE_phi(phi_rank)+WK_TAE_phi(phi_rank);
%     WK_TAE_phi(phi_rank)=WK_TAE_phi(phi_rank)/omega_TAE^2;
%     Wpot_TAE_phi(phi_rank)=nTAE*0.5*DPHI*DX*DX*sum(sum(Rpos_XZsmall_map.*...
%         ((BstarX_XZ_map.^2+BstarZ_XZ_map.^2+Bstarphi_XZ_map.^2)/mu0)));
  
    Wpot_TAE_phi(phi_rank)=W_TAE_phi(phi_rank)-WK_TAE_phi(phi_rank);

    
    %     Bstar_map_phi(phi_rank,:,:)=Btot_PR_map(:,pTAE_inf:pTAE_sup);
    BsX_map_phi(phi_rank,:,:)=bX_PR_map(:,pTAE_inf:pTAE_sup);
    BsZ_map_phi(phi_rank,:,:)=bZ_PR_map(:,pTAE_inf:pTAE_sup);
    Bsphi_map_phi(phi_rank,:,:)=bphi_PR_map(:,pTAE_inf:pTAE_sup);
%    psi_map_phi(phi_rank,:,:)=psi_star_PR_map(:,pTAE_inf:pTAE_sup);
    A_map_phi(phi_rank,:,:)=A_PR_map(:,pTAE_inf:pTAE_sup);
    iA_map_phi(phi_rank,:,:)=iA_PR_map(:,pTAE_inf:pTAE_sup);
    
    Epot_map_phi(phi_rank,:,:)=Phi_PR_map(:,pTAE_inf:pTAE_sup);
    iEpot_map_phi(phi_rank,:,:)=iPhi_PR_map(:,pTAE_inf:pTAE_sup);
    Efield_X_map_phi(phi_rank,:,:)=EX_PR_map(:,pTAE_inf:pTAE_sup);
    Efield_Z_map_phi(phi_rank,:,:)=EZ_PR_map(:,pTAE_inf:pTAE_sup);
    
    dA_map_phi(phi_rank,:,:)=dA_PR_map(:,pTAE_inf:pTAE_sup);
    dEpot_map_phi(phi_rank,:,:)=dPhi_PR_map(:,pTAE_inf:pTAE_sup);
   
end

W_TAE=sum(W_TAE_phi(1:end-1))
W_TAE_evol(f)=W_TAE;
    
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
