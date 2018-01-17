
% f=80
% time=2*((f-1)/(NB_FRAME))*2*pi/omega_TAE

W_TAE=0;
W_TAE_phi=zeros(NB_PHI,1);
Wmag_TAE_phi=zeros(NB_PHI,1);
WP_TAE_phi=zeros(NB_PHI,1);
Wpot_TAE_phi=zeros(NB_PHI,1);
WK_TAE_phi=zeros(NB_PHI,1);

load(filenameB);
load(filenameE);

for phi_rank=1:NB_PHI-1
    phi_nTAE=2*pi*(phi_rank-1)/(NB_PHI-1)
    map_WTAE_phi_rank;
    Wmag_TAE_phi(phi_rank)=nTAE*0.5*DPHI*DX*DX*sum(sum(Rpos_XZsmall_map.*...
        ((BsX_XZ_map.^2+BsZ_XZ_map.^2+Bsphi_XZ_map.^2)/mu0)));%+...

    WK_TAE_phi(phi_rank)=nTAE*0.5*DPHI*DX*DX*sum(sum(Rpos_XZsmall_map.*...
        (mH*Ni0).*(ksi_dot_XZ_map_frame_rank).^2));
    WP_TAE_phi(phi_rank)=(WK_TAE_phi(phi_rank)-Wmag_TAE_phi(phi_rank));

    W_TAE_phi(phi_rank)=Wmag_TAE_phi(phi_rank)+WK_TAE_phi(phi_rank);
  
    Wpot_TAE_phi(phi_rank)=W_TAE_phi(phi_rank)-WK_TAE_phi(phi_rank);


end

W_TAE=sum(WK_TAE_phi(1:end-1))
W_TAE_evol(f)=W_TAE;
    

