
if APPLY_RMP_FIELD==1
    
    load('BS_AUG_theta_psi_phi_2016-03-08.mat');
    BR_RMP=BR;
    BZ_RMP=BZ;
    Bphi_RMP=Bphi;
    if CALCULATE_TRUE_PPHI==1
        dAphi_dphi=permute(dAphi_dphi,[3 1 2]);
    end
    
    BR_RMP=permute(BR_RMP,[3 1 2]);
    BZ_RMP=permute(BZ_RMP,[3 1 2]);
    Bphi_RMP=permute(Bphi_RMP,[3 1 2]);
    
    
    NB_PHI=size(BR_RMP,1)       % Number of toroidal positions incl. 0 and 2pi
    NB_THETA_RMP=size(BR_RMP,2) % Number of poloidal positions incl. 0 and 2pi
    NB_PSI=size(BR_RMP,3)       % Number of radial positions
    DPSIVAL=round(Nradial/NB_PSI)  % estimating radial interval for 3D fields
    
    DTHETA=2*pi/(NB_THETA_RMP-1)
    DPHI=2*pi/(NB_PHI-1);   %Angle between toroidal positions
    
    scale_phi=2*pi*((0:(NB_PHI-1))/(NB_PHI-1));
    scale_theta=2*pi*((0:(NB_THETA_RMP-1))/(NB_THETA_RMP-1));
    scale_psi=1:DPSIVAL:Nradial;
    NB_PHI_DATA_HALF=round(0.5*(NB_PHI-1)); %Number of half the toroidal positions (not used here)
    
    clear BR BZ Bphi;
    
end
    
filename=strcat(DATA_FOLDER,'tokamak_PR_map.mat');
load(filename,'Nradial');
filename=strcat(DATA_FOLDER,'B_fields.mat');
load(filename,'Btor_PR_map');

Bavg=mean(Btor_PR_map(:,1))

% corrected theta maps for complete interpolation
QNB_THETA=round(0.25*NB_THETA);
HQNB_THETA=round(0.5*QNB_THETA);
theta_low_XZsmall_map=theta_XZsmall_map;
theta_up_XZsmall_map=theta_XZsmall_map;
theta_up_XZsmall_map(find(theta_XZsmall_map<QNB_THETA*DTHETA))=theta_XZsmall_map(find(theta_XZsmall_map<QNB_THETA*DTHETA))+2*pi;
theta_low_XZsmall_map(find(theta_XZsmall_map>(NB_THETA-QNB_THETA-2)*DTHETA))=theta_XZsmall_map(find(theta_XZsmall_map>(NB_THETA-QNB_THETA-2)*DTHETA))-2*pi;

% for correction of Eperp
BMAX=max(max(Btot_XZ_map))
BMIN=min(min(Btot_XZ_map(Btot_XZ_map>1)))

SIMULATION_RADIAL_LIMIT=Nradial-4
