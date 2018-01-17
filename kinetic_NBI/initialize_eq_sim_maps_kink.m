
tau_sim=time_scale(end)
%TPRECISE=3.0;
%NB_TIME_STAMPS=length(time_scale)
%DELTA_TIME=(1e-9)/TPRECISE;
%h=DELTA_TIME;
%NB_TIME_STEPS=round(0.1*tau_sim/DELTA_TIME)
%RECORD_PRECISION=round(NB_TIME_STEPS/NB_TIME_STAMPS)

% no more than one time stamp in one loop
% ten time steps in one loop
%time_scale=(1:NB_TIME_STAMPS)*DELTA_TIME*RECORD_PRECISION*10;

%tau_sim=2.2e-4;

TPRECISE=1.0;
TIME_STAMP_PRECISION=10;
RECORD_PRECISION=TIME_STAMP_PRECISION*20;
DELTA_TIME=(1e-9)/TPRECISE;
h=DELTA_TIME;
NB_TIME_STEPS=round(0.1*tau_sim/DELTA_TIME)
% no more than one time stamp in one loop
NB_TIME_STAMPS=round(NB_TIME_STEPS/RECORD_PRECISION)
% ten time steps in one loop
time_scale=(1:NB_TIME_STAMPS)*DELTA_TIME*RECORD_PRECISION*10;



gamma_vD_omega_psi_map=zeros(NB_BINS_OMEGA,NB_BINS_PSI);
gamma_vD_theta_psi_map=zeros(NB_BINS_THETA,NB_BINS_PSI);
gamma_vD_phi_psi_map=zeros(129,NB_BINS_PSI);
gamma_mu_omega_psi_map=zeros(NB_BINS_OMEGA,NB_BINS_PSI);
gamma_mu_theta_psi_map=zeros(NB_BINS_THETA,NB_BINS_PSI);

%run('calculate_collapse_frozen_drift_speed_maps')

    

    
theta_XZsmall_map=zeros(size_X,size_Z);
theta_XZsmall_map(:,:)=theta_XZ_map(Xinf:Xsup,Zinf:Zsup);
    








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

SIMULATION_RADIAL_LIMIT=Nradial-6

