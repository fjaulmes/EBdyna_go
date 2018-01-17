
% close all;
pause(0.2);

REINIT_ALL_MAPS=0;
REINIT_SIMULATION_DATA=1;

if REINIT_ALL_MAPS==1
    clear all
    initialize_folder_names;
    filename=strcat(DATA_FOLDER,'physics_constants.mat');
    load(filename);
    filename=strcat(DATA_FOLDER,'tokamak_map_dimensions.mat');
    load(filename);

    mid_X_large=mid_X;
    mid_Z_large=mid_Z;    
    NB_PHI=129;
    DPHI=2*pi/(NB_PHI-1);
    NB_PHI_DATA_HALF=round(0.5*(NB_PHI-1));
    calculate_pre_collapse_drift_speed_maps;
    clear all;
    REINIT_SIMULATION_DATA=1;
end

if REINIT_SIMULATION_DATA==1
    clear all;
    clc
    format compact
    initialize_folder_names_test;
    
    filename=strcat(DATA_FOLDER,'physics_constants.mat');
    load(filename);
    filename=strcat(DATA_FOLDER,'XZsmall_fields_tokamak_pre_collapse.mat');
    load(filename);
    filename=strcat(DATA_FOLDER,'motions_map_dimensions.mat');
    load(filename);
    filename=strcat(DATA_FOLDER,'psi_star_evol.mat');
    load(filename,'psi_star_2D_evol_lin');
    % clear psi_star_2D_evol_interp;  % no need for this high precision

    filename=strcat(DATA_FOLDER,'flux_geometry.mat');
    load(filename, 'dr_X_PR_map');
    load(filename, 'dr_Z_PR_map');

    filename=strcat(DATA_FOLDER,'Epot_psi_star_dot_evol.mat');
    load(filename,'Epot_evol');
    filename=strcat(DATA_FOLDER,'psi_profiles_kadomtsev.mat');
    load(filename, 'psi_star_initial');
    
    REINIT_SIMULATION_DATA=1;
   
end

mHe=mD
ZHe=1

TOROIDAL_FIELD_DIRECTION=sign(mean(mean(Bphi_XZsmall_map)))
PSI_CORE_SIGN=sign(psi_scale(1))
PSI_STAR_SIGN=sign(psi_star_initial(round(0.5*size_r)))

% psi_star_map_phi_evol=interp1(1:1001,psi_star_2D_evol_interp,1:101,'cubic');

% initialize completely the maps for pre-collapse calculations
FAST_SAWTOOTH=2;
tau_cr=4e-4;

TPRECISE=2;
TIME_STAMP_PRECISION=20*TPRECISE;
DELTA_TIME=(1e-9)/TPRECISE;
h=DELTA_TIME;
NB_TIME_STEPS=round(0.05*(tau_cr/FAST_SAWTOOTH)/DELTA_TIME)
% one time stamp in one loop
NB_TIME_STAMPS=round(NB_TIME_STEPS/TIME_STAMP_PRECISION);
% ten time steps in one loop
time_scale=(1:NB_TIME_STAMPS)*DELTA_TIME*TIME_STAMP_PRECISION*20;

%simulation options
SAVE_DATA_FILE=1;
DISPLAY_OUTPUTS=1;
CALCULATE_TRUE_PPHI=1;
EKIN0_FAC=1
SAVENAME=strcat('fewG_Ekin',num2str(EKIN0_FAC*8));
SAVENAME=strcat(SAVENAME,'_collapse_041214_fc2h');
SAVENAME=strcat(SAVENAME,num2str(TPRECISE))


if SAVE_DATA_FILE==0
    disp('no data file for this simulation')
end

if DISPLAY_OUTPUTS==0
    disp('no graphics will be displayed during this simulation')
else
    filename=strcat(DATA_FOLDER,'tokamak_PR_map.mat');
    load(filename);
    [XXsmall ZZsmall]=meshgrid(scale_X,scale_Z);
    finesse_data_X=reshape((Rpos_PR_map(:,1:size_r)-R0),NP*size_r,1);
    finesse_data_Z=reshape(Z_PR_map(:,1:size_r),NP*size_r,1);
end


%run('calculate_collapse_frozen_drift_speed_maps')
% psi_star_map_phi_evol=interp1(1:1001,psi_star_map_phi_evol_interp,1:101,'cubic');
% run('calculate_rotB_vDcurv');
% rotational_b_pll=(rotB_phi.*bphi_XZ_map+rotB_Z.*bZ_XZ_map+rotB_X.*bX_XZ_map)./Btot_XZ_map;

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


Nalphas_psi_pos=28;
Nalphas_lambda=Nalphas_psi_pos*16;
Nalphas_simulated=Nalphas_lambda*8;

alphas_pos_z=0.6*[-0.53 -0.51 -0.49 -0.46 -0.42 -0.37 -0.33 -0.3 -0.28 -0.24 -0.18 -0.12 -0.08 -0.02 0.02 0.08 0.12 0.18 0.24 0.28 0.3 0.33 0.37 0.42 0.46 0.49 0.51 0.53]
alphas_pos_z=[alphas_pos_z alphas_pos_z alphas_pos_z alphas_pos_z];
alphas_pos_z=[alphas_pos_z alphas_pos_z  ];
alphas_pos_x=alphas_pos_z*0-0.01;

alphas_pos_x=[alphas_pos_x alphas_pos_x];
alphas_pos_z=[alphas_pos_z alphas_pos_z];

alphas_pos_x=[alphas_pos_x alphas_pos_x alphas_pos_x alphas_pos_x alphas_pos_x alphas_pos_x alphas_pos_x alphas_pos_x];
alphas_pos_z=[alphas_pos_z alphas_pos_z alphas_pos_z alphas_pos_z alphas_pos_z alphas_pos_z alphas_pos_z alphas_pos_z];

alphas_pos_z=alphas_pos_z' + Z_axis;
alphas_pos_x=alphas_pos_x' + Raxis-R0;

alphas_pos_phi=ones(1,Nalphas_psi_pos)*0.39;
alphas_pos_phi=[alphas_pos_phi 2*alphas_pos_phi 3*alphas_pos_phi 4*alphas_pos_phi];
alphas_pos_phi=[alphas_pos_phi alphas_pos_phi+0.39*4 alphas_pos_phi+0.39*8 alphas_pos_phi+0.39*12];

alphas_pos_phi=[alphas_pos_phi alphas_pos_phi alphas_pos_phi alphas_pos_phi alphas_pos_phi alphas_pos_phi alphas_pos_phi alphas_pos_phi];
alphas_pos_phi=alphas_pos_phi';

alphas_Bfield=interp2(scale_X,scale_Z,Btot_XZ_map',alphas_pos_x,alphas_pos_z,'*linear');

alphas_Ekin=EKIN0_FAC*ones(1,Nalphas_psi_pos)*8*1e3;  % keV
alphas_Ekin=[alphas_Ekin alphas_Ekin alphas_Ekin alphas_Ekin];
alphas_Ekin=[alphas_Ekin alphas_Ekin alphas_Ekin alphas_Ekin];
alphas_Ekin=alphas_Ekin';

alphas_vpll=(ones(Nalphas_lambda,1)*0.1);
alphas_vpll=[-9*alphas_vpll.*sqrt(alphas_Ekin*2*eV/mHe) ; -6*alphas_vpll.*sqrt(alphas_Ekin*2*eV/mHe) ; -1*alphas_vpll.*sqrt(alphas_Ekin*2*eV/mHe) ; -0.5*alphas_vpll.*sqrt(alphas_Ekin*2*eV/mHe) ; ...
 0.5*alphas_vpll.*sqrt(alphas_Ekin*2*eV/mHe) ; 1*alphas_vpll.*sqrt(alphas_Ekin*2*eV/mHe) ; 6*alphas_vpll.*sqrt(alphas_Ekin*2*eV/mHe) ; 9*alphas_vpll.*sqrt(alphas_Ekin*2*eV/mHe)];

alphas_Ekin=[alphas_Ekin ; alphas_Ekin ; alphas_Ekin ; alphas_Ekin ;  alphas_Ekin ; alphas_Ekin ; alphas_Ekin ; alphas_Ekin];
alphas_mm=(alphas_Ekin-0.5*((mHe/eV)*alphas_vpll.^2))./alphas_Bfield;
alphas_ejected=zeros(Nalphas_simulated,1);
alphas_Eperp=alphas_Bfield.*alphas_mm;

alphas_Epot=zeros(Nalphas_simulated,1);
alphas_momentum=zeros(Nalphas_simulated,1);


Bavg=mean(Btot_XZ_map(round(0.8*mid_X),:));

%initial accurate energy values
alphas_vE_sq=zeros(Nalphas_simulated,1);

alphas_Epll=0.5*(mHe/eV)*alphas_vpll.^2;

%initial position values
alphas_mm0=alphas_mm;
X0=alphas_pos_x;
Z0=alphas_pos_z;
phi0=alphas_pos_phi;


alphas_Eperp=alphas_Ekin-alphas_Epll;
alphas_vperp=sqrt(2*alphas_Eperp*eV/mHe);
alphas_rhoL=(mHe/eV)*(1/ZHe)*alphas_vperp./alphas_Bfield;


disp('**************************************************************');

%initializing simulation arrays
Xpos_output=zeros(NB_TIME_STAMPS,Nalphas_simulated);
Zpos_output=zeros(NB_TIME_STAMPS,Nalphas_simulated);
phipos_output=zeros(NB_TIME_STAMPS,Nalphas_simulated);
psipos_output=zeros(NB_TIME_STAMPS,Nalphas_simulated);
vD_output=zeros(NB_TIME_STAMPS,Nalphas_simulated);
vparallel_output=zeros(NB_TIME_STAMPS,Nalphas_simulated);
Eperp_output=zeros(NB_TIME_STAMPS,Nalphas_simulated);
vEsq_output=zeros(NB_TIME_STAMPS,Nalphas_simulated);
theta_output=zeros(NB_TIME_STAMPS,Nalphas_simulated);
Ekin_output=zeros(NB_TIME_STAMPS,Nalphas_simulated);
vphi_output=zeros(NB_TIME_STAMPS,Nalphas_simulated);
psi_star_output=zeros(NB_TIME_STAMPS,Nalphas_simulated);
pphi_output=zeros(NB_TIME_STAMPS,Nalphas_simulated);
% kin_energy_corr_output=zeros(NB_TIME_STAMPS,Nalphas_simulated);
psi_value_output=zeros(NB_TIME_STAMPS,Nalphas_simulated);
rhoL_output=zeros(NB_TIME_STAMPS,Nalphas_simulated);
Epot_output=zeros(NB_TIME_STAMPS,Nalphas_simulated);
Bfield_output=zeros(NB_TIME_STAMPS,Nalphas_simulated);
vEradial_output=zeros(NB_TIME_STAMPS,Nalphas_simulated);
% pphi0_output=zeros(NB_TIME_STAMPS,Nalphas_simulated);
Etot_output=zeros(NB_TIME_STAMPS,Nalphas_simulated);



alphas_psi=zeros(Nalphas_simulated,1);
alphas_theta=zeros(Nalphas_simulated,1);
vExB_X=zeros(Nalphas_simulated,1);
vExB_Z=zeros(Nalphas_simulated,1);
vExB_phi=zeros(Nalphas_simulated,1);
bX=zeros(Nalphas_simulated,1);
bZ=zeros(Nalphas_simulated,1);
bphi=zeros(Nalphas_simulated,1);
alphas_psi_star=zeros(Nalphas_simulated,1);
alphas_psi_value=zeros(Nalphas_simulated,1);
alphas_psiH_value=zeros(Nalphas_simulated,1);
alphas_psi_value_corr=zeros(Nalphas_simulated,1);
INDEX_LIST_OMEGA=ones(Nalphas_simulated,1);
alphas_grad_Phi=zeros(Nalphas_simulated,1);
% alphas_kin_energy_corr=zeros(Nalphas_simulated,1);
EX=zeros(Nalphas_simulated,1);
EZ=zeros(Nalphas_simulated,1);
Ephi=zeros(Nalphas_simulated,1);
alphas_grad_psi_star=zeros(Nalphas_simulated,1);

% intermediate calculations arrays
delta_E=zeros(Nalphas_simulated,1);
alphas_vpll_int=alphas_vpll;
vparallel_tilde=alphas_vpll;


% initial field and position
alphas_psi_value=interp2(scale_X,scale_Z,psi_XZsmall_map',alphas_pos_x,alphas_pos_z,'*linear');
if PSI_CORE_SIGN<0
    alphas_psi_value=max(alphas_psi_value,-psi_global);
else
    alphas_psi_value=min(alphas_psi_value,psi_global);
end
alphas_psi=interp1(psi_scale,1:257,alphas_psi_value);

alphas_theta=interp2(scale_X,scale_Z,theta_XZsmall_map',alphas_pos_x,alphas_pos_z,'*nearest');


frame_rank_next_prev=0;
time=0;

% Initialize drift values
frame_rank_precise=1;
frame_rank_precise_prev=frame_rank_precise

frame_rank_prev=10*floor((frame_rank_precise-1)/10)+1;
frame_rank_next=10*ceil((frame_rank_precise-1)/10)+1;

frame_rank=frame_rank_prev;

load_Emaps_frame_rank;
load_Bmaps_frame_rank;

% patch_Btot_map_phi;

% psi star map
psi_star_omega_map_half(:,:)=psi_star_2D_evol_lin(round(frame_rank/10)+1,:,:);

% Using symmetry to reconstruct a poloidal turn
psi_star_omega_map=zeros(size_r,NB_THETA);
psi_star_omega_map(:,1:round(0.5*NB_THETA))=psi_star_omega_map_half(:,:);
psi_star_omega_map(:,round(0.5*NB_THETA):NB_THETA)=psi_star_omega_map_half(:,round(0.5*NB_THETA):-1:1);
psi_star_omega_map_rank=psi_star_omega_map;

E_potential_omega_map(:,:)=Epot_evol(round((frame_rank-1)/10+1),:,:);


bX_map_phi_next=bX_map_phi;
bZ_map_phi_next=bZ_map_phi;
% bphi_map_phi_next=bphi_map_phi;
B_PR_map_phi_next=Btot_map_phi;

%     Efield_X_map_phi_next=zeros(NB_PHI,NB_THETA,size_r);
%     Efield_Z_map_phi_next=zeros(NB_PHI,NB_THETA,size_r);
%     Epot_map_phi_next=zeros(NB_PHI,NB_THETA,size_r);


Efield_X_map_phi_next=FAST_SAWTOOTH*Efield_X_map_phi;
Efield_Z_map_phi_next=FAST_SAWTOOTH*Efield_Z_map_phi;
Epot_map_phi_next=FAST_SAWTOOTH*E_potential_omega_map;

psi_star_map_next=psi_star_omega_map_rank;


grad_Phi_map_phi_next=FAST_SAWTOOTH*grad_Phi_tor_map_phi;
% Efield_X_map_phi_next=FAST_SAWTOOTH*Efield_X_map_phi;
% Efield_Z_map_phi_next=FAST_SAWTOOTH*Efield_Z_map_phi;
% Efield_phi_map_phi_next=FAST_SAWTOOTH*Efield_phi_map_phi;
grad_psi_star_map_phi_next=grad_psi_star_map_phi;

    
frame_rank_next=frame_rank_prev;

% estimating psi and theta values
alphas_pos_phi_wrapped=wrap2pi(alphas_pos_phi);
interpolate_theta_psi_fromXZ;
% sorting out which particles are inside mixing region and which are out
INNER_PART=find(alphas_psi<=(size_r-1));
OUTER_PART=find(alphas_psi>(size_r-1));
[IL3D_1 IL3D_2 IL3D_3 IL3D_4 IL3D_5 IL3D_6 IL3D_7 IL3D_8 slopex slopey slopez] = ...
        build_3Dinterp_indexarrays(scale_phi, scale_theta, scale_psi, DPHI,DTHETA,1,alphas_pos_phi_wrapped(INNER_PART), alphas_theta(INNER_PART), alphas_psi(INNER_PART));

if (DISPLAY_OUTPUTS==1)
    PN=12+Nalphas_lambda*1
    phi_rank=round(alphas_pos_phi_wrapped(PN)/DPHI)+1;
    Xpos=alphas_pos_x(PN)
    Zpos=alphas_pos_z(PN)
    Xpos_prev=Xpos;
    Zpos_prev=Zpos;
end

delta_t_coef=100*(h*FAST_SAWTOOTH/tau_cr);
update_E_b_fields;


% direction of the B field at local positions
% amplitude of the B field and potential at local positions
update_Gfields_collapse;

norm_angle=rand(Nalphas_simulated,1)*2*pi;
uX=sqrt(1./(1+(bX./bphi).^2));
uZ=0*uX;
uphi=-(bX./bphi).*uX;
unorm=sqrt(uX.^2+uZ.^2+uphi.^2);

wX=1./sqrt((uX./uphi).^2+1+((1./bZ).*(uX./uphi).*(bphi-bX)).^2);
wZ=(wX./bZ).*((uX./uphi).*bphi-bX);
wphi=-(uX./uphi).*wX;
wnorm=sqrt(wX.^2+wZ.^2+wphi.^2);
wX=wX./wnorm;
wZ=wZ./wnorm;
wphi=wphi./wnorm;

% normal vector N = (cos) u + (sin) w
% that N.b=0 precisely
NX=(cos(norm_angle).*uX+sin(norm_angle).*wX);
NZ=((cos(norm_angle).*uZ+sin(norm_angle).*wZ));
Nphi=(cos(norm_angle).*uphi+sin(norm_angle).*wphi);


v0_X=alphas_vpll.*bX+alphas_vperp.*NX;
v0_Z=alphas_vpll.*bZ+alphas_vperp.*NZ;
v0_phi=alphas_vpll.*bphi+alphas_vperp.*Nphi;


v_X=v0_X;
v_Z=v0_Z;
v_phi=v0_phi;
alphas_Ekin0=alphas_Ekin;


%initialize Bfield properly
alphas_mm_part=alphas_mm;
% update_GT_3D_collapse;

interpolate_theta_psi_fromXZ;
INNER_PART=find(alphas_psi<=(size_r-1));
OUTER_PART=find(alphas_psi>(size_r-1));


alphas_Omega=(ZHe*eV/mHe)*alphas_Bfield;
epsilon=0.5*h*alphas_Omega;
alphas_determinant=1+(epsilon.^2);

% previous mid point speed values
v_X_prev=v0_X-epsilon.*(bphi.*v0_Z-bZ.*v0_phi+EX);
v_Z_prev=v0_Z-epsilon.*(bX.*v0_phi-bphi.*v0_X+EZ);
v_phi_prev=v0_phi-epsilon.*(bZ.*v0_X-bX.*v0_Z+Ephi);

% since the above is pretty approximate, we rescale for total energy
alphas_vtot_sq=(2*alphas_Ekin*eV/mHe);
vtot_recalc_sq=(v_X.^2+v_Z.^2+v_phi.^2);%
Ekin_corr=0.5*(1+sqrt(alphas_vtot_sq./vtot_recalc_sq));
v_X_prev=v_X_prev.*(Ekin_corr);
v_Z_prev=v_Z_prev.*(Ekin_corr);
v_phi_prev=v_phi_prev.*(Ekin_corr);

v_phi_prev_prev=v_phi_prev;
% v_phi_prev=v_phi;

alphas_omega=wrap2pi(alphas_theta-alphas_pos_phi_wrapped);


if ~isempty(OUTER_PART)
    alphas_grad_Phi(OUTER_PART)=zeros(size(OUTER_PART,1),1);
    alphas_Epot(OUTER_PART)=zeros(size(OUTER_PART,1),1);
end
% alphas_grad_Phi(INNER_PART)=lininterp3( grad_Phi_tor_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
alphas_Epot(INNER_PART)=interp2_omega_map(scale_psi(1:size_r),scale_theta,psi_star_omega_map,1,DTHETA,alphas_psi, alphas_omega,INNER_PART);

% Canonical angular momentum initial value
alphas_psi_value=interp2_XZ(interp_x,interp_z,psi_XZsmall_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
alphas_psiH_value=interp2_XZ(interp_x,interp_z,psiH_XZsmall_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);

alphas_psi_star(INNER_PART)=interp2(scale_psi(1:size_r),scale_theta,psi_star_omega_map',alphas_psi(INNER_PART), alphas_omega(INNER_PART),'*linear');
alphas_psi_value_corr(INNER_PART)=alphas_psiH_value(INNER_PART)+alphas_psi_star(INNER_PART);
alphas_psi_value_corr(OUTER_PART)=alphas_psi_value(OUTER_PART);

alphas_pphi0=(mHe/eV)*alphas_Rpos.*v_phi-(ZHe)*alphas_psi_value_corr;
% alphas_pphi0_half=alphas_pphi0;
% 
% alphas_pphi0_half(INNER_PART)=alphas_pphi0_half(INNER_PART)-(0.5*h*ZHe).*(alphas_Rpos(INNER_PART)).*(alphas_grad_Phi(INNER_PART));



% for eneergy evolution
alphas_Epot_prev=alphas_Epot;
alphas_psi_star_prev=alphas_psi_star;
alphas_Bfield_prev=alphas_Bfield;


%trapping parameter
radial_pos=(a/257)*interp2(scale_X,scale_Z,radial_XZsmall_map',alphas_pos_x,alphas_pos_z,'*linear');
Bavg=mean(alphas_Bfield);
alphas_kappa=sqrt((alphas_Ekin*R0+Bavg*alphas_mm.*(radial_pos-R0))./(2*alphas_mm.*radial_pos*Bavg));
alphas_lambda=Bavg*alphas_mm./alphas_Ekin;

alphas_Eperp=alphas_Ekin-alphas_Epll;
alphas_vperp=sqrt(2*alphas_Eperp*eV/mHe);
alphas_rhoL=(mHe/eV)*(1/ZHe)*alphas_vperp./alphas_Bfield;


alphas_Etot=alphas_Ekin0+ZHe*alphas_Epot;

alphas_Ekin_half=alphas_Ekin;
alphas_Ekin_prev=alphas_Ekin;

adapt_speed_Ekin_G;

interpolate_theta_psi_fromXZ;
INNER_PART=find(alphas_psi<=(size_r-1));
OUTER_PART=find(alphas_psi>(size_r-1));
[IL3D_1 IL3D_2 IL3D_3 IL3D_4 IL3D_5 IL3D_6 IL3D_7 IL3D_8 slopex slopey slopez] = ...
    build_3Dinterp_indexarrays(scale_phi, scale_theta, scale_psi, DPHI,DTHETA,1,alphas_pos_phi_wrapped(INNER_PART), alphas_theta(INNER_PART), alphas_psi(INNER_PART));

% amplitude of the B field and potential gradients at local positions
if ~isempty(OUTER_PART)
    alphas_Bfield(OUTER_PART)=interp2_XZ(interp_x,interp_z,Btot_XZ_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4,OUTER_PART);
    alphas_grad_Phi(OUTER_PART)=zeros(size(OUTER_PART,1),1);
    alphas_grad_psi_star(OUTER_PART)=zeros(size(OUTER_PART,1),1);
    alphas_Epot(OUTER_PART)=zeros(size(OUTER_PART,1),1);
    alphas_psi_star(OUTER_PART)=zeros(size(OUTER_PART,1),1);
end
alphas_Bfield(INNER_PART)=lininterp3( Btot_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
%alphas_grad_Phi(INNER_PART)=lininterp3( grad_Phi_tor_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
alphas_grad_psi_star(INNER_PART)=lininterp3( grad_psi_star_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
alphas_Epot(INNER_PART)=interp2_omega_map(scale_psi(1:size_r),scale_theta,psi_star_omega_map,1,DTHETA,alphas_psi, alphas_omega,INNER_PART);
if PSI_STAR_SIGN>0
    alphas_psi_star(INNER_PART)=max(interp2(scale_psi(1:size_r),scale_theta,psi_star_omega_map',alphas_psi(INNER_PART), alphas_omega(INNER_PART),'*linear'),0);
else
    alphas_psi_star(INNER_PART)=min(interp2(scale_psi(1:size_r),scale_theta,psi_star_omega_map',alphas_psi(INNER_PART), alphas_omega(INNER_PART),'*linear'),0);
end
% alphas_Epot_step_prev=alphas_Epot;
alphas_psi_star_half_prev=alphas_psi_star;
alphas_psi_star_part_prev=alphas_psi_star;
alphas_psi_star_prev=alphas_psi_star;
alphas_Ekin_prev=alphas_Ekin;
delta_ps=alphas_psi_star-alphas_psi_star_prev;
delta_ps_prev=delta_ps;

psi_star_omega_map_half_step=psi_star_omega_map;
%delta_ps_prev=delta_ps;

alphas_omega0=alphas_omega;
alphas_Rpos_int=alphas_Rpos;
alphas_vphi_grad_psi_star=v_phi.*alphas_grad_psi_star;
% half_delta_t_coef=100*(0.5*h*FAST_SAWTOOTH/tau_cr);

load(strcat(DATA_FOLDER,'q_profile.mat'),'q_initial_profile','psi_rank_q1');
alphas_r=interp1(1:Nradial,radial_r_value_flux,alphas_psi);
delta_psi_q1=40
psi_core=round(psi_rank_q1-delta_psi_q1)
r_q1=interp1(1:257,radial_r_value_flux,psi_rank_q1);
r_core=interp1(1:257,radial_r_value_flux,psi_core);

% psiq1=interp1(q_initial_profile,1:Nradial,1);
RED_PART=find(alphas_r<=r_core);
BLUE_PART=find(alphas_r>r_core);

% benchmarking time between two data recordings
disp('**************************************************************');
tic
Xpos_int=alphas_pos_x(PN);
Zpos_int=alphas_pos_z(PN);

time=0;

%%
for time_step=1:NB_TIME_STEPS
   
    time_step_integration_GT_collapse;
    time_step_integration_GT_collapse;
    time_step_integration_GT_collapse;
    time_step_integration_GT_collapse;
    time_step_integration_GT_collapse;
    time_step_integration_GT_collapse;
    time_step_integration_GT_collapse;
    time_step_integration_GT_collapse;
    time_step_integration_GT_collapse;
    time_step_integration_GT_collapse;

	time_step_integration_GT_collapse;
    time_step_integration_GT_collapse;
    time_step_integration_GT_collapse;
    time_step_integration_GT_collapse;
    time_step_integration_GT_collapse;
    time_step_integration_GT_collapse;
    time_step_integration_GT_collapse;
    time_step_integration_GT_collapse;
    time_step_integration_GT_collapse;
    time_step_integration_GT_collapse;

    if (mod(time_step-1,TIME_STAMP_PRECISION)==0)
        time
%         figure(1)
%         subplot(2,1,1)
%         imagesc(E_potential_omega_map'); colorbar
%         
%         subplot(2,1,2)
%         hist(alphas_Epot(INNER_PART),40)
%         
%         pause(0.1)
%         
        adapt_speed_Ekin_G;
        v_X_prev=v_X;
        v_Z_prev=v_Z;
        v_phi_prev=v_phi;
        alphas_vtot_sq=(2*alphas_Ekin*eV/mHe);
        vtot_recalc_sq=(v_X.^2+v_Z.^2+v_phi.^2);%
        Ekin_corr=sqrt(alphas_vtot_sq./vtot_recalc_sq);
        v_X_step=v_X.*Ekin_corr;
        v_Z_step=v_Z.*Ekin_corr;
        v_phi_step=v_phi.*Ekin_corr;

		alphas_psi_value=interp2_XZ(interp_x,interp_z,psi_XZsmall_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
        alphas_psiH_value=interp2_XZ(interp_x,interp_z,psiH_XZsmall_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
        
        time_stamp=ceil((time_step)/TIME_STAMP_PRECISION);
        alphas_omega=wrap2pi(alphas_theta-alphas_pos_phi_wrapped);
        %alphas_psi_star(INNER_PART)=interp2(1:size_r,scale_theta,psi_star_omega_map',alphas_psi(INNER_PART), alphas_omega(INNER_PART),'*linear');
        if ~isempty(OUTER_PART)
            alphas_psi_value_corr(OUTER_PART)=alphas_psi_value(OUTER_PART);
        end
        alphas_psi_value_corr(INNER_PART)=alphas_psiH_value(INNER_PART)+alphas_psi_star(INNER_PART);
        
        % rescaling the value of pphi0 to our knowledge of the toroidal
        % speed ans psi value
        if CALCULATE_TRUE_PPHI==0
            alphas_pphi0=(mHe/eV)*alphas_Rpos.*v_phi-(ZHe)*alphas_psi_value_corr;
        end
        % amplitude of the B field and potential at half time step local positions
%         update_Gfields_collapse;

        alphas_vpll=v_X_step.*bX+v_Z_step.*bZ+v_phi_step.*bphi;
        vExB_X=(EZ.*bphi-Ephi.*bZ)./alphas_Bfield;
        vExB_Z=(Ephi.*bX-EX.*bphi)./alphas_Bfield;
        vExB_phi=(EX.*bZ-EX.*bX)./alphas_Bfield;
        alphas_vE_sq=vExB_X.^2+vExB_Z.^2+vExB_phi.^2;
        alphas_Epll=0.5*(mHe/eV)*alphas_vpll.^2;
        alphas_Eperp=max(0.5*(mHe/eV)*(v_X_step.^2+v_Z_step.^2+v_phi_step.^2)-alphas_Epll,0);
        alphas_vperp=sqrt(2*alphas_Eperp*eV/mHe);
        alphas_rhoL=(mHe/eV)*(1/ZHe)*alphas_vperp./alphas_Bfield;
        dr_X=interp2(scale_theta,1:Nradial,dr_X_PR_map', alphas_omega,alphas_psi,'*linear');
        dr_Z=interp2(scale_theta,1:Nradial,dr_Z_PR_map', alphas_omega,alphas_psi,'*linear');
        outcast=find(alphas_psi>simulation_size_r+22);
        if (~isempty(outcast))
            recast=randi(Nalphas_simulated,size(outcast,1),1);
            reposition_lost_particles_3DG;
            alphas_ejected(outcast)=1;
            disp(strcat('number of ejected particles =  ',num2str(size(outcast,1))));
        end
		pos_X_gc=(mHe/eV)*(1/ZHe)*(v_Z_step.*bphi-v_phi_step.*bZ)./alphas_Bfield;
        pos_Z_gc=(mHe/eV)*(1/ZHe)*(v_phi_step.*bX-v_X_step.*bphi)./alphas_Bfield;

        Xpos_output(time_stamp,:)=alphas_pos_x+pos_X_gc;
        Zpos_output(time_stamp,:)=alphas_pos_z+pos_Z_gc;
        phipos_output(time_stamp,:)=alphas_pos_phi;
        psipos_output(time_stamp,:)=alphas_psi;
        vparallel_output(time_stamp,:)=alphas_vpll;
        Eperp_output(time_stamp,:)=alphas_Eperp;
        vEsq_output(time_stamp,:)=alphas_vE_sq;
        theta_output(time_stamp,:)=alphas_theta;
        Ekin_output(time_stamp,:)=alphas_Ekin;
        vphi_output(time_stamp,:)=v_phi_step;
        pphi_output(time_stamp,:)=alphas_pphi0;
        psi_star_output(time_stamp,:)=alphas_psi_star;
        psi_value_output(time_stamp,:)=alphas_psi_value_corr;
        rhoL_output(time_stamp,:)=alphas_rhoL;
        Epot_output(time_stamp,:)=alphas_Epot;
        Bfield_output(time_stamp,:)=alphas_Bfield;
        vEradial_output(time_stamp,:)=vExB_X.*dr_X+vExB_Z.*dr_Z;
%         pphi0_output(time_stamp,:)=alphas_pphi0;
%         alphas_vpll(1:4)
        Etot_output(time_stamp,:)=alphas_Etot;
		

    end
%     
    if (DISPLAY_OUTPUTS==1)
        if (mod(time_step+25,50)==0)
            PART_POP=find(~alphas_ejected);
            Xpos_int=alphas_pos_x(PN);
            Zpos_int=alphas_pos_z(PN);
            
%             figure(2);
%             
%             phi_rank=round(alphas_pos_phi_wrapped(PN)/DTHETA)+1;
%             minus_phi_rank=min(NB_THETA-phi_rank+1,NB_THETA);
%             Xpos=alphas_pos_x(PN);
%             Zpos=alphas_pos_z(PN);
%             %             vEXphi(:,:)=vExB_X_map_phi(phi_rank,:,:);
%             %             vEZphi(:,:)=vExB_Z_map_phi(phi_rank,:,:);
%             %             vEmapPR=(vEXphi.*dr_X_PR_map(:,1:size_r)+vEZphi.*dr_Z_PR_map(:,1:size_r));
%             
%             contour(scale_X+R0,scale_Z-Z_axis,(psi_norm_XZsmall_map*5*1e-5)',[257 257]*5*1e-5,'b','LineWidth',4);colorbar;
%             xlim([-0.3 0.35]+R0);
%             ylim([-0.45 0.45]);
%             axis xy
%             pause(0.1)
        end
    end
    
    if (SAVE_DATA_FILE==1) && ((mod(time_step,50)==0)||(time_step==1))
        toc
        if (DISPLAY_OUTPUTS==1)
%%            
            figure(2);
            clf;
            set(gca,'fontsize',20)
%             contour(scale_X+R0,scale_Z-Z_axis,(psi_norm_XZsmall_map*5*1e-5)',[257 257]*5*1e-5,'b','LineWidth',4);colorbar;
            hold on
            plot(alphas_pos_x(RED_PART)+R0,alphas_pos_z(RED_PART),'r.');
            plot(alphas_pos_x(BLUE_PART)+R0,alphas_pos_z(BLUE_PART),'b.');
            
            phi_rank=round(alphas_pos_phi_wrapped(PN)/DTHETA)+1;
            minus_phi_rank=min(NB_THETA-phi_rank+1,NB_THETA);
            Xpos=alphas_pos_x(PN);
            Zpos=alphas_pos_z(PN);
            %             vEXphi(:,:)=vExB_X_map_phi(phi_rank,:,:);
            %             vEZphi(:,:)=vExB_Z_map_phi(phi_rank,:,:);
            %             vEmapPR=(vEXphi.*dr_X_PR_map(:,1:size_r)+vEZphi.*dr_Z_PR_map(:,1:size_r));
            
            clear psi_star_map_phi_rank;
            psi_star_map_phi_rank(:,:)=[psi_star_omega_map(:,minus_phi_rank:NB_THETA)  psi_star_omega_map(:,1:minus_phi_rank-1)];
            psi_star_map_phi_rank=psi_star_map_phi_rank';
            psidata=reshape(psi_star_map_phi_rank(:,1:size_r),NB_THETA*size_r,1);
            
            psistarmap=griddata(finesse_data_X,finesse_data_Z,psidata,XXsmall,ZZsmall,'cubic');
            psistarmap=psistarmap';
            psistarmap(isnan(psistarmap))=0;
            contour(scale_X+R0,scale_Z-Z_axis,psistarmap',(-8.0:0.1:0.5)*1e-3,'linewidth',2);colorbar;
            xlim([-0.25 0.3]+R0);
            ylim([-0.35 0.35]);
            axis xy
            
%             disp('alphas_psi_star(PN)=')
%             disp(alphas_psi_star(PN))
            
            hold on;
            
            plot([ Xpos_int Xpos]+R0,[ Zpos_int Zpos],'k','LineWidth',7);
            plot([ Xpos_int Xpos]+R0,[ Zpos_int Zpos],'g','LineWidth',4);
%             plot([Xpos_prev Xpos_int Xpos],[Zpos_prev Zpos_int Zpos],'k','LineWidth',4);
            Xpos_prev=Xpos;
            Zpos_prev=Zpos;
            %%
            pause(0.4);
            if (time_step<10000)
                frame_name='0';
            else
                frame_name='';
            end
            if (time_step<1000)
                frame_name=strcat(frame_name,'0');
            end
            if (time_step<100)
                frame_name=strcat(frame_name,'0');
            end
            frame_name=strcat(frame_name,num2str(time_step));
            filename='cartoon\t';
            filename=strcat(filename,frame_name,'.bmp');
            disp(filename);
            title(strcat('phi rank #',num2str(phi_rank)));
           
            F=getframe;
            [im,map] = frame2im(F);    %Return associated image data
            
            imwrite(im,filename,'bmp');
        end

        
        tic
        clear outcast recast
        % correct the particles that are out of simulation domain
        % by giving them a position randomly in the initial distribution
        outcast=find(alphas_psi>simulation_size_r+22);
        recast=randi(Nalphas_simulated,size(outcast,1),1);
        if (~isempty(outcast))
            reposition_lost_particles_3DG;
            alphas_ejected(outcast)=1;
        end

        disp('------------------------------------------------')
        disp(strcat('time_step # ',num2str(time_step),' ; time = ',num2str(time)));
        disp(strcat('time_stamp # ',num2str(time_stamp),' ; time = ',num2str(time)));
        disp(strcat('number of repositionned particles =  ',num2str(size(outcast,1))));
        disp(strcat('total number of ejected particles =  ',num2str(size(find(alphas_ejected),1))));

        save (SAVENAME ,'Epot_output','Etot_output','Bfield_output','Xpos_output','Zpos_output','phipos_output','theta_output',...
            'psi_value_output','psi_star_output','vparallel_output','Ekin_output','Eperp_output','rhoL_output','psipos_output',...
            'vphi_output','vEradial_output','pphi_output','time_scale','alphas_ejected');
        disp('------------------------------------------------')
    end



end

%%
toc

Xpos_output=Xpos_output';
Zpos_output=Zpos_output';
phipos_output=phipos_output';
vparallel_output=vparallel_output';
Eperp_output=Eperp_output';
vEsq_output=vEsq_output';
theta_output=theta_output';
psipos_output=psipos_output';
Ekin_output=Ekin_output';
vphi_output=vphi_output';
psi_star_output=psi_star_output';
pphi_output=pphi_output';
% kin_energy_corr_output=kin_energy_corr_output';
psi_value_output=psi_value_output';
rhoL_output=rhoL_output';
Epot_output=Epot_output';
Bfield_output=Bfield_output';
vEradial_output=vEradial_output';
% pphi0_output=pphi0_output';
Etot_output=Etot_output';

Xpos_outputG=Xpos_output;
Zpos_outputG=Zpos_output;
phipos_outputG=phipos_output;
vparallel_outputG=vparallel_output;
Eperp_outputG=Eperp_output;
vEsq_outputG=vEsq_output;
theta_outputG=theta_output;
psipos_outputG=psipos_output;
Ekin_outputG=Ekin_output;
vphi_outputG=vphi_output;
psi_star_outputG=psi_star_output;
pphi_outputG=pphi_output;
psi_value_outputG=psi_value_output;
rhoL_outputG=rhoL_output;
Epot_outputG=Epot_output;
Bfield_outputG=Bfield_output;
vEradial_outputG=vEradial_output;
% pphi0_outputG=pphi0_output;
Etot_outputG=Etot_output;

alphas_ejected_G=alphas_ejected;
time_scale_G=time_scale;

if SAVE_DATA_FILE==1
     save (SAVENAME ,'Epot_outputG','Etot_outputG','Bfield_outputG','Xpos_outputG','Zpos_outputG','phipos_outputG','theta_outputG',...
         'psi_value_outputG','psi_star_outputG','vparallel_outputG','Ekin_outputG','Eperp_outputG','rhoL_outputG','psipos_outputG',...
         'vphi_outputG','vEradial_outputG','pphi_outputG','time_scale_G','alphas_ejected_G','alphas_omega0','alphas_omega','alphas_mm0','alphas_mm_part');
end
disp('**************************************************************');
disp('done');
Nalphas_simulated


