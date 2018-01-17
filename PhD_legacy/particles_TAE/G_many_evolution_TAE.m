
close all;
pause(0.2);

REINIT_ALL_MAPS=0;
REINIT_SIMULATION_DATA=1;

if REINIT_ALL_MAPS==1
    clear all
    initialize_folder_names_test;
    filename=strcat(DATA_FOLDER,'physics_constants.mat');
    load(filename);
    filename=strcat(DATA_FOLDER,'tokamak_map_dimensions.mat');
    load(filename);
    filename=strcat(DATA_FOLDER,'TAE_data.mat');
    load(filename);

    DPHI=TAE_angle/(NB_PHI-1);
    size_r=pTAE_sup-pTAE_inf+1;

    mid_X_large=mid_X;
    mid_Z_large=mid_Z;    
    calculate_TAE_drift_speed_maps;
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
    filename=strcat(DATA_FOLDER,'XZsmall_fields_tokamak_TAE.mat');
    load(filename);
    filename=strcat(DATA_FOLDER,'motions_map_dimensions.mat');
    load(filename);
    filename=strcat(DATA_FOLDER,'psi_profiles_kadomtsev.mat');
    load(filename,'psi_star_initial');

    filename=strcat(DATA_FOLDER,'tokamak_PR_map.mat');
    load(filename, 'Rpos_PR_map');
    load(filename, 'radial_r_value_flux');
    load(filename,'BX_PR_map')
    load(filename, 'BZ_PR_map')
    filename=strcat(DATA_FOLDER,'B_fields.mat');
    load(filename,'Bpol_PR_map');
    load(filename,'Btor_PR_map');

    filename=strcat(DATA_FOLDER,'TAE_data.mat');
    load(filename);
    
    Btot_PR_map_ini=sqrt(Bpol_PR_map.^2+Btor_PR_map.^2);
    TAE_angle=2*pi/nTAE
    NB_PHI
    DPHI=TAE_angle/(NB_PHI-1);
    size_r=size_r_TAE;
    scale_phi=TAE_angle*(0:NB_PHI-1)/(NB_PHI-1);
    scale_psi=pTAE_inf:pTAE_sup;
    TAE_amplitude=10.0;
    
    NB_PART_RESCALE=5e15

    REINIT_SIMULATION_DATA=1;

end


initialize_TAE_sim_maps;


TOROIDAL_FIELD_DIRECTION=sign(mean(mean(Bphi_XZsmall_map)))
PSI_CORE_SIGN=sign(psi_scale(1))
PSI_STAR_SIGN=sign(psi_star_initial(round(0.5*size_r)))

% psi_star_map_phi_evol=interp1(1:1001,psi_star_2D_evol_interp,1:101,'cubic');

% initialize completely the maps for pre-collapse calculations
NB_OSCILLATIONS=2;

TPRECISE=4;
TIME_STAMP_PRECISION=10*TPRECISE;

TAE_PERIOD=(2*pi)/omega_TAE
SIMULATION_TIME=TAE_PERIOD*NB_OSCILLATIONS
TIME_STEP_REF_SIZE=TAE_PERIOD/(NB_FRAME)/20/10

DELTA_TIME=TIME_STEP_REF_SIZE/TPRECISE;
h=DELTA_TIME;
NB_TIME_STEPS=round(0.05*SIMULATION_TIME/DELTA_TIME)
INTER_FRAME_TIME=TAE_PERIOD/NB_FRAME

% one time stamp in one loop
NB_TIME_STAMPS=round(NB_TIME_STEPS/TIME_STAMP_PRECISION);
TIME_STAMP_SIZE=round(0.5*20*NB_TIME_STEPS/NB_TIME_STAMPS)
% ten time steps in one loop
time_scale=(1:NB_TIME_STAMPS)*DELTA_TIME*TIME_STAMP_PRECISION*20;

TAE_ANGLE=(2*pi)/nTAE

CALCULATE_VD_POWER=1;



%simulation options
REINIT_PERP_SPEED=0;
SAVE_DATA_FILE=1;
CALCULATE_TRUE_PPHI=0;



% INPUTFILE='initialG_flatD_90keV_pre_collapse.mat'
% SAVENAME='flatD_90keV_collapse_Glisa_fc2h2_G260214.mat'

% INPUTFILE='initialG_D_MB_pre_collapse_b2.mat'
% SAVENAME='D_MB_b2_collapse_Glisa_fc2h2_G260214.mat'

INPUTFILE='initialG_alphas_TAE_counter.mat'
SAVENAME='alphas_TAE_vperp0_lisa_o2h4_G220714.mat'

mHe
ZHe

load(INPUTFILE);

alphas_Ekin=real(alphas_Ekin);
alphas_mm=max(alphas_mm,0);

% PART_POP_TAE=find((abs(alphas_vpll)<8.0*1e6).*(abs(alphas_vpll)>7*1e6));
% alphas_pos_x=alphas_pos_x(PART_POP_TAE);
% alphas_pos_z=alphas_pos_z(PART_POP_TAE);
% alphas_pos_phi=alphas_pos_phi(PART_POP_TAE);
% alphas_Ekin=alphas_Ekin(PART_POP_TAE);
% alphas_mm=alphas_mm(PART_POP_TAE);
% alphas_vpll=alphas_vpll(PART_POP_TAE);
% alphas_psi=alphas_psi(PART_POP_TAE);
% v_X=v_X(PART_POP_TAE);
% v_Z=v_Z(PART_POP_TAE);
% v_phi=v_phi(PART_POP_TAE);
% alphas_pphi0=alphas_pphi0(PART_POP_TAE);



Nalphas_simulated=length(alphas_vpll)





if SAVE_DATA_FILE==0
    disp('no data file for this simulation')
end

% if DISPLAY_OUTPUTS==0
    disp('no graphics will be displayed during this simulation')
% else
%    filename=strcat(DATA_FOLDER,'tokamak_PR_map.mat');
%    load(filename);
%    [XXsmall ZZsmall]=meshgrid(scale_X,scale_Z);
%    finesse_data_X=reshape((Rpos_PR_map(:,1:size_r)-R0),NP*size_r,1);
%    finesse_data_Z=reshape(Z_PR_map(:,1:size_r),NP*size_r,1);
% end


%run('calculate_collapse_frozen_drift_speed_maps')
% psi_star_map_phi_evol=psi_star_2D_evol;

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




alphas_Bfield=interp2(scale_X,scale_Z,Btot_XZ_map',alphas_pos_x,alphas_pos_z,'*linear');
alphas_ejected=zeros(Nalphas_simulated,1);
alphas_eject_posX=zeros(Nalphas_simulated,1);
alphas_eject_posZ=zeros(Nalphas_simulated,1);
alphas_eject_vpll=zeros(Nalphas_simulated,1);
alphas_Eperp=alphas_Bfield.*alphas_mm;

alphas_Epot=zeros(Nalphas_simulated,1);
part2wave_power=zeros(Nalphas_simulated,1);
part2wave_vD_power=zeros(Nalphas_simulated,1);

%alphas_Ekin=alphas_Ekin+0.5*(mHe/eV)*alphas_vdrift_sq;

Bavg=mean(alphas_Bfield);

%initial accurate energy values
alphas_vE_sq=zeros(Nalphas_simulated,1);

alphas_Epll=0.5*(mHe/eV)*alphas_vpll.^2;
alphas_Ekin=alphas_Eperp+alphas_Epll;

%initial position values
alphas_Ekin0=alphas_Ekin;
alphas_mm0=alphas_mm;
X0=alphas_pos_x;
Z0=alphas_pos_z;
phi0=alphas_pos_phi;


alphas_Eperp=max(alphas_Ekin-alphas_Epll,0);
alphas_vperp=sqrt(2*alphas_Eperp*eV/mHe);
alphas_rhoL=(mHe/eV)*(1/ZHe)*alphas_vperp./alphas_Bfield;


disp('**************************************************************');

%initializing simulation arrays
alphas_grad_psi_star=zeros(Nalphas_simulated,1);
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
EX=zeros(Nalphas_simulated,1);
EZ=zeros(Nalphas_simulated,1);
Ephi=zeros(Nalphas_simulated,1);

part2wave_power_evol=zeros(NB_TIME_STAMPS,1);
part2wave_vD_power_evol=zeros(NB_TIME_STAMPS,1);
WTAE_evol=zeros(NB_TIME_STAMPS,1);
gamma_TAE_evol=zeros(NB_TIME_STAMPS,1);
Ekin_part_evol=zeros(NB_TIME_STAMPS,1);
Etot_part_evol=zeros(NB_TIME_STAMPS,1);
power_exchange_evol=zeros(NB_TIME_STAMPS,Nalphas_simulated);
pphi_evol=zeros(NB_TIME_STAMPS,Nalphas_simulated);

gamma_TAE_record=zeros(TIME_STAMP_SIZE,1);
WTAE_record=zeros(TIME_STAMP_SIZE,1);
Ekin_record=zeros(TIME_STAMP_SIZE,1);
Etot_record=zeros(TIME_STAMP_SIZE,1);
Pexchange_part_record=zeros(Nalphas_simulated,1);
pphi_record=zeros(Nalphas_simulated,1);

% intermediate calculations arrays
delta_E=zeros(Nalphas_simulated,1);



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

% Initialize drift values
if ~exist('time')
    time=0;
    frame_rank_precise=1;
end
% time_oscill=mod(time,TAE_PERIOD);
% if time<SIMULATION_TIME
%     frame_rank_precise=(NB_FRAME+1)*(time_oscill/TAE_PERIOD);
% else
%     frame_rank_precise=1
% end
frame_rank_prev=floor(frame_rank_precise);
frame_rank_next=ceil(frame_rank_precise);
frame_rank=frame_rank_prev;

load_Emaps_frame_rank;
load_Bmaps_frame_rank;

BsX_map_phi_next=BsX_map_phi;
BsZ_map_phi_next=BsZ_map_phi;
Btot_map_phi=sqrt((BsX_map_phi_next+BpolX_map_phi_ini).^2+(BsZ_map_phi_next+BpolZ_map_phi_ini).^2+Btor_map_phi_ini.^2);
B_PR_map_phi_next=Btot_map_phi;

% psi_star_map_next=psi_map_phi;
Efield_X_map_phi_next=Efield_X_map_phi;
Efield_Z_map_phi_next=Efield_Z_map_phi;
Epot_map_phi_next=Epot_map_phi;

frame_rank_next=frame_rank_prev;

    
frame_rank_next=frame_rank_prev;

% estimating psi and theta values
alphas_pos_phi_wrapped=wrap2pi(alphas_pos_phi);
alphas_pos_phi_TAE=wrapTAEangle(alphas_pos_phi_wrapped,TAE_angle);
interpolate_theta_psi_fromXZ;
% sorting out which particles are inside mixing region and which are out
INNER_PART=find((alphas_psi<=(pTAE_sup-1)).*(alphas_psi>=(pTAE_inf+1)));
OUTER_PART=find((alphas_psi>(pTAE_sup-1))+(alphas_psi<(pTAE_inf+1))>0);
[IL3D_1 IL3D_2 IL3D_3 IL3D_4 IL3D_5 IL3D_6 IL3D_7 IL3D_8 slopex slopey slopez] = ...
        build_3Dinterp_indexarrays(scale_phi, scale_theta, scale_psi, DPHI,DTHETA,1,alphas_pos_phi_TAE(INNER_PART), alphas_theta(INNER_PART), alphas_psi(INNER_PART));

delta_t_coef=h/INTER_FRAME_TIME;
f_counter=0;
update_E_b_fields;

% amplitude of the B field and potential at local positions
update_Gfields_TAE;

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

if (REINIT_PERP_SPEED==1)
	v0_X=alphas_vpll.*bX+alphas_vperp.*NX;
	v0_Z=alphas_vpll.*bZ+alphas_vperp.*NZ;
	v0_phi=alphas_vpll.*bphi+alphas_vperp.*Nphi;

	v_X=v0_X;
	v_Z=v0_Z;
	v_phi=v0_phi;
else
	v0_X=v_X;
	v0_Z=v_Z;
	v0_phi=v_phi;
end

alphas_Ekin0=alphas_Ekin;
alphas_mm_part=alphas_mm;

interpolate_theta_psi_fromXZ;
INNER_PART=find((alphas_psi<=(pTAE_sup-1)).*(alphas_psi>=(pTAE_inf+1)));
OUTER_PART=find((alphas_psi>(pTAE_sup-1))+(alphas_psi<(pTAE_inf+1))>0);


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
alphas_Epot(INNER_PART)=lininterp3( E_potential_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);

% Canonical angular momentum initial value
alphas_psi_value=interp2_XZ(interp_x,interp_z,psi_XZsmall_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);

alphas_psi_star=-(kTAE/omega_TAE)*alphas_Epot.*alphas_Rpos;
% alphas_psi_star(INNER_PART)=lininterp3(psi_star_phi_map,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
alphas_psi_value_corr(INNER_PART)=alphas_psi_value(INNER_PART)+alphas_psi_star(INNER_PART);
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
alphas_lambda=Bavg*alphas_mm0./alphas_Ekin;

alphas_Eperp=alphas_Ekin-alphas_Epll;
alphas_vperp=sqrt(2*alphas_Eperp*eV/mHe);
alphas_rhoL=(mHe/eV)*(1/ZHe)*alphas_vperp./alphas_Bfield;


alphas_Etot=alphas_Ekin0+ZHe*alphas_Epot;

alphas_Ekin_half=alphas_Ekin;
alphas_Ekin_prev=alphas_Ekin;

adapt_speed_Ekin_G;

interpolate_theta_psi_fromXZ;
INNER_PART=find((alphas_psi<=(pTAE_sup-1)).*(alphas_psi>=(pTAE_inf+1)));
OUTER_PART=find((alphas_psi>(pTAE_sup-1))+(alphas_psi<(pTAE_inf+1))>0);
[IL3D_1 IL3D_2 IL3D_3 IL3D_4 IL3D_5 IL3D_6 IL3D_7 IL3D_8 slopex slopey slopez] = ...
    build_3Dinterp_indexarrays(scale_phi, scale_theta, scale_psi, DPHI,DTHETA,1,alphas_pos_phi_TAE(INNER_PART), alphas_theta(INNER_PART), alphas_psi(INNER_PART));

% amplitude of the B field and potential gradients at local positions
if ~isempty(OUTER_PART)
    alphas_Bfield(OUTER_PART)=interp2_XZ(interp_x,interp_z,Btot_XZ_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4,OUTER_PART);
    alphas_Epot(OUTER_PART)=zeros(size(OUTER_PART,1),1);
    alphas_psi_star(OUTER_PART)=zeros(size(OUTER_PART,1),1);
end
% alphas_Bfield=interp2_XZ(interp_x,interp_z,Btot_XZ_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
alphas_Bfield(INNER_PART)=lininterp3( Btot_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
alphas_Epot(INNER_PART)=lininterp3( E_potential_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
% alphas_psi_star(INNER_PART)=lininterp3(psi_star_phi_map,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
alphas_psi_star=-(kTAE/omega_TAE)*alphas_Epot.*alphas_Rpos;

% alphas_Epot_step_prev=alphas_Epot;
Btot_map_phi_prev=Btot_map_phi_ini;
alphas_Bfield_gc=alphas_Bfield;
alphas_Bfield_gc_prev=alphas_Bfield;
delta_B=0*alphas_Bfield_gc;

% alphas_psi_star_half_prev=alphas_psi_star;
% alphas_psi_star_part_prev=alphas_psi_star;
% alphas_psi_star_prev=alphas_psi_star;
% alphas_Ekin_prev=alphas_Ekin;
% delta_ps=alphas_psi_star-alphas_psi_star_prev;
% delta_ps_prev=delta_ps;

alphas_Rpos_int=alphas_Rpos;


% benchmarking time between two data recordings
disp('**************************************************************');
tic

gamma_TAE=0;

if CALCULATE_VD_POWER==1
    vD0=((2*alphas_Ekin)-alphas_mm_part.*alphas_Bfield)/ZHe./alphas_Bfield;
    vDX=interp2_XZ(interp_x,interp_z,vD_X_XZ_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
    vDZ=interp2_XZ(interp_x,interp_z,vD_Z_XZ_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
    vDphi=interp2_XZ(interp_x,interp_z,vD_phi_XZ_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
    part2wave_vD_power=-ZHe*vD0.*(vDX.*EX+vDZ.*EZ+vDphi.*Ephi); 
end

% INNER_PART=find(alphas_psi<=(size_r-1));
% OUTER_PART=find(alphas_psi>(size_r-1));
% [IL3D_1 IL3D_2 IL3D_3 IL3D_4 IL3D_5 IL3D_6 IL3D_7 IL3D_8 slopex slopey slopez] = ...
%     build_3Dinterp_indexarrays(scale_phi, scale_theta, scale_psi, DPHI,DTHETA,1,alphas_pos_phi_wrapped(INNER_PART), alphas_theta(INNER_PART), alphas_psi(INNER_PART));

W_TAE=(TAE_amplitude)^2;
P_VD_PART_TAE=0;

alphas_mm_gc=alphas_mm;
Btot_map_phi_prev=Btot_map_phi_ini;
alphas_Bfield_gc=alphas_Bfield;
alphas_Bfield_gc_prev=alphas_Bfield_gc;
alphas_Ekin_part=alphas_Ekin;
Epart_tot_vD_rel_ini=0;
Epart_tot_rel_ini=0;
Ekin_tot_part=NB_PART_RESCALE*eV*sum(alphas_Ekin);
time_stamp=1;
record_step=0;

for time_step=1:NB_TIME_STEPS
    
    time_step_integration_G_TAE;
    time_step_integration_G_TAE;
    time_step_integration_G_TAE;
    time_step_integration_G_TAE;
    time_step_integration_G_TAE;
    time_step_integration_G_TAE;
    time_step_integration_G_TAE;
    time_step_integration_G_TAE;
    time_step_integration_G_TAE;
    time_step_integration_G_TAE;

	time_step_integration_G_TAE;
    time_step_integration_G_TAE;
    time_step_integration_G_TAE;
    time_step_integration_G_TAE;
    time_step_integration_G_TAE;
    time_step_integration_G_TAE;
    time_step_integration_G_TAE;
    time_step_integration_G_TAE;
    time_step_integration_G_TAE;
    time_step_integration_G_TAE;


	outcast=find(alphas_psi>Nradial-2);
	if (~isempty(outcast))
		alphas_eject_posX(outcast)=alphas_pos_x(outcast);
		alphas_eject_posZ(outcast)=alphas_pos_z(outcast);
		alphas_eject_vpll(outcast)=alphas_vpll(outcast);
		alphas_pos_x(outcast)=0;
		alphas_pos_z(outcast)=0;
        alphas_ejected(outcast)=1;
        disp(strcat('number of ejected particles (reset to 0,0) =  ',num2str(size(outcast,1))));
	end
    
    if (mod(time_step,TIME_STAMP_PRECISION)==0)
        W_TAE
        gamma_TAE
        record_step=0;

        adapt_speed_Ekin_G;
%         v_X_step=(1.5*v_X-0.5*v_X_prev);
%         v_Z_step=(1.5*v_Z-0.5*v_Z_prev);
%         v_phi_step=(1.5*v_phi-0.5*v_phi_prev);
        v_X_prev=v_X;
        v_Z_prev=v_Z;
        v_phi_prev=v_phi;
		
		alphas_psi_value=interp2_XZ(interp_x,interp_z,psi_XZsmall_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
       
        time_stamp=ceil((time_step)/TIME_STAMP_PRECISION);

        if ~isempty(OUTER_PART)
            alphas_psi_value_corr(OUTER_PART)=alphas_psi_value(OUTER_PART);
        end
        alphas_psi_value_corr(INNER_PART)=alphas_psi_value(INNER_PART)+alphas_psi_star(INNER_PART);
        
        % rescaling the value of pphi0 to our knowledge of the toroidal
        % speed ans psi value
        if CALCULATE_TRUE_PPHI==0
			alphas_pphi0=(mHe/eV)*alphas_Rpos.*v_phi-(ZHe)*alphas_psi_value_corr;
        end
		
        % amplitude of the B field and potential at half time step local positions
        update_Gfields_TAE;

        alphas_vpll=v_X_step.*bX+v_Z_step.*bZ+v_phi_step.*bphi;
        vExB_X=(EZ.*bphi-Ephi.*bZ)./alphas_Bfield;
        vExB_Z=(Ephi.*bX-EX.*bphi)./alphas_Bfield;
        vExB_phi=(EX.*bZ-EX.*bX)./alphas_Bfield;
        alphas_vE_sq=vExB_X.^2+vExB_Z.^2+vExB_phi.^2;
        alphas_Epll=0.5*(mHe/eV)*alphas_vpll.^2;
        alphas_Eperp=max(0.5*(mHe/eV)*(v_X_step.^2+v_Z_step.^2+v_phi_step.^2)-alphas_Epll,0);
		alphas_mm=alphas_Eperp./alphas_Bfield;
        outcast=find(alphas_psi>Nradial-2);
        alphas_vperp=sqrt(2*alphas_Eperp*eV/mHe);
        alphas_rhoL=(mHe/eV)*(1/ZHe)*alphas_vperp./alphas_Bfield;
        if (~isempty(outcast))
            recast=randi(Nalphas_simulated,size(outcast,1),1);
            reposition_lost_particles_3DG;
            alphas_ejected(outcast)=1;
            disp(strcat('number of ejected particles =  ',num2str(size(outcast,1))));
        end
        
        part2wave_power_evol(time_stamp)=P_PART_TAE;
        part2wave_vD_power_evol(time_stamp)=P_VD_PART_TAE;
        WTAE_evol(time_stamp)=W_TAE;
        Ekin_part_evol(time_stamp)=NB_PART_RESCALE*eV*mean(Ekin_record);
        Etot_part_evol(time_stamp)=NB_PART_RESCALE*eV*mean(Etot_record);;
        gamma_TAE_evol(time_stamp)=mean(gamma_TAE_record);
        Epart_tot_vD_rel_evol(time_stamp)=Epart_tot_vD_rel_ini;
        Epart_tot_rel_evol(time_stamp)=Epart_tot_rel_ini;
        power_exchange_evol(time_stamp,:)=Pexchange_part_record/TIME_STAMP_SIZE;
        pphi_evol(time_stamp,:)=pphi_record/TIME_STAMP_SIZE;
        Pexchange_part_record=zeros(Nalphas_simulated,1);
        pphi_record=zeros(Nalphas_simulated,1);

    end
    
    %     end
    
    if (SAVE_DATA_FILE==1) && (mod(time_step,200)==0)
        toc
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
        disp(strcat('number of repositionned particles =  ',num2str(size(outcast,1))));
        disp(strcat('total number of ejected particles =  ',num2str(size(find(alphas_ejected),1))));
        
        pos_X_gc=(mHe/eV)*(1/ZHe)*(v_Z_step.*bphi-v_phi_step.*bZ)./alphas_Bfield;
        pos_Z_gc=(mHe/eV)*(1/ZHe)*(v_phi_step.*bX-v_X_step.*bphi)./alphas_Bfield;
		pos_X_gc=alphas_pos_x+pos_X_gc;
		pos_Z_gc=alphas_pos_z+pos_Z_gc;

        % filename=strcat('flatD_90keV_collapse_Glisa_fc2h2_G260214record_',num2str(time_step),'.mat')
    save (SAVENAME,'mHe','ZHe','time','frame_rank_precise','pos_X_gc','pos_Z_gc','alphas_pos_x','alphas_pos_z',...
        'alphas_eject_posX','alphas_eject_posZ','alphas_eject_vpll','alphas_pos_phi','alphas_omega','v_X','v_Z','v_phi','alphas_vpll','alphas_mm',...
        'alphas_vE_sq','alphas_psi','alphas_pphi0','alphas_psi_value_corr','alphas_Etot','alphas_Ekin','alphas_grad_Phi','alphas_ejected',...
        'Epart_tot_rel_evol','Epart_tot_vD_rel_evol','WTAE_evol','Ekin_part_evol','Etot_part_evol','gamma_TAE_evol','power_exchange_evol');
        disp('------------------------------------------------')
    end
    
    
    
end

%%
pos_X_gc=(mHe/eV)*(1/ZHe)*(v_Z_step.*bphi-v_phi_step.*bZ)./alphas_Bfield;
pos_Z_gc=(mHe/eV)*(1/ZHe)*(v_phi_step.*bX-v_X_step.*bphi)./alphas_Bfield;
pos_X_gc=alphas_pos_x+pos_X_gc;
pos_Z_gc=alphas_pos_z+pos_Z_gc;

toc


if SAVE_DATA_FILE==1
    save (SAVENAME,'mHe','ZHe','time','frame_rank_precise','pos_X_gc','pos_Z_gc','alphas_pos_x','alphas_pos_z',...
        'alphas_eject_posX','alphas_eject_posZ','alphas_eject_vpll','alphas_pos_phi','alphas_omega','v_X','v_Z','v_phi','alphas_vpll','alphas_mm',...
        'alphas_vE_sq','alphas_psi','alphas_pphi0','alphas_psi_value_corr','alphas_Etot','alphas_Ekin','alphas_grad_Phi','alphas_ejected',...
        'Epart_tot_rel_evol','Epart_tot_vD_rel_evol','WTAE_evol','Ekin_part_evol','Etot_part_evol','gamma_TAE_evol','power_exchange_evol',...
        'pphi_evol');
end
disp('**************************************************************');
disp('done');
Nalphas_simulated
