time_step_values=zeros(NB_PROCESS,1);
gamma_TAE_values=zeros(NB_PROCESS,1);
gamma_TAE_vD_values=zeros(NB_PROCESS,1);
gamma_TAE_local=0;
gamma_TAE_vD_local=0;
nTAE

TOROIDAL_FIELD_DIRECTION=sign(mean(mean(Bphi_XZsmall_map)))
PSI_CORE_SIGN=sign(psi_scale(1))
% PSI_STAR_SIGN=sign(psi_star_initial(round(0.5*size_r)))

% psi_star_map_phi_evol=interp1(1:1001,psi_star_2D_evol_interp,1:101,'cubic');

% initialize completely the maps for pre-collapse calculations





if SAVE_DATA_FILE==0
    disp('no data file for this simulation')
end




%run('calculate_collapse_frozen_drift_speed_maps')
% psi_star_map_phi_evol=interp1(1:1001,psi_star_map_phi_evol_interp,1:101,'cubic');
% run('calculate_rotB_vDcurv');
% rotational_b_pll=(rotB_phi.*bphi_XZ_map+rotB_Z.*bZ_XZ_map+rotB_X.*bX_XZ_map)./Btot_XZ_map;


% for correction of Eperp
BMAX=max(max(Btot_XZ_map))
BMIN=min(min(Btot_XZ_map(Btot_XZ_map>1)))


if LOAD_DATA_FILE==0
    REINIT_PERP_SPEED=1;
    
    Nalphas_psi_pos=10;
    Nalphas_lambda=Nalphas_psi_pos*8;
    Nalphas_simulated=Nalphas_lambda*8;
    
    alphas_pos_z=2.4*[  -0.47  -0.45  -0.43 -0.4 -0.33   0.33  0.4 0.43 0.45 0.47 ]
    alphas_pos_z=[alphas_pos_z alphas_pos_z alphas_pos_z alphas_pos_z];
    % alphas_pos_z=[alphas_pos_z alphas_pos_z  ];
    alphas_pos_x=alphas_pos_z*0-0.01;
    
    alphas_pos_x=[alphas_pos_x alphas_pos_x];
    alphas_pos_z=[alphas_pos_z alphas_pos_z];
    
    alphas_pos_x=[alphas_pos_x alphas_pos_x alphas_pos_x alphas_pos_x alphas_pos_x alphas_pos_x alphas_pos_x alphas_pos_x];
    alphas_pos_z=[alphas_pos_z alphas_pos_z alphas_pos_z alphas_pos_z alphas_pos_z alphas_pos_z alphas_pos_z alphas_pos_z];
    
    alphas_pos_z=alphas_pos_z' + Z_axis;
    alphas_pos_x=alphas_pos_x' + Raxis-R0;
    
    alphas_pos_phi=ones(1,Nalphas_psi_pos)*0.39;
    alphas_pos_phi=[alphas_pos_phi 3*alphas_pos_phi ];
    alphas_pos_phi=[alphas_pos_phi alphas_pos_phi+0.39*4 alphas_pos_phi+0.39*8 alphas_pos_phi+0.39*12];
    
    alphas_pos_phi=[alphas_pos_phi alphas_pos_phi alphas_pos_phi alphas_pos_phi alphas_pos_phi alphas_pos_phi alphas_pos_phi alphas_pos_phi];
    alphas_pos_phi=alphas_pos_phi';
    
    alphas_Bfield=interp2(scale_X,scale_Z,Btot_XZ_map',alphas_pos_x,alphas_pos_z,'*linear');
    
    alphas_Ekin=EKIN0_FAC*ones(1,Nalphas_psi_pos)*8*1e3;  % keV
    alphas_Ekin=[alphas_Ekin alphas_Ekin];
    alphas_Ekin=[alphas_Ekin alphas_Ekin alphas_Ekin alphas_Ekin];
    alphas_Ekin=alphas_Ekin';
    
    alphas_vpll=(ones(Nalphas_lambda,1)*0.1);
    alphas_vpll=[-9*alphas_vpll.*sqrt(alphas_Ekin*2*eV/mHe) ; -8.5*alphas_vpll.*sqrt(alphas_Ekin*2*eV/mHe) ; -8*alphas_vpll.*sqrt(alphas_Ekin*2*eV/mHe) ; -7.5*alphas_vpll.*sqrt(alphas_Ekin*2*eV/mHe) ; ...
        7.5*alphas_vpll.*sqrt(alphas_Ekin*2*eV/mHe) ; 8*alphas_vpll.*sqrt(alphas_Ekin*2*eV/mHe) ; 8.5*alphas_vpll.*sqrt(alphas_Ekin*2*eV/mHe) ; 9*alphas_vpll.*sqrt(alphas_Ekin*2*eV/mHe)];
    
    alphas_Ekin=[alphas_Ekin ; alphas_Ekin ; alphas_Ekin ; alphas_Ekin ;  alphas_Ekin ; alphas_Ekin ; alphas_Ekin ; alphas_Ekin];
    
    alphas_mm=(alphas_Ekin-0.5*((mHe/eV)*alphas_vpll.^2))./alphas_Bfield;
    
else
    load(INPUTFILE);
    
    alphas_Ekin=real(alphas_Ekin);
    alphas_mm=max(alphas_mm,0);
    
%     PART_POP_TAE=find((abs(alphas_vpll)<9*1e6).*(abs(alphas_vpll)>5*1e6));
%     PART_POP_TAE=find((abs(alphas_pphi0)<20).*(abs(alphas_pphi0)>14).*(alphas_psi<0.9*pTAE_sup));
    % PART_POP_TAE=find((alphas_psi<2.0*pTAE_sup));
	%half of the population
	%PART_POP_TAE=PART_POP_TAE(1:2:end);
    
    PARTICLES_SPLIT=floor(Nalphas_simulated/NB_PROCESS)
	PART_POP_TAE=(1:1:NB_PROCESS*PARTICLES_SPLIT);
	Nalphas_simulated=length(PART_POP_TAE)
%     FILTER_PART=rand(length(PART_POP_TAE),1)*(length(PART_POP_TAE)-1)+1;
%     FILTER_PART=unique(round(FILTER_PART),'stable');
    PART_POP_TAE=PART_POP_TAE((PROCESS_NUMBER-1)*PARTICLES_SPLIT+1:(PROCESS_NUMBER)*PARTICLES_SPLIT);
    alphas_pos_x=alphas_pos_x(PART_POP_TAE);
    alphas_pos_z=alphas_pos_z(PART_POP_TAE);
    alphas_pos_phi=alphas_pos_phi(PART_POP_TAE);
    alphas_Ekin=alphas_Ekin(PART_POP_TAE);
    alphas_mm=alphas_mm(PART_POP_TAE);
    alphas_vpll=alphas_vpll(PART_POP_TAE);
    alphas_psi=alphas_psi(PART_POP_TAE);
	if ~exist('alphas_weight')
		alphas_weight=ones(length(PART_POP_TAE),1);
	else
		alphas_weight=alphas_weight(PART_POP_TAE);
	end
    v_X=v_X(PART_POP_TAE);
    v_Z=v_Z(PART_POP_TAE);
    v_phi=v_phi(PART_POP_TAE);
    %alphas_pphi0=alphas_pphi0(PART_POP_TAE);
 	%r_avg=r_avg(PART_POP_TAE);
    alphas_Bfield=interp2(scale_X,scale_Z,Btot_XZ_map',alphas_pos_x,alphas_pos_z,'*linear');
    
    Nalphas_simulated=length(alphas_vpll)
    
end

% trying to fix double precision issue with the outputs
alphas_pos_x=double(alphas_pos_x);
alphas_pos_z=double(alphas_pos_z);
alphas_pos_phi=double(alphas_pos_phi);
alphas_Ekin=double(alphas_Ekin);
alphas_mm=double(alphas_mm);
alphas_vpll=double(alphas_vpll);
alphas_psi=double(alphas_psi);
alphas_weight=double(alphas_weight);
v_X=double(v_X);
v_Z=double(v_Z);
v_phi=double(v_phi);
alphas_Bfield=double(alphas_Bfield);


alphas_ejected=zeros(Nalphas_simulated,1);
alphas_Eperp=alphas_Bfield.*alphas_mm;

alphas_Epot=zeros(Nalphas_simulated,1);
part2wave_power=zeros(Nalphas_simulated,1);
part2wave_vD_power=zeros(Nalphas_simulated,1);


%Bavg=mean(Btot_XZ_map(round(0.8*mid_X),:));
Bavg=mean(Btor_PR_map(:,1));

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
Xpos_output=double(zeros(NB_TIME_STAMPS,Nalphas_simulated));
Xpos_part_output=double(zeros(NB_TIME_STAMPS,Nalphas_simulated));
Zpos_output=double(zeros(NB_TIME_STAMPS,Nalphas_simulated));
phipos_output=double(zeros(NB_TIME_STAMPS,Nalphas_simulated));
psipos_output=double(zeros(NB_TIME_STAMPS,Nalphas_simulated));
vD_output=double(zeros(NB_TIME_STAMPS,Nalphas_simulated));
vparallel_output=double(zeros(NB_TIME_STAMPS,Nalphas_simulated));
Eperp_output=double(zeros(NB_TIME_STAMPS,Nalphas_simulated));
vEsq_output=double(zeros(NB_TIME_STAMPS,Nalphas_simulated));
theta_output=double(zeros(NB_TIME_STAMPS,Nalphas_simulated));
Ekin_output=double(zeros(NB_TIME_STAMPS,Nalphas_simulated));
Ekin_gc_output=double(zeros(NB_TIME_STAMPS,Nalphas_simulated));
vphi_output=double(zeros(NB_TIME_STAMPS,Nalphas_simulated));
psi_star_output=double(zeros(NB_TIME_STAMPS,Nalphas_simulated));
pphi_output=double(zeros(NB_TIME_STAMPS,Nalphas_simulated));
% kin_energy_corr_output=double(zeros(NB_TIME_STAMPS,Nalphas_simulated));
psi_value_output=double(zeros(NB_TIME_STAMPS,Nalphas_simulated));
rhoL_output=double(zeros(NB_TIME_STAMPS,Nalphas_simulated));
Epot_output=double(zeros(NB_TIME_STAMPS,Nalphas_simulated));
Bfield_output=double(zeros(NB_TIME_STAMPS,Nalphas_simulated));
vEradial_output=double(zeros(NB_TIME_STAMPS,Nalphas_simulated));
% pphi0_output=double(zeros(NB_TIME_STAMPS,Nalphas_simulated));
mm_output=double(zeros(NB_TIME_STAMPS,Nalphas_simulated));
Etot_output=double(zeros(NB_TIME_STAMPS,Nalphas_simulated));
part2wave_power_output=double(zeros(NB_TIME_STAMPS,Nalphas_simulated));
part2wave_vD_power_output=double(zeros(NB_TIME_STAMPS,Nalphas_simulated));

part2wave_power_evol=zeros(NB_GYRO_STAMPS,1);
part2wave_vD_power_evol=zeros(NB_GYRO_STAMPS,1);
WTAE_evol=zeros(NB_GYRO_STAMPS,1);
gamma_TAE_evol=zeros(NB_GYRO_STAMPS,1);
gamma_TAE_vD_evol=zeros(NB_GYRO_STAMPS,1);
Ekin_part_evol=zeros(NB_GYRO_STAMPS,1);
Etot_th_part_evol=zeros(NB_GYRO_STAMPS,1);
Etot_part_evol=zeros(NB_GYRO_STAMPS,1);
Epart_tot_vD_rel_evol=zeros(NB_GYRO_STAMPS,1);
Epart_tot_rel_evol=zeros(NB_GYRO_STAMPS,1);

power_exchange_evol=double(zeros(NB_TIME_STAMPS,Nalphas_simulated));
pphi_evol=double(zeros(NB_TIME_STAMPS,Nalphas_simulated));

gamma_TAE_record=0;
% WTAE_record=0;
Ekin_record=0;
Etot_record=0;
Etot_th_record=0;
Pexchange_part_record=zeros(Nalphas_simulated,1);
pphi_record=zeros(Nalphas_simulated,1);



alphas_psi=zeros(Nalphas_simulated,1);
alphas_theta=zeros(Nalphas_simulated,1);
alphas_psi_gc=zeros(Nalphas_simulated,1);
alphas_theta_gc=zeros(Nalphas_simulated,1);
alphas_r_value=zeros(Nalphas_simulated,1);
alphas_s_value=zeros(Nalphas_simulated,1);
alphas_q_value=zeros(Nalphas_simulated,1);
alphas_Apll_value=zeros(Nalphas_simulated,1);

v_phi_corr_integ=zeros(Nalphas_simulated,1);
vExB_X=zeros(Nalphas_simulated,1);
vExB_Z=zeros(Nalphas_simulated,1);
vExB_phi=zeros(Nalphas_simulated,1);
bX=zeros(Nalphas_simulated,1);
bZ=zeros(Nalphas_simulated,1);
bphi=zeros(Nalphas_simulated,1);
BX_gc=zeros(Nalphas_simulated,1);
BZ_gc=zeros(Nalphas_simulated,1);
Bphi_gc=zeros(Nalphas_simulated,1);
gB_X=zeros(Nalphas_simulated,1);
gB_Z=zeros(Nalphas_simulated,1);
alphas_psi_star=zeros(Nalphas_simulated,1);
alphas_psi_value=zeros(Nalphas_simulated,1);
alphas_psiH_value=zeros(Nalphas_simulated,1);
alphas_psi_value_corr=zeros(Nalphas_simulated,1);
INDEX_LIST_OMEGA=ones(Nalphas_simulated,1);
% alphas_kin_energy_corr=zeros(Nalphas_simulated,1);
EX=zeros(Nalphas_simulated,1);
EZ=zeros(Nalphas_simulated,1);
Ephi=zeros(Nalphas_simulated,1);
alphas_ejected=zeros(Nalphas_simulated,1);

NON_EJECTED_POP=(~alphas_ejected);
disp(strcat('number of ejected particles =  ',num2str(size(alphas_ejected,1))));
EJECTED_POP=find(alphas_ejected);
			


if (CALCULATE_TRUE_PPHI==1)
    alphas_grad_psi_star=zeros(Nalphas_simulated,1);
    alphas_grad_Phi=zeros(Nalphas_simulated,1);
    alphas_vphi_grad_psi_star=zeros(Nalphas_simulated,1);
end

% intermediate calculations arrays
delta_E=zeros(Nalphas_simulated,1);
dE=zeros(Nalphas_simulated,1);
dE_prev=zeros(Nalphas_simulated,1);
delta_E_th=zeros(Nalphas_simulated,1);
delta_pphi_th=zeros(Nalphas_simulated,1);
alphas_vpll_int=alphas_vpll;
vparallel_tilde=alphas_vpll;


% initial field and position
alphas_psi_value=interp2(scale_X,scale_Z,psi_XZsmall_map',alphas_pos_x,alphas_pos_z,'*linear');
if PSI_CORE_SIGN<0
    alphas_psi_value=max(alphas_psi_value,-psi_global);
else
    alphas_psi_value=min(alphas_psi_value,psi_global);
end
alphas_psi=interp1(psi_scale,1:Nradial,alphas_psi_value);

alphas_theta=interp2(scale_X,scale_Z,theta_XZsmall_map',alphas_pos_x,alphas_pos_z,'*nearest');


frame_rank_next_prev=0;
time=0;
time_oscill=0;

% Initialize fields values
frame_rank_precise=1;
frame_rank_precise_prev=frame_rank_precise

frame_rank_prev=floor((frame_rank_precise-1))+1;
frame_rank_next=ceil((frame_rank_precise-1))+1;

frame_rank=frame_rank_prev;

load_Emaps_frame_rank;
load_Bmaps_frame_rank;


BsX_map_phi_next=BsX_map_phi;
BsZ_map_phi_next=BsZ_map_phi;
Bsphi_map_phi_next=Bsphi_map_phi;

% BsX_map_phi_next=bsX_map_phi.*Bstar_map_phi;
% BsZ_map_phi_next=bsZ_map_phi.*Bstar_map_phi;
% bphi_map_phi_next=bphi_map_phi;
Btot_map_phi=sqrt((BsX_map_phi_next+BpolX_map_phi_ini).^2+(BsZ_map_phi_next+BpolZ_map_phi_ini).^2+(Bsphi_map_phi+Btor_map_phi_ini).^2);
B_PR_map_phi_next=Btot_map_phi;
% psi_star_map_next=(kTAE/omega_TAE)*Epot_map_phi;

%     Efield_X_map_phi_next=zeros(NB_PHI,NB_THETA,size_r);
%     Efield_Z_map_phi_next=zeros(NB_PHI,NB_THETA,size_r);
%     Epot_map_phi_next=zeros(NB_PHI,NB_THETA,size_r);


Efield_X_map_phi_next=Efield_X_map_phi;
Efield_Z_map_phi_next=Efield_Z_map_phi;
% Epot_map_phi_next=Epot_map_phi;

% grad_psi_star_map_phi_next=grad_psi_star_map_phi;


% grad_Phi_map_phi_next=FAST_SAWTOOTH*grad_Phi_tor_map_phi;
    
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

% if (DISPLAY_OUTPUTS==1)
%     PN=12+Nalphas_lambda*3
%     phi_rank=round(alphas_pos_phi_wrapped(PN)/DPHI)+1;
%     Xpos=alphas_pos_x(PN)
%     Zpos=alphas_pos_z(PN)
%     Xpos_prev=Xpos;
%     Zpos_prev=Zpos;
% end

delta_t_coef=h/INTER_FRAME_TIME;

f_counter=0;
update_E_b_fields;
% set initial perturbed fields to 0

% bsX_map_phi=0*bsX_map_phi;
% bsZ_map_phi=0*bsZ_map_phi;
% Btot_map_phi=Btot_map_phi_ini;
% Efield_X_map_phi=0*Efield_X_map_phi;
% Efield_Z_map_phi=0*Efield_Z_map_phi;
% psi_star_phi_map=0*psi_map_phi;
% Epot_map_phi=0*Epot_map_phi;

alphas_iEpot=zeros(Nalphas_simulated,1);
alphas_iA=zeros(Nalphas_simulated,1);
calculate_potentials;

% direction of the B field at local positions
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

if REINIT_PERP_SPEED==1
    
    v0_X=alphas_vpll.*bX+alphas_vperp.*NX;
    v0_Z=alphas_vpll.*bZ+alphas_vperp.*NZ;
    v0_phi=alphas_vpll.*bphi+alphas_vperp.*Nphi;
    
    
    v_X=v0_X;
    v_Z=v0_Z;
    v_phi=v0_phi;
    alphas_Ekin0=alphas_Ekin;
else
    v0_X=v_X;
    v0_Z=v_Z;
    v0_phi=v_phi;
    alphas_Ekin0=alphas_Ekin;
  
end

%initialize Bfield properly
alphas_mm_part=alphas_mm;
% update_GT_3D_collapse;

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
v_tot_corr_integ=vtot_recalc_sq*0;
Ekin_corr=0.5*(1+sqrt(alphas_vtot_sq./vtot_recalc_sq));
v_X_prev=v_X_prev.*(Ekin_corr);
v_Z_prev=v_Z_prev.*(Ekin_corr);
v_phi_prev=v_phi_prev.*(Ekin_corr);

v_phi_prev_prev=v_phi_prev;
% v_phi_prev=v_phi;
v_phi_step=(1.5*v_phi-0.5*v_phi_prev);


calculate_potentials;

if ~isempty(OUTER_PART)
    alphas_grad_Phi(OUTER_PART)=zeros(size(OUTER_PART,1),1);
    alphas_Epot(OUTER_PART)=zeros(size(OUTER_PART,1),1);
end
% alphas_grad_Phi(INNER_PART)=lininterp3( grad_Phi_tor_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
% alphas_Epot(INNER_PART)=lininterp3( E_potential_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);



% Canonical angular momentum initial value
alphas_psi_value=interp2_XZ(interp_x,interp_z,psi_XZsmall_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
% alphas_psiH_value=interp2_XZ(interp_x,interp_z,psi_XZsmall_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);

alphas_psi_star=-bphi.*alphas_Apll_value.*alphas_Rpos;


% alphas_psi_star(INNER_PART)=lininterp3(psi_star_phi_map,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
alphas_psi_value_corr(INNER_PART)=alphas_psi_value(INNER_PART)+alphas_psi_star(INNER_PART);
alphas_psi_value_corr(OUTER_PART)=alphas_psi_value(OUTER_PART);

% alphas_pphi0_half=alphas_pphi0;
% 
% alphas_pphi0_half(INNER_PART)=alphas_pphi0_half(INNER_PART)-(0.5*h*ZHe).*(alphas_Rpos(INNER_PART)).*(alphas_grad_Phi(INNER_PART));



% for eneergy evolution
alphas_Epot_gc_prev=alphas_Epot;
alphas_psi_star_prev=alphas_psi_star;
alphas_Bfield_prev=alphas_Bfield;


% radial_pos=(a/257)*interp2(scale_X,scale_Z,radial_XZsmall_map',alphas_pos_x,alphas_pos_z,'*linear');
% Bavg=mean(alphas_Bfield);
% alphas_kappa=sqrt((alphas_Ekin*R0+Bavg*alphas_mm.*(radial_pos-R0))./(2*alphas_mm.*radial_pos*Bavg));

%trapping parameter
alphas_lambda=Bavg*alphas_mm./alphas_Ekin;

alphas_Eperp=alphas_Ekin-alphas_Epll;
alphas_vperp=sqrt(2*alphas_Eperp*eV/mHe);
alphas_rhoL=(mHe/eV)*(1/ZHe)*alphas_vperp./alphas_Bfield;


alphas_Etot=alphas_Ekin0+ZHe*alphas_Epot;

alphas_Ekin_half=alphas_Ekin;
alphas_Ekin_prev=alphas_Ekin;

v_X_step=v_X;
v_Z_step=v_Z;
v_phi_step=v_phi;
v_X_next=v_X;
v_Z_next=v_Z;
v_phi_next=v_phi;
delta_vphi=v_phi*0;
adapt_speed_Ekin_G;

interpolate_theta_psi_fromXZ;
INNER_PART=find((alphas_psi<=(pTAE_sup-1)).*(alphas_psi>=(pTAE_inf+1)));
OUTER_PART=find((alphas_psi>(pTAE_sup-1))+(alphas_psi<(pTAE_inf+1))>0);
[IL3D_1 IL3D_2 IL3D_3 IL3D_4 IL3D_5 IL3D_6 IL3D_7 IL3D_8 slopex slopey slopez] = ...
    build_3Dinterp_indexarrays(scale_phi, scale_theta, scale_psi, DPHI,DTHETA,1,alphas_pos_phi_TAE(INNER_PART), alphas_theta(INNER_PART), alphas_psi(INNER_PART));

% amplitude of the B field and potential gradients at local positions
if ~isempty(OUTER_PART)
    alphas_Bfield(OUTER_PART)=interp2_XZ(interp_x,interp_z,Btot_XZ_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4,OUTER_PART);
    alphas_grad_Phi(OUTER_PART)=zeros(size(OUTER_PART,1),1);
    alphas_grad_psi_star(OUTER_PART)=zeros(size(OUTER_PART,1),1);
    alphas_Epot(OUTER_PART)=zeros(size(OUTER_PART,1),1);
%     alphas_psi_star(OUTER_PART)=zeros(size(OUTER_PART,1),1);
end
alphas_Bfield(INNER_PART)=lininterp3( Btot_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);

% alphas_Bfield=interp2_XZ(interp_x,interp_z,Btot_XZ_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);

%alphas_grad_Phi(INNER_PART)=lininterp3( grad_Phi_tor_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
%alphas_grad_psi_star(INNER_PART)=lininterp3( grad_psi_star_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
% alphas_grad_psi_star(INNER_PART)=alphas_grad_psi_star(INNER_PART)./alphas_Rpos(INNER_PART);
%alphas_Epot(INNER_PART)=lininterp3( E_potential_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);

%alphas_psi_star=-(kTAE/omega_TAE)*alphas_Epot.*alphas_Rpos;


% alphas_Epot_step_prev=alphas_Epot;
alphas_psi_value_corr(INNER_PART)=alphas_psi_value(INNER_PART)+alphas_psi_star(INNER_PART);
alphas_psi_value_corr(OUTER_PART)=alphas_psi_value(OUTER_PART);
alphas_pphi0=(mHe/eV)*alphas_Rpos.*v_phi-(ZHe)*alphas_psi_value_corr;
alphas_Bfield_prev=alphas_Bfield;

alphas_psi_star_half_prev=alphas_psi_star;
alphas_psi_star_part_prev=alphas_psi_star;
alphas_psi_star_prev=alphas_psi_star;
alphas_Ekin_prev=alphas_Ekin;
delta_ps=alphas_psi_star-alphas_psi_star_prev;
delta_ps_prev=delta_ps;
delta_B=alphas_Bfield-alphas_Bfield_prev;
delta_pphi=alphas_pphi0*0;
delta_pphi_prev=delta_pphi;
delta_pphi_raw=delta_pphi;
%delta_ps_prev=delta_ps;

alphas_Rpos_int=alphas_Rpos;
% alphas_vphi_grad_psi_star=v_phi.*alphas_grad_psi_star;
% half_delta_t_coef=100*(0.5*h*FAST_SAWTOOTH/tau_cr);

% benchmarking time between two data recordings
disp('**************************************************************');
tic

time=0;
gamma_TAE=0;

if CALCULATE_VD_POWER==1
    vD0=((2*alphas_Ekin)-alphas_mm_part.*alphas_Bfield)/ZHe;
    vDX=interp2_XZ(interp_x,interp_z,vD_X_XZ_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
    vDZ=interp2_XZ(interp_x,interp_z,vD_Z_XZ_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
    vDphi=interp2_XZ(interp_x,interp_z,vD_phi_XZ_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
    part2wave_vD_power=-ZHe*vD0.*(vDX.*EX+vDZ.*EZ+vDphi.*Ephi); 
end

W_TAE=(TAE_amplitude)^2;
P_VD_PART_TAE=0;
P_PART_TAE=0;

alphas_pphi0_prev=alphas_pphi0;
delta_Ekin_vD=zeros(Nalphas_simulated,1);
alphas_Epot_gc=alphas_Epot;
alphas_Epot_prev=alphas_Epot;
alphas_Epot_gc_part_prev=alphas_Epot_gc;
delta_Epot=alphas_Epot_gc-alphas_Epot_gc_part_prev;
alphas_Etot_th=alphas_Etot;
alphas_Etot_prev=alphas_Etot;
alphas_Etot_th_prev=alphas_Etot_th;
alphas_Ekin_th=alphas_Ekin;
alphas_mm_gc=alphas_mm;
Btot_map_phi_prev=Btot_map_phi_ini;
alphas_Bfield_gc=alphas_Bfield;
alphas_Bfield_gc_prev=alphas_Bfield_gc;
Btot_map_phi_prev=Btot_map_phi;
%E_potential_map_phi_prev=E_potential_map_phi;
%calculate_B_gc_value_prev;
v_phi_gc=v_phi;
bphi_gc=bphi;
alphas_Ekin_part=alphas_Ekin;
% Ekin_vD=alphas_Ekin;
Epart_tot_vD_rel_ini=0;
% Epart_tot_vD_rel_ini_prev=Epart_tot_vD_rel_ini;
Epart_tot_rel_ini=0;
% Ekin_tot_part_prev=0;
Ekin_tot_part=NB_PART_RESCALE*eV*sum(alphas_Ekin);
% flag_update_vphi_output=0;
time_stamp=1;
time_stamp_gyro=1;
record_step=0;
% part2wave_vD_power_prev=part2wave_vD_power;
update_G_3D_TAE;

quasi_lin_flag=0;

time_step=1;
save_simulation_data;