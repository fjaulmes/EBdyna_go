    load(INPUTFILE);
 
if ~exist('alphas_momentum')
    alphas_momentum=alphas_Ekin*0+MEAN_Q1_ROTATION;
else
    %alphas_momentum=max(alphas_momentum(PART_POP_COLLAPSE),MEAN_Q1_ROTATION);
    %alphas_momentum=alphas_Ekin*0+mean(alphas_momentum);
end
 
    alphas_Ekin=real(alphas_Ekin);
    alphas_mm=max(alphas_mm,0);
	
    PARTICLES_SPLIT=floor(Nalphas_simulated/NB_PROCESS)
    
    PART_POP_COLLAPSE=((PROCESS_NUMBER-1)*PARTICLES_SPLIT+1:1:(PROCESS_NUMBER)*PARTICLES_SPLIT);
    alphas_pos_x=alphas_pos_x(PART_POP_COLLAPSE);
    alphas_pos_z=alphas_pos_z(PART_POP_COLLAPSE);
    alphas_pos_phi=alphas_pos_phi(PART_POP_COLLAPSE);
    alphas_Ekin=alphas_Ekin(PART_POP_COLLAPSE);
    alphas_mm=alphas_mm(PART_POP_COLLAPSE);
    alphas_vpll=alphas_vpll(PART_POP_COLLAPSE);
    alphas_psi=alphas_psi(PART_POP_COLLAPSE);
    alphas_momentum=alphas_momentum(PART_POP_COLLAPSE);
	if ~exist('alphas_weight')
		alphas_weight=ones(length(PART_POP_COLLAPSE),1);
	else
		alphas_weight=alphas_weight(PART_POP_COLLAPSE);
	end
    v_X=v_X(PART_POP_COLLAPSE);
    v_Z=v_Z(PART_POP_COLLAPSE);
    v_phi=v_phi(PART_POP_COLLAPSE);
    %alphas_pphi0=alphas_pphi0(PART_POP_COLLAPSE);
 	%r_avg=r_avg(PART_POP_COLLAPSE);
    alphas_Bfield=interp2(scale_X,scale_Z,Btot_XZ_map',alphas_pos_x,alphas_pos_z,'*linear');
    
    Nalphas_simulated=length(alphas_vpll)
	


%load(INPUTFILE);
mHe
ZHe








alphas_Bfield=interp2(scale_X,scale_Z,Btot_XZ_map',alphas_pos_x,alphas_pos_z,'*linear');
alphas_ejected=zeros(Nalphas_simulated,1);
alphas_eject_posX=zeros(Nalphas_simulated,1);
alphas_eject_posZ=zeros(Nalphas_simulated,1);
alphas_eject_vpll=zeros(Nalphas_simulated,1);
alphas_Eperp=alphas_Bfield.*alphas_mm;

alphas_Epot=zeros(Nalphas_simulated,1);
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



% intermediate calculations arrays
delta_E=zeros(Nalphas_simulated,1);



% initial field and position
alphas_psi_value=interp2(scale_X,scale_Z,psi_XZsmall_map',alphas_pos_x,alphas_pos_z,'*linear');
if PSI_CORE_SIGN<0
    alphas_psi_value=max(alphas_psi_value,-psi_global);
else
    alphas_psi_value=min(alphas_psi_value,psi_global);
end
error('Fix hardcoded radial points!')
alphas_psi=interp1(psi_scale,1:513,alphas_psi_value);

alphas_theta=interp2(scale_X,scale_Z,theta_XZsmall_map',alphas_pos_x,alphas_pos_z,'*nearest');


frame_rank_next_prev=0;

% Initialize drift values
if ~exist('time')
    time=0;
    frame_rank_precise=1;
end
frame_rank_prev=10*floor((frame_rank_precise-1)/10)+1;
frame_rank_next=10*ceil((frame_rank_precise-1)/10)+1;
frame_rank_precise_prev=frame_rank_precise;

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


% grad_Phi_map_phi_next=FAST_SAWTOOTH*grad_Phi_tor_map_phi;
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

% if (DISPLAY_OUTPUTS==1)
%    PN=3;
%    phi_rank=round(alphas_pos_phi_wrapped(PN)/DPHI)+1;
%    Xpos=alphas_pos_x(PN);
%    Zpos=alphas_pos_z(PN);
%    Xpos_prev=Xpos;
%    Zpos_prev=Zpos;
% end

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
% alphas_Epot(INNER_PART)=lininterp3( E_potential_PR_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
alphas_Epot(INNER_PART)=interp2_omega_map(scale_psi(1:size_r),scale_theta,E_potential_omega_map,1,DTHETA,alphas_psi, alphas_omega,INNER_PART);

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
alphas_lambda=Bavg*alphas_mm0./alphas_Ekin;

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
% alphas_grad_Phi(INNER_PART)=lininterp3( grad_Phi_tor_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
alphas_grad_psi_star(INNER_PART)=lininterp3( grad_psi_star_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
alphas_Epot(INNER_PART)=interp2_omega_map(scale_psi(1:size_r),scale_theta,E_potential_omega_map,1,DTHETA,alphas_psi, alphas_omega,INNER_PART);
% alphas_Epot(INNER_PART)=lininterp3( E_potential_PR_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
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
psi_star_omega_map_half_step=psi_star_omega_map;
delta_ps_prev=delta_ps;

alphas_omega0=alphas_omega;
alphas_Rpos_int=alphas_Rpos;
alphas_vphi_grad_psi_star=v_phi.*alphas_grad_psi_star;
% half_delta_t_coef=100*(0.5*h*FAST_SAWTOOTH/tau_cr);

