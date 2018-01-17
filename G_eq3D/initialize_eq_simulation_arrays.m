

% load the distribution file with all particles
load(DISTNAME);

alphas_Ekin=real(alphas_Ekin);

% splitting equally the particles among the processors

% can take a subdivision of the total
FRACTION_POP=1;

PARTICLES_SPLIT=floor(Nalphas_simulated/NB_PROCESS)
PART_POP_JOB=(1:NB_PROCESS*PARTICLES_SPLIT);
Nalphas_simulated=length(PART_POP_JOB)

PART_POP_JOB=PART_POP_JOB((PROCESS_NUMBER-1)*PARTICLES_SPLIT+1:FRACTION_POP:(PROCESS_NUMBER)*PARTICLES_SPLIT);
alphas_pos_x=alphas_pos_x(PART_POP_JOB);
alphas_pos_z=alphas_pos_z(PART_POP_JOB);
alphas_pos_phi=alphas_pos_phi(PART_POP_JOB);
alphas_Ekin=alphas_Ekin(PART_POP_JOB);
alphas_vpll=alphas_vpll(PART_POP_JOB);

if exist('v_X')
    v_X=v_X(PART_POP_JOB);
    v_Z=v_Z(PART_POP_JOB);
    v_phi=v_phi(PART_POP_JOB);
end


alphas_Bfield=interp2(scale_X,scale_Z,Btot_XZ_map',alphas_pos_x,alphas_pos_z,'*cubic');
alphas_mm=(alphas_Ekin-0.5*((mHe/eV)*alphas_vpll.^2))./alphas_Bfield;
alphas_mm=max(alphas_mm,0);

Nalphas_simulated=length(alphas_vpll)


disp('verifying kinetic energy consistency:')
max(alphas_Ekin)
			
if ~exist('alphas_momentum')
    alphas_momentum=alphas_Ekin*0;
	ACTIVATE_CENTRIFUG_TERM=0
else
    alphas_momentum=alphas_momentum(PART_POP_JOB);
	if sum(alphas_momentum)~=0
		ACTIVATE_CENTRIFUG_TERM=1
	else
		ACTIVATE_CENTRIFUG_TERM=0
	end
end
if ~exist('alphas_current')
    alphas_current=alphas_Ekin*0;
else
    alphas_current=alphas_current(PART_POP_JOB);
end
if ~exist('alphas_weight')
    alphas_weight=alphas_Ekin*0+1;
else
    alphas_weight=alphas_weight(PART_POP_JOB);
end

% enforcing double precision
% seems there was an issue with some of the data file

alphas_Ekin=double(alphas_Ekin);

alphas_mm=double(alphas_mm);

alphas_pos_x=double(alphas_pos_x);
alphas_pos_z=double(alphas_pos_z);
alphas_pos_phi=double(alphas_pos_phi);
alphas_vpll=double(alphas_vpll);

alphas_momentum=double(alphas_momentum);
alphas_current=double(alphas_current);
alphas_weight=double(alphas_weight);



alphas_ejected=zeros(Nalphas_simulated,1);
alphas_Eperp=alphas_Bfield.*alphas_mm;






alphas_Epll=0.5*(mHe/eV)*alphas_vpll.^2;
alphas_Ekin=max(alphas_Ekin,alphas_Epll+alphas_Eperp);


%initial position values
alphas_Ekin0=alphas_Ekin;
alphas_mm0=alphas_mm;
X0=alphas_pos_x;
Z0=alphas_pos_z;
phi0=alphas_pos_phi;
alphas_Rpos=alphas_pos_x;
alphas_Rpos_int=alphas_Rpos;

alphas_Eperp=max(alphas_Ekin-alphas_Epll,0);
alphas_vperp=sqrt(2*alphas_Eperp*eV/mHe);
alphas_rhoL=(mHe/eV)*(1/ZHe)*alphas_vperp./alphas_Bfield;


disp('**************************************************************');

% initializing simulation arrays
% outputs recorded in the precession file

phipos_output=zeros(NB_TIME_STAMPS,Nalphas_simulated);
% theta_output=zeros(NB_TIME_STAMPS,Nalphas_simulated);
vpll_output=zeros(NB_TIME_STAMPS,Nalphas_simulated);
Xpos_gc_output=zeros(NB_TIME_STAMPS,Nalphas_simulated);
Zpos_gc_output=zeros(NB_TIME_STAMPS,Nalphas_simulated);


alphas_psi=zeros(Nalphas_simulated,1);
alphas_theta=zeros(Nalphas_simulated,1);
bX=zeros(Nalphas_simulated,1);
bX_corr=zeros(Nalphas_simulated,1);
bX_coriolis=zeros(Nalphas_simulated,1);
bX_tot=zeros(Nalphas_simulated,1);
bZ=zeros(Nalphas_simulated,1);
bphi=zeros(Nalphas_simulated,1);
alphas_psi_value=zeros(Nalphas_simulated,1);
alphas_psi_value_corr=zeros(Nalphas_simulated,1);
if CALCULATE_TRUE_PPHI==1
    alphas_dAphi_dphi=zeros(Nalphas_simulated,1);
end


% initial field and position
alphas_psi_value=interp2(scale_X,scale_Z,psi_XZsmall_map',alphas_pos_x,alphas_pos_z,'*linear');
alphas_psi_value=max(alphas_psi_value,-psi_global);
alphas_psi=interp1(psi_scale,1:Nradial,alphas_psi_value);

alphas_theta=interp2(scale_X,scale_Z,theta_XZsmall_map',alphas_pos_x,alphas_pos_z,'*nearest');



time=0;





% estimating psi and theta values
interpolate_theta_psi_fromXZ;
alphas_pos_phi_wrapped=wrap2pi(alphas_pos_phi);

% direction of the B field at local positions
bX=interp2_XZ(interp_x,interp_z,bX_XZ_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
bZ=interp2_XZ(interp_x,interp_z,bZ_XZ_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
bphi=TOROIDAL_FIELD_DIRECTION*sqrt(1-(bX.^2+bZ.^2));
	
% attribute gyro phase if not included in initial data
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

if ~exist('v_X')

	v0_X=alphas_vpll.*bX+alphas_vperp.*NX;
	v0_Z=alphas_vpll.*bZ+alphas_vperp.*NZ;
	v0_phi=alphas_vpll.*bphi+alphas_vperp.*Nphi;

	v0_X=double(v0_X);
	v0_Z=double(v0_Z);
	v0_phi=double(v0_phi);
		
	v_X=v0_X;
	v_Z=v0_Z;
	v_phi=v0_phi;
else
	v0_X=v_X;
	v0_Z=v_Z;
	v0_phi=v_phi;
end


alphas_Omega=(ZHe*eV/mHe)*alphas_Bfield;
epsilon=0.5*h*alphas_Omega;
alphas_determinant=1+(epsilon.^2);

% previous mid point speed values
v_X_prev=v0_X+epsilon.*(-bphi.*v0_Z+bZ.*v0_phi);
v_Z_prev=v0_Z+epsilon.*(bphi.*v0_X-bX.*v0_phi);
v_phi_prev=v0_phi+epsilon.*(-bZ.*v0_X+bX.*v0_Z);

v_phi_prev_prev=v_phi_prev;

%initialize Bfield properly
alphas_mm_part=alphas_mm;
% update_GT_3D_collapse;

interpolate_theta_psi_fromXZ;

% amplitude of the B field and potential at half time step local positions
alphas_Bfield=interp2_XZ(interp_x,interp_z,Btot_XZ_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
if APPLY_RMP_FIELD==1
    BR_tilde=alphas_Bfield*0;
    BZ_tilde=alphas_Bfield*0;
    Bphi_tilde=alphas_Bfield*0;
    [IL3D_1 IL3D_2 IL3D_3 IL3D_4 IL3D_5 IL3D_6 IL3D_7 IL3D_8 slopex slopey slopez] = ...
        build_3Dinterp_indexarrays(scale_phi, scale_theta, scale_psi, DPHI,DTHETA,DPSIVAL,alphas_pos_phi_wrapped, alphas_theta, alphas_psi);
    BR_tilde=lininterp3( BR_RMP,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
    BZ_tilde=lininterp3( BZ_RMP,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
    Bphi_tilde=lininterp3( Bphi_RMP,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
    if CALCULATE_TRUE_PPHI==1
        alphas_dAphi_dphi=alphas_mm*0;
        alphas_Delta_pphi=alphas_mm*0;
        alphas_dAphi_dphi=lininterp3( dAphi_dphi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
    end
end


alphas_omega=wrap2pi(alphas_theta-alphas_pos_phi_wrapped);

alphas_psi_value_corr=alphas_psi_value;

% Canonical angular momentum evolution
alphas_pphi0=(mHe/eV)*alphas_Rpos.*v_phi-(ZHe)*alphas_psi_value_corr;



%trapping parameter
alphas_lambda=Bavg*alphas_mm./alphas_Ekin;

alphas_Eperp=alphas_Ekin-alphas_Epll;
alphas_vperp=sqrt(2*alphas_Eperp*eV/mHe);
alphas_rhoL=(mHe/eV)*(1/ZHe)*alphas_vperp./alphas_Bfield;



alphas_Ekin_half=alphas_Ekin;
alphas_Ekin_prev=alphas_Ekin;
alphas_pphi0_prev=alphas_pphi0;



interpolate_theta_psi_fromXZ;
