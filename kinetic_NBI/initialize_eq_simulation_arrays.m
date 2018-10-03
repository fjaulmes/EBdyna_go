Nradial=length(psi_scale);
SIMULATION_RADIAL_LIMIT=Nradial-2


load(DISTNAME);
disp('verifying kinetic energy consistency:')
max(alphas_Ekin)
			
if ~exist('alphas_momentum')
    alphas_momentum=alphas_Ekin*0;
else
    alphas_momentum=alphas_Ekin*0+mean(alphas_momentum);
	disp('changing toroidal momentum (rad/s) to average value:')
	disp(mean(alphas_momentum))
end

alphas_Bfield=interp2(scale_X,scale_Z,Btot_XZ_map',alphas_pos_x,alphas_pos_z,'*linear');

alphas_mm=(alphas_Ekin-0.5*((mHe/eV)*alphas_vpll.^2))./alphas_Bfield;
alphas_ejected=zeros(Nalphas_simulated,1);
alphas_Eperp=alphas_Bfield.*alphas_mm;



Bavg=mean(Btot_XZ_map(round(0.8*mid_X),:));



alphas_Epll=0.5*(mHe/eV)*alphas_vpll.^2;
%alphas_Eperp=abs(alphas_Ekin-alphas_Epll);
alphas_Ekin=max(alphas_Ekin,alphas_Epll+alphas_Eperp);


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

phipos_output=zeros(NB_TIME_STAMPS,Nalphas_simulated);
% theta_output=zeros(NB_TIME_STAMPS,Nalphas_simulated);
vpll_output=zeros(NB_TIME_STAMPS,Nalphas_simulated);
Xpos_gc_output=zeros(NB_TIME_STAMPS,Nalphas_simulated);
Zpos_gc_output=zeros(NB_TIME_STAMPS,Nalphas_simulated);


alphas_psi=zeros(Nalphas_simulated,1);
alphas_theta=zeros(Nalphas_simulated,1);
bX=zeros(Nalphas_simulated,1);
bZ=zeros(Nalphas_simulated,1);
bphi=zeros(Nalphas_simulated,1);
alphas_psi_value=zeros(Nalphas_simulated,1);
alphas_psi_value_corr=zeros(Nalphas_simulated,1);




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

alphas_Omega=(ZHe*eV/mHe)*alphas_Bfield;
epsilon=0.5*h*alphas_Omega;
alphas_determinant=1+(epsilon.^2);

% previous mid point speed values
v_X_prev=v0_X+epsilon.*(-bphi.*v0_Z+bZ.*v0_phi);
v_Z_prev=v0_Z+epsilon.*(bphi.*v0_X-bX.*v0_phi);
v_phi_prev=v0_phi+epsilon.*(-bZ.*v0_X+bX.*v0_Z);

v_phi_prev_prev=v_phi_prev;
% v_phi_prev=v_phi;

%initialize Bfield properly
alphas_mm_part=alphas_mm;
% update_GT_3D_collapse;

interpolate_theta_psi_fromXZ;

% amplitude of the B field and potential at half time step local positions
alphas_Bfield=interp2_XZ(interp_x,interp_z,Btot_XZ_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);

% Canonical angular momentum evolution
alphas_omega=wrap2pi(alphas_theta-alphas_pos_phi_wrapped);

alphas_psi_value_corr=alphas_psi_value;

alphas_pphi0=(mHe/eV)*alphas_Rpos.*v_phi-(ZHe)*alphas_psi_value_corr;



%trapping parameter
%radial_pos=(a/257)*interp2(scale_X,scale_Z,radial_XZsmall_map',alphas_pos_x,alphas_pos_z,'*linear');
Bavg=mean(alphas_Bfield);
%alphas_kappa=sqrt((alphas_Ekin*R0+Bavg*alphas_mm.*(radial_pos-R0))./(2*alphas_mm.*radial_pos*Bavg));
alphas_lambda=Bavg*alphas_mm./alphas_Ekin;

alphas_Eperp=alphas_Ekin-alphas_Epll;
alphas_vperp=sqrt(2*alphas_Eperp*eV/mHe);
alphas_rhoL=(mHe/eV)*(1/ZHe)*alphas_vperp./alphas_Bfield;




alphas_Ekin_half=alphas_Ekin;
alphas_Ekin_prev=alphas_Ekin;
alphas_pphi0_prev=alphas_pphi0;



interpolate_theta_psi_fromXZ;
