
% direction of the B field at local positions

bX=interp2_XZ(interp_x,interp_z,bX_XZ_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
bZ=interp2_XZ(interp_x,interp_z,bZ_XZ_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);

bphi=TOROIDAL_FIELD_DIRECTION*sqrt(1-(bX.^2+bZ.^2));

Bfield_corr=alphas_Bfield;
bZ_corr=bZ;
bX_corr=bX;
bphi_corr=bphi;


% cylindrical correction to the magnetic field

v_phi_next_tilde=1.5*v_phi-0.5*v_phi_prev;

% Bfield_corr=alphas_Bfield;
% bZ_corr=bZ;
% bX_corr=bX;
% bphi_corr=bphi;

Bfield_Z=bZ.*alphas_Bfield;
Bfield_Z=Bfield_Z-(mHe/(ZHe*eV))*v_phi_next_tilde./(alphas_Rpos);
Bfield_corr=alphas_Bfield.*(sqrt(bX.^2+(Bfield_Z./alphas_Bfield).^2+bphi.^2));
bZ_corr=Bfield_Z./Bfield_corr;
bX_corr=(bX.*alphas_Bfield)./Bfield_corr;
bphi_corr=(bphi.*alphas_Bfield)./Bfield_corr;
% 


alphas_Omega=(ZHe*eV/mHe)*Bfield_corr;
epsilon=0.5*h*alphas_Omega;
epsilon_sq=epsilon.^2;
alphas_determinant=1+(epsilon_sq);


% Matrix for iteration to next speed (no E)

M31=2*epsilon.*(bZ_corr)+2*epsilon_sq.*bphi_corr.*bX_corr;
M32=2*epsilon.*(-bX_corr)+2*epsilon_sq.*bphi_corr.*bZ_corr;
M33=(1-epsilon_sq)+2*epsilon_sq.*(bphi_corr).^2;

% get the pure Bfield part of motion

v_phi_next=M31.*v_X+M32.*v_Z+M33.*v_phi;


v_phi_next=v_phi_next./alphas_determinant;

% v_slope=(v_phi_next-v_phi_prev)/4;
% v_phi_next_tilde=v_phi+v_slope;
v_phi_next_tilde=0.5*(v_phi+v_phi_next);
% v_phi_next_tilde=calc_vstep_value(v_phi_prev,v_phi,v_phi_next);

Bfield_Z=bZ.*alphas_Bfield;
Bfield_Z=Bfield_Z-(mHe/(ZHe*eV))*v_phi_next_tilde./(alphas_Rpos);
Bfield_corr=alphas_Bfield.*(sqrt(bX.^2+(Bfield_Z./alphas_Bfield).^2+bphi.^2));
bZ_corr=Bfield_Z./Bfield_corr;
bX_corr=(bX.*alphas_Bfield)./Bfield_corr;
bphi_corr=(bphi.*alphas_Bfield)./Bfield_corr;






alphas_Omega=(ZHe*eV/mHe)*Bfield_corr;
epsilon=0.5*h*alphas_Omega;
epsilon_sq=epsilon.^2;
alphas_determinant=1+(epsilon_sq);


% Matrix for iteration to next speed (no E)

M11=(1-epsilon_sq)+2*epsilon_sq.*(bX_corr).^2;
%M11=M11./alphas_determinant;
M12=2*epsilon.*(bphi_corr)+2*epsilon_sq.*bX_corr.*bZ_corr;
%M12=M12./alphas_determinant;
M13=2*epsilon.*(-bZ_corr)+2*epsilon_sq.*bX_corr.*bphi_corr;
%M13=M13./alphas_determinant;

M21=2*epsilon.*(-bphi_corr)+2*epsilon_sq.*bZ_corr.*bX_corr;
%M21=M21./alphas_determinant;
M22=(1-epsilon_sq)+2*epsilon_sq.*(bZ_corr).^2;
%M22=M22./alphas_determinant;
M23=2*epsilon.*(bX_corr)+2*epsilon_sq.*bZ_corr.*bphi_corr;
%M23=M23./alphas_determinant;

M31=2*epsilon.*(bZ_corr)+2*epsilon_sq.*bphi_corr.*bX_corr;
%M31=M31./alphas_determinant;
M32=2*epsilon.*(-bX_corr)+2*epsilon_sq.*bphi_corr.*bZ_corr;
%M32=M32./alphas_determinant;
M33=(1-epsilon_sq)+2*epsilon_sq.*(bphi_corr).^2;
%M33=M33./alphas_determinant;
% 

% get the pure Bfield part of motion

v_X_next=M11.*v_X+M12.*v_Z+M13.*v_phi;
v_Z_next=M21.*v_X+M22.*v_Z+M23.*v_phi;
v_phi_next=M31.*v_X+M32.*v_Z+M33.*v_phi;


v_X_next=v_X_next./alphas_determinant;
v_Z_next=v_Z_next./alphas_determinant;
v_phi_next=v_phi_next./alphas_determinant;

