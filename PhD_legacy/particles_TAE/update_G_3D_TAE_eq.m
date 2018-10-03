
% cylindrical correction to the magnetic field

v_phi_tilde=1.5*v_phi_prev-0.5*v_phi_prev_prev;

Bfield_corr=alphas_Bfield;
bZ_corr=bZ;
bX_corr=bX;
bphi_corr=bphi;

Bfield_Z=bZ.*alphas_Bfield;
Bfield_Z=Bfield_Z-(mHe/(ZHe*eV))*v_phi_tilde./(alphas_Rpos);
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

%M11=(1-epsilon_sq)+2*epsilon_sq.*(bX_corr).^2;
%M12=2*epsilon.*(bphi_corr)+2*epsilon_sq.*bX_corr.*bZ_corr;
%M13=2*epsilon.*(-bZ_corr)+2*epsilon_sq.*bX_corr.*bphi_corr;

%M21=2*epsilon.*(-bphi_corr)+2*epsilon_sq.*bZ_corr.*bX_corr;
%M22=(1-epsilon_sq)+2*epsilon_sq.*(bZ_corr).^2;
%M23=2*epsilon.*(bX_corr)+2*epsilon_sq.*bZ_corr.*bphi_corr;

M31=2*epsilon.*(bZ_corr)+2*epsilon_sq.*bphi_corr.*bX_corr;
M32=2*epsilon.*(-bX_corr)+2*epsilon_sq.*bphi_corr.*bZ_corr;
M33=(1-epsilon_sq)+2*epsilon_sq.*(bphi_corr).^2;
% 

% get the pure Bfield part of motion

%v_X=M11.*v_X_prev+M12.*v_Z_prev+M13.*v_phi_prev;
%v_Z=M21.*v_X_prev+M22.*v_Z_prev+M23.*v_phi_prev;
v_phi=M31.*v_X_prev+M32.*v_Z_prev+M33.*v_phi_prev;


% v_X=v_X./alphas_determinant;
% v_Z=v_Z./alphas_determinant;
v_phi_tilde=0.5*(v_phi_prev+v_phi./alphas_determinant);




% now complete itereation with most accurate time step value

Bfield_corr=alphas_Bfield;
bZ_corr=bZ;
bX_corr=bX;
bphi_corr=bphi;

Bfield_Z=bZ.*alphas_Bfield;
Bfield_Z=Bfield_Z-(mHe/(ZHe*eV))*v_phi_tilde./(alphas_Rpos);
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

v_X=M11.*v_X_prev+M12.*v_Z_prev+M13.*v_phi_prev;
v_Z=M21.*v_X_prev+M22.*v_Z_prev+M23.*v_phi_prev;
v_phi=M31.*v_X_prev+M32.*v_Z_prev+M33.*v_phi_prev;


v_X=v_X./alphas_determinant;
v_Z=v_Z./alphas_determinant;
v_phi=v_phi./alphas_determinant;

