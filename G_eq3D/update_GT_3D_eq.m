%Finds new velocity (n+1/2) 
% Needs vphi for cylindrical correction.
% First estimates with linear interpolation
% Then a step only for vphi for best approximate
% Lastly the final step for v (every component)

%% Find B-field
% direction of the B field at local positions

Bfield_X_eq=interp2_XZ(interp_x,interp_z,BpolX_initial_XZsmall_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
Bfield_Z_eq=interp2_XZ(interp_x,interp_z,BpolZ_initial_XZsmall_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
Bfield_phi_eq=interp2_XZ(interp_x,interp_z,Bphi_XZsmall_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);

if APPLY_RMP_FIELD==1
    [IL3D_1, IL3D_2, IL3D_3, IL3D_4, IL3D_5, IL3D_6, IL3D_7, IL3D_8, slopex, slopey, slopez] = ...
        build_3Dinterp_indexarrays(scale_phi, scale_theta, scale_psi, DPHI,DTHETA,DPSIVAL,alphas_pos_phi_wrapped, alphas_theta, alphas_psi);
    BR_tilde=lininterp3( BR_RMP,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
    BZ_tilde=lininterp3( BZ_RMP,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
    Bphi_tilde=lininterp3( Bphi_RMP,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
    Bfield_X=Bfield_X_eq+BR_tilde;
    Bfield_Z=Bfield_Z_eq+BZ_tilde;
    Bfield_phi=Bfield_phi_eq+Bphi_tilde;
else
    Bfield_X=Bfield_X_eq;
    Bfield_Z=Bfield_Z_eq;
    Bfield_phi=Bfield_phi_eq;
end

alphas_Bfield=sqrt(Bfield_X.^2+Bfield_Z.^2+Bfield_phi.^2);

%% Estimate v_phi (n+1/2)
v_phi_tilde=1.5*v_phi_prev-0.5*v_phi_prev_prev;

%% Cylindrical correction B-field
Bfield_Z_tilde=Bfield_Z-(mHe/(ZHe*eV))*v_phi_tilde./(alphas_Rpos); %New B_z

Bfield_corr=sqrt(Bfield_X.^2+Bfield_Z_tilde.^2+Bfield_phi.^2); % New |B|

% unit vectors
bZ_corr=Bfield_Z_tilde./Bfield_corr; 
bX_corr=Bfield_X./Bfield_corr;
bphi_corr=Bfield_phi./Bfield_corr;

%% Iteration of v_phi for cylindrical correction
% determinant for inversion matrix
alphas_Omega=(ZHe*eV/mHe)*Bfield_corr;
epsilon=0.5*h*alphas_Omega;
epsilon_sq=epsilon.^2;
alphas_determinant=1+(epsilon_sq);

% Matrix for iteration to next speed (no E)
M31=2*epsilon.*(bZ_corr)+2*epsilon_sq.*bphi_corr.*bX_corr;
M32=2*epsilon.*(-bX_corr)+2*epsilon_sq.*bphi_corr.*bZ_corr;
M33=(1-epsilon_sq)+2*epsilon_sq.*(bphi_corr).^2;

% find the next value of vphi
v_phi=M31.*v_X_prev+M32.*v_Z_prev+M33.*v_phi_prev;
v_phi=v_phi./alphas_determinant;

% improved vphi step value
v_phi_tilde=0.5*v_phi+0.5*v_phi_prev;

%% Cylindrical correction B-field (final)
Bfield_Z_tilde=Bfield_Z-(mHe/(ZHe*eV))*v_phi_tilde./(alphas_Rpos); %New B_z

Bfield_corr=sqrt(Bfield_X.^2+Bfield_Z_tilde.^2+Bfield_phi.^2); % New |B|

% unit vectors
bZ_corr=Bfield_Z_tilde./Bfield_corr; 
bX_corr=Bfield_X./Bfield_corr;
bphi_corr=Bfield_phi./Bfield_corr;

%% Matrix for v step (n+1/2)
% determinant for inversion matrix
alphas_Omega=(ZHe*eV/mHe)*Bfield_corr;
epsilon=0.5*h*alphas_Omega;
epsilon_sq=epsilon.^2;
alphas_determinant=1+(epsilon_sq);

% Matrix for iteration to next speed (no E)
M11=(1-epsilon_sq)+2*epsilon_sq.*(bX_corr).^2;
M12=2*epsilon.*(bphi_corr)+2*epsilon_sq.*bX_corr.*bZ_corr;
M13=2*epsilon.*(-bZ_corr)+2*epsilon_sq.*bX_corr.*bphi_corr;

M21=2*epsilon.*(-bphi_corr)+2*epsilon_sq.*bZ_corr.*bX_corr;
M22=(1-epsilon_sq)+2*epsilon_sq.*(bZ_corr).^2;
M23=2*epsilon.*(bX_corr)+2*epsilon_sq.*bZ_corr.*bphi_corr;

M31=2*epsilon.*(bZ_corr)+2*epsilon_sq.*bphi_corr.*bX_corr;
M32=2*epsilon.*(-bX_corr)+2*epsilon_sq.*bphi_corr.*bZ_corr;
M33=(1-epsilon_sq)+2*epsilon_sq.*(bphi_corr).^2;

%% Update v (n+1/2)
% get the pure Bfield part of motion
v_X=M11.*v_X_prev+M12.*v_Z_prev+M13.*v_phi_prev;
v_Z=M21.*v_X_prev+M22.*v_Z_prev+M23.*v_phi_prev;
v_phi=M31.*v_X_prev+M32.*v_Z_prev+M33.*v_phi_prev;

% next velocity half time step values
v_X=v_X./alphas_determinant;
v_Z=v_Z./alphas_determinant;
v_phi=v_phi./alphas_determinant;