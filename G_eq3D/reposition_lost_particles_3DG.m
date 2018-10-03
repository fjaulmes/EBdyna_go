
% we only give them a new position (further away from magnetic axis)
% and keep their initial energy
alphas_pos_x(outcast)=max(min(X0(recast)-(Raxis-R0),scale_X(end-3)),scale_X(4));
alphas_pos_z(outcast)=max(min(Z0(recast),scale_Z(end-3)),scale_Z(4));
alphas_pos_phi(outcast)=rand(size(recast))*2*pi;


%  evaluate quickly new Bfield;
alphas_Bfield=interp2(scale_X,scale_Z,Btot_XZ_map',alphas_pos_x,alphas_pos_z,'*linear');

alphas_Epll(outcast)=0.5*(mHe/eV)*(alphas_vpll(outcast)).^2;
alphas_Eperp(outcast)=max(alphas_Ekin(outcast)-alphas_Epll(outcast),0);
alphas_mm_part(outcast)=alphas_Eperp(outcast)./alphas_Bfield(outcast);
alphas_vperp(outcast)=sqrt(2*alphas_Eperp(outcast)*eV/mHe);


% normal vector N = (cos) u + (sin) w
% from reused initial positions
NX(outcast)=cos(norm_angle(recast)).*uX(recast)+sin(norm_angle(recast)).*wX(recast);
NZ(outcast)=cos(norm_angle(recast)).*uZ(recast)+sin(norm_angle(recast)).*wZ(recast);
Nphi(outcast)=cos(norm_angle(recast)).*uphi(recast)+sin(norm_angle(recast)).*wphi(recast);

v_X(outcast)=alphas_vpll(outcast).*bX(outcast)+alphas_vperp(outcast).*NX(outcast);
v_Z(outcast)=alphas_vpll(outcast).*bZ(outcast)+alphas_vperp(outcast).*NZ(outcast);
v_phi(outcast)=alphas_vpll(outcast).*bphi(outcast)+alphas_vperp(outcast).*Nphi(outcast);

% taking the previous speed as initial speed
v_X_prev(outcast)=v_X(outcast);
v_Z_prev(outcast)=v_Z(outcast);
v_phi_prev(outcast)=v_phi(outcast);


interpolate_theta_psi_fromXZ;

