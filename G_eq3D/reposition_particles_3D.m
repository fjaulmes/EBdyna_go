function  [x,v]=reposition_particles_3D(par,dim,maps,x,v,outcast)
% Keep energy, but give a new position
      recast=randi(size(x,1),size(outcast,1),1);
x(outcast,1)=max(min(X0(recast)-(Raxis-R0),scale_X(end-3)),scale_X(4));
x(outcast,2)=max(min(Z0(recast),scale_Z(end-3)),scale_Z(4));
rng('shuffle')
x(outcast,3)=rand(x(outcast,3))*2*pi;


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
end