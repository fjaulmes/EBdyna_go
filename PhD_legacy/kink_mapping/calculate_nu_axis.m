nu_axis=zeros(NR,2);
psi_star_3=zeros(NR,2);
z1_axis=zeros(NR,2);
z2_axis=zeros(NR,2);


NU_OFFSET=0;

for (r=1:NR)
    r_value=scale_r(r);
    rv2=r_value^2;
    omega=1;
    omega_value=scale_omega(omega);
    z1_axis(r,1)=a1*((rv2-rx^2)*cos(alpha)+rx*(rx-sign_ksi0*r_value*cos(omega_value))*sin(alpha));
    z2_axis(r,1)=a2*(rx^2-rv2);
    omega=Nomega;
    omega_value=scale_omega(omega);
    z1_axis(r,2)=a1*((rv2-rx^2)*cos(alpha)+rx*(rx-sign_ksi0*r_value*cos(omega_value))*sin(alpha));
    z2_axis(r,2)=a2*(rx^2-rv2);
    
end
den=z1_axis+z2_axis+NU_OFFSET;
nu_axis=z1_axis.*z2_axis./(den);
nu_axis(isnan(nu_axis))=0;