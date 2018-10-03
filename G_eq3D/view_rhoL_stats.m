alphas_Epll=0.5*(mHe/eV)*alphas_vpll.^2;
alphas_Eperp=max(alphas_Ekin-alphas_Epll,0);
alphas_Bfield=interp2(scale_X,scale_Z,Btot_XZ_map',alphas_pos_x,alphas_pos_z);
alphas_vperp=sqrt(2*alphas_Eperp*eV/mHe);
alphas_rhoL=(mHe/eV)*(1/ZHe)*alphas_vperp./alphas_Bfield;

hold on
plot(alphas_lambda(COUNTER_PASSING),alphas_rhoL(COUNTER_PASSING)/rmix,'b.')
plot(alphas_lambda(ALL_TRAPPED),alphas_rhoL(ALL_TRAPPED)/rmix,'g.')