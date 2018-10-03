% alphas_Ekin_half=0.5*alphas_Ekin+0.5*alphas_Ekin_prev;
alphas_vtot_sq=(2*alphas_Ekin*eV/mHe);
vtot_recalc_sq=(v_X.^2+v_Z.^2+v_phi.^2);%
vtot_recalc_sq_prev=vtot_recalc_sq;

%delta_Ekin=ECOEF_P*(alphas_vtot_sq-vtot_recalc_sq);
%vtot_integ=vtot_integ+ECOEF_INTEG*delta_Ekin;   
%vtot_recalc_sq=vtot_recalc_sq+delta_Ekin+vtot_integ;

Ekin_corr=(1-ECOEF)+(ECOEF)*sqrt(alphas_vtot_sq./vtot_recalc_sq_prev);
% Ekin_corr=sqrt(alphas_vtot_sq./vtot_recalc_sq);

%half point corrected values
v_X=v_X.*(Ekin_corr);
v_Z=v_Z.*(Ekin_corr);
v_phi=v_phi.*(Ekin_corr);

v_X_prev=v_X;
v_Z_prev=v_Z;
v_phi_prev=v_phi;

%v_X_prev=v_X_prev.*(Ekin_corr);
%v_X_prev=v_X_prev.*(Ekin_corr);
%v_X_prev=v_X_prev.*(Ekin_corr);

%v_X_step=(0.5*v_X+0.5*v_X_prev);
%v_Z_step=(0.5*v_Z+0.5*v_Z_prev);
%v_phi_step=(0.5*v_phi+0.5*v_phi_prev);