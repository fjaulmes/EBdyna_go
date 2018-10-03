alphas_Ekin_half=0.5*alphas_Ekin+0.5*alphas_Ekin_prev;

alphas_vtot_sq=(2*alphas_Ekin_half*eV/mHe);
vtot_recalc_sq=(v_X.^2+v_Z.^2+v_phi.^2);%
Ekin_corr=0.5*(1+sqrt(alphas_vtot_sq./vtot_recalc_sq));
%Ekin_corr=sqrt(alphas_vtot_sq./vtot_recalc_sq);

%half point corrected values
v_X=v_X.*(Ekin_corr);
v_Z=v_Z.*(Ekin_corr);
v_phi=v_phi.*(Ekin_corr);
% v_phi=0.5*(v_phi+v_phi.*(Ekin_corr));

%adjusting the time step value
%v_X_prev=v_X;
%v_Z_prev=v_Z;
%v_phi_prev=v_phi;


