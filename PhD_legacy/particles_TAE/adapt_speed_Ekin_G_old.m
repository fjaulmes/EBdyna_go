% alphas_Ekin_half=0.5*alphas_Ekin+0.5*alphas_Ekin_prev;
if (~exist('v_X_next'))
  v_X_next=v_X;
  v_Z_next=v_Z;
  v_phi_next=v_phi;
end
if (~exist('v_X_prev'))
  v_X_prev=v_X;
  v_Z_prev=v_Z;
  v_phi_prev=v_phi;
end
if (~exist('v_X_step'))
  v_X_step=v_X;
  v_Z_step=v_Z;
  v_phi_step=v_phi;
end

alphas_vtot_sq=(2*alphas_Ekin*eV/mHe);
vtot_recalc_sq=(v_X_step.^2+v_Z_step.^2+v_phi_step.^2);%
Ekin_corr=0.5*(1+sqrt(alphas_vtot_sq./vtot_recalc_sq));
% Ekin_corr=sqrt(alphas_vtot_sq./vtot_recalc_sq);

%half point corrected values
%adjusting the next time step value
v_X_next=v_X_next.*(Ekin_corr);
v_Z_next=v_Z_next.*(Ekin_corr);
v_phi_next=v_phi_next.*(Ekin_corr);
% v_phi=0.5*(v_phi+v_phi.*(Ekin_corr));

%adjusting the time step value
alphas_vtot_sq=(2*alphas_Ekin*eV/mHe);
vtot_recalc_sq=(v_X_step.^2+v_Z_step.^2+v_phi_step.^2);%
Ekin_corr=0.5*(1+sqrt(alphas_vtot_sq./vtot_recalc_sq));
% Ekin_corr=sqrt(alphas_vtot_sq./vtot_recalc_sq);

v_X=v_X.*(Ekin_corr);
v_Z=v_Z.*(Ekin_corr);
v_phi=v_phi.*(Ekin_corr);

v_X_step=0.5*(v_X+v_X_next);
v_Z_step=0.5*(v_Z+v_Z_next);
v_phi_step=0.5*(v_phi+v_phi_next);

