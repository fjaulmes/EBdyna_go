
% f1=(ksi_dot_recalc-rx_dot_recalc_evol)./ksi_dot_recalc;
% f2=rx_dot_recalc_evol./ksi_dot_recalc;

% f1=u_sep./ksi_dot;
% f2=rx_dot./ksi_dot;
% f1=f1.*min(Delta_SP./(2*pi*r_core),1);
% f2=f2.*Delta_SP./(2*pi*rx_evol_lin);
% 
% 
% f1(TRANSITION_FRAME)=(0.3*f1(TRANSITION_FRAME-1)+0.7*f1(TRANSITION_FRAME+1));
% f2(TRANSITION_FRAME)=(0.3*f2(TRANSITION_FRAME-1)+0.7*f2(TRANSITION_FRAME+1));
% f1(TRANSITION_FRAME-1)=0.5*(f1(TRANSITION_FRAME)+f1(TRANSITION_FRAME-2));
% f2(TRANSITION_FRAME-1)=0.5*(f2(TRANSITION_FRAME)+f2(TRANSITION_FRAME-2));

E1(1:TRANSITION_FRAME)=interp1(ksi0_evol_lin_ini(1:TRANSITION_FRAME),E1ini(1:TRANSITION_FRAME),min(ksi_recalc_evol(1:TRANSITION_FRAME),0.99*rmix));
E2(1:TRANSITION_FRAME)=interp1(ksi0_evol_lin_ini(1:TRANSITION_FRAME),E2ini(1:TRANSITION_FRAME),min(ksi_recalc_evol(1:TRANSITION_FRAME),0.99*rmix));
Efinal(1:TRANSITION_FRAME)=interp1(ksi0_evol_lin_ini(1:TRANSITION_FRAME),Efinalini(1:TRANSITION_FRAME),min(ksi_recalc_evol(1:TRANSITION_FRAME),0.99*rmix));
volume3_recalc(1:TRANSITION_FRAME)=interp1(ksi0_evol_lin_ini(1:TRANSITION_FRAME),volume3_evol_lin(1:TRANSITION_FRAME),min(ksi_recalc_evol(1:TRANSITION_FRAME),0.99*rmix));
dvolume3_recalc=gradient(volume3_recalc,1:length(volume1_evol_lin));

deltaE_prev=deltaE_vA;
deltaE=0.5*(E1+E2)-Efinal;

Bstar_diff_evol(1:TRANSITION_FRAME)=interp1(ksi0_evol_lin_ini(1:TRANSITION_FRAME),Bstar1_evol_lin(1:TRANSITION_FRAME)-Bstar3_evol_lin(1:TRANSITION_FRAME),min(ksi_recalc_evol(1:TRANSITION_FRAME),0.99*rmix));

jpll_evol=Bstar_diff_evol./(2*mu0*delta_evol);
joule_heating_new=eta1*(jpll_evol.^2.*(2*pi*R0*delta_evol*delta_e));
joule_heating_new=min(deltaE.*dvolume3_recalc,joule_heating_new);
joule_heating=0.5*joule_heating_new+0.5*joule_heating;
joule_heating=min(deltaE.*dvolume3_recalc,joule_heating);
% joule_heating=joule_heating_new;
deltaE_joules=max(deltaE.*dvolume3_recalc-joule_heating,0);
deltaE_new=deltaE_joules./dvolume3_recalc;
deltaE_vA=0.95*deltaE_new+0.05*deltaE_prev;

% deltaE(1:TRANSITION_FRAME)=max((f1(1:length(E1)).*E1(1:TRANSITION_FRAME)+f2(1:length(E1)).*E2(1:TRANSITION_FRAME))-Efinal(1:TRANSITION_FRAME),0);
% deltaE(1:TRANSITION_FRAME)=(f1(1:length(E1)).*E1(1:TRANSITION_FRAME)+f2(1:length(E1)).*E2(1:TRANSITION_FRAME));%-Efinal(1:TRANSITION_FRAME),0);
% deltaE=0.5*(deltaE+deltaE_prev);

% the volume energy is spli in two for the velocity
vA_out=sqrt(2*deltaE_vA)/sqrt(Ni1*mbulk);
% vA_out=cos(delta_theta0(1:TRANSITION_FRAME)-thetac_evol).*vA_out;


% smoothing vA?

% vA_out_error=(Delta_SP.*ksi_dot_recalc(1:TRANSITION_FRAME)./(2*delta_evol))-vA_out;
% vA_out_error_integ=vA_out_error_integ+0.05*vA_out_error;
% vA_out=vA_out+0.05*vA_out_error+0.02*vA_out_error_integ;


% vA_out=0.5*(vA_out+Delta_SP.*ksi_dot_recalc(1:TRANSITION_FRAME)./(2*delta_evol));
%correct the last point
vA_out(end)=vA_out(end-1);

% vA_out=0.5*(vA_out+0.5*Delta_SP.*ksi_dot_recalc(1:TRANSITION_FRAME)./delta_evol);

%vA_out(TRANSITION_FRAME-3:TRANSITION_FRAME)=0.5*(vA_out(TRANSITION_FRAME-3:TRANSITION_FRAME)+vA_out_ini(TRANSITION_FRAME-3:TRANSITION_FRAME));

% tau_star=0.5*Delta_SP./vA_out;
% 
% recalculate_layer_width;
% vA_out=0.5*(vA_out+Delta_SP.*ksi_dot_recalc(1:TRANSITION_FRAME)./(2*delta_evol));
% vA_out(end)=vA_out(end-1);
% 
tau_star=0.5*(Delta_SP./vA_out+2*delta_evol./ksi_dot_recalc(1:length(vA_out)));
% tau_star=2*Delta_SP./vA_out;
recalculate_layer_width;

% vA_out=0.5*(vA_out+Delta_SP.*ksi_dot_recalc(1:TRANSITION_FRAME)./(2*delta_evol));
% vA_out(end)=vA_out(end-1);
% tau_star=0.5*Delta_SP./vA_out;
% recalculate_layer_width;
% 
% vA_out=0.5*(vA_out+Delta_SP.*ksi_dot_recalc(1:TRANSITION_FRAME)./(2*delta_evol));
% vA_out(end)=vA_out(end-1);
% tau_star=0.5*Delta_SP./vA_out;
% recalculate_layer_width;
% 
% vA_out=0.5*(vA_out+Delta_SP.*ksi_dot_recalc(1:TRANSITION_FRAME)./(2*delta_evol));
% vA_out(end)=vA_out(end-1);
% tau_star=0.5*Delta_SP./vA_out;
% recalculate_layer_width;
