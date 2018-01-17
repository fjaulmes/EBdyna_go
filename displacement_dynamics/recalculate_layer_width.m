delta_evol=((C0/omegape)*(1./tau_star+(vthe1./l_out)).^(1/2)).*(tau_star.^(1/2));

%compromise for convergence
% delta_evol=0.5*(delta_evol+Delta_SP.*ksi_dot_recalc(1:TRANSITION_FRAME)./(2*vA_out));

delta_error=(Delta_SP.*ksi_dot_recalc(1:TRANSITION_FRAME)./(2*vA_out))-delta_evol;
delta_error_integ=delta_error_integ+0.05*delta_error;
delta_evol=delta_evol+0.5*delta_error+0.2*delta_error_integ;

delta_evol=max(delta_evol,delta_e);
delta_evol=min(delta_evol,MAX_DELTA_ALLOWED);

delta_avg=mean(delta_evol(8:end-2));

%assumption of non increasing layer width
% [min_delta min_pos ]=min(delta_evol);
% delta_evol(min_pos:end)=min_delta;
