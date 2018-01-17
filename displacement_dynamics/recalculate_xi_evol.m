% ksi_dot_recalc=zeros(1,length(time_scale_lin));



DELTA_TIME=(time_scale_lin(2)-time_scale_lin(1))*crash_duration;
% ksi_dot_recalc=gradient(ksi0_evol_lin_ini,time_scale_lin*crash_duration);
ksi_dot_recalc(1:TRANSITION_FRAME)=g_xi(1:TRANSITION_FRAME);

% g_value=interp1(ksi0_evol_lin_ini(1:TRANSITION_FRAME),g_xi,rmix);
% ksi_dot_recalc(TRANSITION_FRAME)=(g_value);
ksi_recalc_evol(TRANSITION_FRAME)=0.5*(rmix+max_ksi0);
ksi_recalc_evol(TRANSITION_FRAME-1)=max(ksi_recalc_evol(TRANSITION_FRAME)-ksi_dot_recalc(TRANSITION_FRAME)*DELTA_TIME,0);
% g_value=interp1(ksi0_evol_lin_ini(1:TRANSITION_FRAME),g_xi,ksi_recalc_evol(TRANSITION_FRAME-1));
% ksi_dot_recalc(TRANSITION_FRAME-1)=(g_value);
% we need to integrate one time with this time step size
% and see if it fits

for t=TRANSITION_FRAME:-1:2
    % approximate derivative but
    % we are not in less than a mm precision anyway....
    ksi_recalc_evol(t-1)=max(ksi_recalc_evol(t)-0.5*(ksi_dot_recalc(t)+ksi_dot_recalc(t-1))*DELTA_TIME,0);
%     if ksi_recalc_evol(t-1)>0.01
%         g_value=interp1(ksi0_evol_lin_ini(1:TRANSITION_FRAME),g_xi,ksi_recalc_evol(t-1));
%     end
%     ksi_dot_recalc(t-1)=g_value;
end

xi_time_ini=ksi_recalc_evol(1)

if xi_time_ini>XI0_INI
	while xi_time_ini>XI0_INI
        % velocity kept tconstant
        % so we want a longer crash to go further
		crash_duration=min(crash_duration+0.001*tau_cr,8*tau_cr);
%         ksi_dot_recalc=gradient(ksi0_evol_lin_ini,time_scale_lin*crash_duration);

%         g_xi=g_xi.*(1+0.01*ksi0_evol_lin_ini(1:TRANSITION_FRAME)/max(ksi0_evol_lin_ini(1:TRANSITION_FRAME)));
%         g_xi=g_xi.*(1.005);
		DELTA_TIME=(time_scale_lin(2)-time_scale_lin(1))*crash_duration;
		for t=TRANSITION_FRAME:-1:2
			% this is a bit tedious but we capture here the full variations of g
			ksi_recalc_evol(t-1)=max(ksi_recalc_evol(t)-0.5*(ksi_dot_recalc(t)+ksi_dot_recalc(t-1))*DELTA_TIME,0);
%             if ksi_recalc_evol(t-1)>0.01
%                 g_value=interp1(ksi0_evol_lin_ini(1:TRANSITION_FRAME),g_xi,ksi_recalc_evol(t-1));
%             end
% 			ksi_dot_recalc(t-1)=g_value;
		end
		xi_time_ini=ksi_recalc_evol(1);
    end
end
if xi_time_ini<=0.1*XI0_INI
	while xi_time_ini<=0.1*XI0_INI
		crash_duration=max(crash_duration-0.001*tau_cr,0.2*tau_cr);
%         ksi_dot_recalc=gradient(ksi0_evol_lin_ini,time_scale_lin*crash_duration);

%         g_xi=g_xi.*(0.995);
        DELTA_TIME=(time_scale_lin(2)-time_scale_lin(1))*crash_duration;
		for t=TRANSITION_FRAME:-1:2
			% this is a bit tedious but we capture here the full variations of g
			ksi_recalc_evol(t-1)=max(ksi_recalc_evol(t)-0.5*(ksi_dot_recalc(t)+ksi_dot_recalc(t-1))*DELTA_TIME,0);
%             if ksi_recalc_evol(t-1)>0.01
%                 g_value=interp1(ksi0_evol_lin_ini(1:TRANSITION_FRAME),g_xi,ksi_recalc_evol(t-1));
%             end
%             ksi_dot_recalc(t-1)=g_value;
		end
		xi_time_ini=ksi_recalc_evol(1);
	end
end


%%
% ksi_recalc_evol(1)=XI0_INI;

for t=1:TRANSITION_FRAME-1
	% to get the begining correct
	ksi_recalc_evol(t+1)=ksi_recalc_evol(t)+0.5*(ksi_dot_recalc(t)+ksi_dot_recalc(t+1))*DELTA_TIME;
end
if abs(ksi_recalc_evol(TRANSITION_FRAME)-rmix)>0.05
	disp('strange value of xi(transition)=')
	disp(ksi_recalc_evol(TRANSITION_FRAME))
	pause
end
% for t=TRANSITION_FRAME+1:TRANSITION_FRAME_XI
% 	ksi_dot_recalc(t)=(max(ksi0_evol_lin_ini)-ksi_recalc_evol(TRANSITION_FRAME))/((TRANSITION_FRAME_XI-TRANSITION_FRAME)*DELTA_TIME);
% end
for t=TRANSITION_FRAME:TRANSITION_FRAME_XI
	% to get the end correct
	ksi_recalc_evol(t+1)=ksi_recalc_evol(t)+0.5*(ksi_dot_recalc(t)+ksi_dot_recalc(t+1))*DELTA_TIME;
end
ksi_recalc_evol(TRANSITION_FRAME_XI-1:end)=0.5*(max(ksi0_evol_lin_ini)+ksi_recalc_evol(TRANSITION_FRAME_XI-1));

rx_recalc_evol=interp1(ksi0_evol_lin_ini(1:TRANSITION_FRAME_XI-1),rx_evol_lin_ini(1:TRANSITION_FRAME_XI-1),ksi_recalc_evol);
rx_recalc_evol(TRANSITION_FRAME:end)=max(rx_evol_lin_ini);
% rx_recalc_evol(1)=r_value_q1_mean;

r_core_recalc=rx_recalc_evol-ksi_recalc_evol;



% finally get back Delta from the most suited angle thetac
Delta_SP=0.5*(Delta_SP+2*(thetac_evol(1:length(vA_out))).*r_core_recalc(1:length(vA_out)));

%this first value makes little sense
Delta_SP(1)=Delta_SP(2);

crash_duration

rx_dot_recalc_evol=gradient(rx_recalc_evol,time_scale_lin*crash_duration);
rx_dot_recalc_evol=max(rx_dot_recalc_evol,0);

% eventually iterate again?
ksi0_evol_lin=ksi_recalc_evol;
r_core=max(r_core_recalc,0);
rx_evol_lin=rx_recalc_evol;
