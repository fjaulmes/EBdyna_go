% ksi_dot_recalc=zeros(1,length(time_scale_lin));


DELTA_TIME=(time_scale_lin(2)-time_scale_lin(1))*crash_duration;
% ksi_dot_recalc=2*delta_evol.*vA_out_normal./(2*delta_theta0(1:length(vA_out)).*rx_evol_lin(1:length(vA_out)));
ksi_dot_recalc=2*delta_evol.*vA_out./Delta_SP(1:length(vA_out));
ksi_dot_recalc(1)=0;
ksi_dot_recalc(2)=0.5*ksi_dot_recalc(3);

% ksi_dot_recalc(TRANSITION_FRAME-2)=0.5*(ksi_dot_recalc(TRANSITION_FRAME-3));
% ksi_dot_recalc(TRANSITION_FRAME-1)=0.5*(ksi_dot_recalc(TRANSITION_FRAME-2));
ksi_dot_recalc(TRANSITION_FRAME+1)=0.5*(ksi_dot_recalc(TRANSITION_FRAME));

% ksi_dot_recalc=1.15*ksi_dot_recalc;
%ultimately, if ksidot is too large it neesds to be changed
% if (sum(ksi_dot_recalc(1:TRANSITION_FRAME+1))*DELTA_TIME+KSI0_INI)>0.99*rmix
%     disp('exceeded ksi limit')
%     ksi_dot_recalc=0.99*ksi_dot_recalc*rmix/(sum(ksi_dot_recalc(1:TRANSITION_FRAME+1))*DELTA_TIME+KSI0_INI);
% end

ksi_dot_recalc(TRANSITION_FRAME+2:length(time_scale_lin))=0;
ksi_recalc_evol=zeros(length(time_scale_lin),1)+KSI0_INI;
for t=2:length(time_scale_lin)
    ksi_recalc_evol(t)=ksi_recalc_evol(t-1)+ksi_dot_recalc(t-1)*DELTA_TIME;
end
crash_duration=crash_duration+0.01*crash_duration*(rmix-(sum(ksi_dot_recalc(1:TRANSITION_FRAME+1))*DELTA_TIME+KSI0_INI));

%adjust the crash duration so that we get to rmix!
% while max(ksi_recalc_evol(1:TRANSITION_FRAME))>rmix
%     DELTA_TIME=(time_scale_lin(2)-time_scale_lin(1))*crash_duration;
%     
%     ksi_recalc_evol=zeros(length(time_scale_lin),1)+KSI0_INI;
%     for t=2:length(time_scale_lin)
%         ksi_recalc_evol(t)=ksi_recalc_evol(t-1)+ksi_dot_recalc(t-1)*DELTA_TIME;
%     end
%     
%     crash_duration=crash_duration+0.01*crash_duration*(rmix-(sum(ksi_dot_recalc(1:TRANSITION_FRAME+1))*DELTA_TIME+KSI0_INI));
% 
% %     crash_duration=0.1*crash_duration+0.9*crash_duration*rmix/(max(ksi_recalc_evol(1:TRANSITION_FRAME)));
%     
% end

DELTA_TIME=(time_scale_lin(2)-time_scale_lin(1))*crash_duration;

ksi_recalc_evol=zeros(length(time_scale_lin),1)+KSI0_INI;
for t=2:length(time_scale_lin)
    ksi_recalc_evol(t)=ksi_recalc_evol(t-1)+ksi_dot_recalc(t-1)*DELTA_TIME;
end

%ultimately, if ksidot is too large it neesds to be changed
% if (sum(ksi_dot_recalc(1:TRANSITION_FRAME+1))*DELTA_TIME+KSI0_INI)>0.99*rmix
%     disp('exceeded ksi limit')
%     ksi_dot_recalc=0.99*ksi_dot_recalc*rmix/(sum(ksi_dot_recalc(1:TRANSITION_FRAME+1))*DELTA_TIME+KSI0_INI);
% end
% ksi_recalc_evol=zeros(length(time_scale_lin),1)+KSI0_INI;
% for t=2:length(time_scale_lin)
%     ksi_recalc_evol(t)=ksi_recalc_evol(t-1)+ksi_dot_recalc(t-1)*DELTA_TIME;
% end
% if max(ksi_recalc_evol(1:TRANSITION_FRAME))>rmix
%     crash_duration=0.5*crash_duration
% end
% 
% ksi_recalc_evol=zeros(length(time_scale_lin),1)+KSI0_INI;
% for t=2:length(time_scale_lin)
%     ksi_recalc_evol(t)=ksi_recalc_evol(t-1)+ksi_dot_recalc(t-1)*DELTA_TIME;
% end

% DELTA_TIME=(time_scale_lin(2)-time_scale_lin(1))*crash_duration;
% % ksi_dot_recalc=2*delta_evol.*vA_out_normal./(2*delta_theta0(1:length(vA_out)).*rx_evol_lin(1:length(vA_out)));
% ksi_dot_recalc=2*delta_evol.*vA_out_normal./Delta_SP(1:length(vA_out));
% ksi_dot_recalc(1)=0;
% ksi_dot_recalc(2)=0.5*ksi_dot_recalc(3);
% ksi_dot_recalc(TRANSITION_FRAME+1)=0.5*(ksi_dot_recalc(TRANSITION_FRAME));
% % ksi_dot_recalc=1.15*ksi_dot_recalc;
% 
% ksi_dot_recalc(TRANSITION_FRAME:length(time_scale_lin))=0;
% ksi_recalc_evol=zeros(length(time_scale_lin),1)+KSI0_INI;
% for t=2:length(time_scale_lin)
%     ksi_recalc_evol(t)=ksi_recalc_evol(t-1)+ksi_dot_recalc(t-1)*DELTA_TIME;
% end
% 
% %adjust the crash duration so that we get to rmix!
% crash_duration=0.5*(crash_duration+crash_duration*rmix/(max(ksi_recalc_evol(1:TRANSITION_FRAME))))


rx_recalc_evol=interp1(ksi0_evol_lin_ini(1:TRANSITION_FRAME_KSI-1),rx_evol_lin_ini(1:TRANSITION_FRAME_KSI-1),min(ksi_recalc_evol,0.99*rmix));
rx_recalc_evol(TRANSITION_FRAME:end)=max(rx_evol_lin_ini);
rx_recalc_evol(TRANSITION_FRAME-1)=0.5*(max(rx_evol_lin_ini)+rx_recalc_evol(TRANSITION_FRAME-2));
ksi_recalc_evol(TRANSITION_FRAME:end)=max(rx_evol_lin_ini);
ksi_recalc_evol=min(ksi_recalc_evol,rmix);

r_core_recalc=rx_recalc_evol-ksi_recalc_evol;
r_core_recalc=r_core_recalc';
r_core_recalc(TRANSITION_FRAME:end)=0.01*r_core_recalc(TRANSITION_FRAME-1);


% finally get back Delta from the most suited angle thetac
Delta_SP=0.5*(Delta_SP+2*(thetac_evol(1:length(vA_out))).*r_core_recalc(1:length(vA_out)));

%this first value makes little sense
Delta_SP(1)=Delta_SP(2);

crash_duration
rx_recalc_evol=(r_core_recalc+ksi_recalc_evol');
rx_dot_recalc_evol=gradient(rx_recalc_evol,time_scale_lin*crash_duration);
rx_dot_recalc_evol=max(rx_dot_recalc_evol,0);

% r_core=max(r_core_recalc,0);
% 
r_core=max(r_core_recalc,0);
rx_evol_lin=rx_recalc_evol;