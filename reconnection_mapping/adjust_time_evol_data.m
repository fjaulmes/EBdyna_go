% now we also need to get back to a smaller radial size mesh
psi_star_2D_evol_result=zeros(Ntime+1,size_r,Nomega);

for t=1:Ntime+1
    for om=1:Nomega
        psi_star_2D_evol_result(t,:,om)=interp1(scale_r,squeeze(psi_star_2D_evol(t,:,om)),radial_r_value(1:size_r),'cubic');
    end
end

psi_star_2D_evol=psi_star_2D_evol_result;

time_lin_transition=time_evol(length(volume1_evol))



% maximal_time_step=min(0.5*(time_evol(3:end)-time_evol(1:end-2)));

%%
% now getting the data on a regular time scale

TIME_STEP=0.01
NB_TIME_STEP=1/TIME_STEP+1
time_scale_lin=(0:NB_TIME_STEP-1)*TIME_STEP;

time_step_transition=interp1(time_scale_lin,1:length(time_scale_lin),time_lin_transition)
time_step_transition=floor(time_step_transition);

volume1_evol_lin=interp1(time_evol(1:length(volume1_evol)),volume1_evol,time_scale_lin(1:time_step_transition),'cubic');
volume2_evol_lin=interp1(time_evol(1:length(volume2_evol)),volume2_evol,time_scale_lin(1:time_step_transition),'cubic');
volume3_evol_lin=interp1(time_evol(1:length(volume3_evol)),volume3_evol,time_scale_lin(1:time_step_transition),'cubic');

Bstar1_evol_lin=interp1(time_evol(1:length(Bstar1_evol)),Bstar1_evol,time_scale_lin(1:time_step_transition),'cubic');
Bstar2_evol_lin=interp1(time_evol(1:length(Bstar2_evol)),Bstar2_evol,time_scale_lin(1:time_step_transition),'cubic');
Bstar3_evol_lin=interp1(time_evol(1:length(Bstar3_evol)),Bstar3_evol,time_scale_lin(1:time_step_transition),'cubic');

P1_evol_lin=interp1(time_evol(1:length(P1_evol)),P1_evol,time_scale_lin(1:time_step_transition),'cubic');
P2_evol_lin=interp1(time_evol(1:length(P2_evol)),P2_evol,time_scale_lin(1:time_step_transition),'cubic');
P3_evol_lin=interp1(time_evol(1:length(P3_evol)),P3_evol,time_scale_lin(1:time_step_transition),'cubic');

PI1_evol_lin=interp1(time_evol(1:length(PI1_evol)),PI1_evol,time_scale_lin(1:time_step_transition),'cubic');
PI2_evol_lin=interp1(time_evol(1:length(PI2_evol)),PI2_evol,time_scale_lin(1:time_step_transition),'cubic');
PI3_evol_lin=interp1(time_evol(1:length(PI3_evol)),PI3_evol,time_scale_lin(1:time_step_transition),'cubic');

PE1_evol_lin=interp1(time_evol(1:length(PE1_evol)),PE1_evol,time_scale_lin(1:time_step_transition),'cubic');
PE2_evol_lin=interp1(time_evol(1:length(PE2_evol)),PE2_evol,time_scale_lin(1:time_step_transition),'cubic');
PE3_evol_lin=interp1(time_evol(1:length(PE3_evol)),PE3_evol,time_scale_lin(1:time_step_transition),'cubic');



ksi0_evol_lin=interp1(time_evol,ksi0_evol,time_scale_lin,'cubic');
rx_evol_lin=interp1(time_evol,rx_evol,time_scale_lin,'cubic');

psi_star_2D_evol_lin=interp1(time_evol,psi_star_2D_evol,time_scale_lin,'cubic');



FILENAME=strcat(DATA_FOLDER,'psi_star_evol.mat')
% save (FILENAME,'psi_star_2D_evol','ksi0_evol','rx_evol','time_evol','size_r','Nomega','rx_evol_interp','psi_star_2D_evol_interp','ksi0_evol_interp','time_scale_precise','time_scale_lin','psi_star_2D_evol_lin','NB_TIME_STEP','TIME_STEP');
save (FILENAME,'psi_star_2D_evol','ksi0_evol','rx_evol','volume1_evol_lin','volume2_evol_lin','volume3_evol_lin',...
    'Bstar1_evol_lin','Bstar2_evol_lin','Bstar3_evol_lin','P1_evol_lin','P2_evol_lin','P3_evol_lin',...
	'PI1_evol_lin','PI2_evol_lin','PI3_evol_lin','PE1_evol_lin','PE2_evol_lin','PE3_evol_lin',...
	'time_evol','size_r','Nomega','time_scale_lin','ksi0_evol_lin','rx_evol_lin','psi_star_2D_evol_lin','NB_TIME_STEP','TIME_STEP');


FILENAME=strcat(DATA_FOLDER,'reconnection_kol_coefs_evol.mat')
save (FILENAME,'der_cont13_1','der_cont13_3','der_cont23_3','der_cont23_2','time_scale_der','xi_scale_der','a1_coef_evol','a2_coef_evol');



% probably no need for thousand frames.....

% TIME_STEP=0.001
% NB_TIME_STEP=1/TIME_STEP+1
% time_scale_precise=(0:NB_TIME_STEP-1)*TIME_STEP;
% 
% ksi0_evol_interp=interp1(time_evol,ksi0_evol,time_scale_precise,'cubic');
% psi_star_2D_evol_interp=interp1(time_evol,psi_star_2D_evol,time_scale_precise,'cubic');
% rx_evol_interp=interp1(time_evol,rx_evol,time_scale_precise,'cubic');
% 
% FILENAME=strcat(DATA_FOLDER,'psi_star_evol_interp.mat')
% save (FILENAME,'psi_star_2D_evol','ksi0_evol','rx_evol','time_evol','size_r','Nomega','rx_evol_interp','psi_star_2D_evol_interp','ksi0_evol_interp','time_scale_precise','time_scale_lin','psi_star_2D_evol_lin','NB_TIME_STEP','TIME_STEP');
