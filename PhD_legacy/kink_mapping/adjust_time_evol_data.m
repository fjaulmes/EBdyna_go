% now we also need to get back to a smaller radial size mesh
Ntime=size(psi_star_2D_evol,1)
psi_star_2D_evol_result=zeros(Ntime,size_r,Nomega);

for t=1:Ntime
    for om=1:Nomega
        psi_star_2D_evol_result(t,:,om)=interp1(scale_r,squeeze(psi_star_2D_evol(t,:,om)),radial_r_value(1:size_r),'cubic');
    end
end

psi_star_2D_evol=psi_star_2D_evol_result;




% maximal_time_step=min(0.5*(time_evol(3:end)-time_evol(1:end-2)));

% now getting the data on a regular time scale

TIME_STEP=0.01
NB_TIME_STEP=1/TIME_STEP+1
time_scale_lin=(0:NB_TIME_STEP-1)*TIME_STEP;

ksi0_evol_lin=interp1(time_evol,ksi0_evol,time_scale_lin,'cubic');
rx_evol_lin=interp1(time_evol,rx_evol,time_scale_lin,'cubic');

psi_star_2D_evol_lin=interp1(time_evol,psi_star_2D_evol,time_scale_lin,'cubic');



FILENAME=strcat(DATA_FOLDER,'psi_star_evol.mat')
% save (FILENAME,'psi_star_2D_evol','ksi0_evol','rx_evol','time_evol','size_r','Nomega','rx_evol_interp','psi_star_2D_evol_interp','ksi0_evol_interp','time_scale_precise','time_scale_lin','psi_star_2D_evol_lin','NB_TIME_STEP','TIME_STEP');
save (FILENAME,'psi_star_2D_evol','ksi0_evol','rx_evol','time_evol','size_r','Nomega','time_scale_lin','ksi0_evol_lin','rx_evol_lin','psi_star_2D_evol_lin','NB_TIME_STEP','TIME_STEP');

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
