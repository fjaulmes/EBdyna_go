


% pitch_pop=atan(r_output_pop./(q_output_pop.*(R0+X_output_pop)));
bX_pop=interp2(scale_X,scale_Z,bZ_XZ_map',X_output_pop,Z_output_pop,'*linear');
bZ_pop=interp2(scale_X,scale_Z,bZ_XZ_map',X_output_pop,Z_output_pop,'*linear');
bphi_pop=sqrt(1-(bX_pop.^2+bZ_pop.^2));

% calculate_gradZ_theta_maps;
gZ_theta_output_pop=interp2(scale_X,scale_Z,gradZ_theta_map',X_output_pop,Z_output_pop,'*nearest');
gZ_theta_avg=mean(gZ_theta_output_pop(:,:),1)';


omega_precess_avg=zeros(Nalphas_simulated,1);
omega_psi_avg=zeros(Nalphas_simulated,1);


disp('time_scale(end_ts)=');
disp(time_scale_corr(end_ts));
delta_time=time_scale_corr(1);


alphas_Eperp=B_avg.*alphas_mm;
alphas_vperp=sqrt(2*alphas_Eperp*eV/mHe);
alphas_rhoL=(mHe/eV)*(1/ZHe)*alphas_vperp./B_avg;


for number=1:Nalphas_simulated
%     number=ALL_TRAPPED(n);
    omega_precess_avg(number)=mean(gZ_theta_output_pop(:,number).*vDZ_output_corr(:,number)-...
    vDphi_output_corr(:,number)./(R0+X_output_pop(:,number)));
    omega_psi_avg(number)=mean(gZ_theta_output_pop(:,number).*(vpll_output_pop(:,number).*bZ_pop(:,number))-...
    (vpll_output_pop(:,number).*bphi_pop(:,number))./(R0+X_output_pop(:,number)));
end

% omega_psi_avg=mean((q_output_pop-1).*vpll_output_pop./(R0+X_output_pop),1);
% omega_psi_avg(COUNTER_PASSING)=mean((q_output_pop(:,COUNTER_PASSING)-1),1)'.*omega_bounce(COUNTER_PASSING);
% omega_psi_avg(CO_PASSING)=mean(-(q_output_pop(:,CO_PASSING)-1),1)'.*omega_bounce(CO_PASSING);
% omega_psi_avg=omega_psi_avg';

% omega_psi_avg=mean((q_output_pop-1).*vpll_output_pop./(R0),1);
% omega_psi_avg=omega_psi_avg';


tau_precess=2*pi./abs(omega_precess_avg);

NO_REDIST=find((tau_precess<0.7e-4).*abs((omega_precess_avg./omega_psi_avg)>0.5));


precession_compared_with_longitudinal_motion_lambda;






