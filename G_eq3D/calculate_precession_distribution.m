




q0=q_initial_profile(1)
Omega_ci=(ZHe*eV)*Bavg/mHe;
alphas_lambda=Bavg*alphas_mm./alphas_Ekin;
alphas_v=sqrt(2*(eV/mHe)*alphas_Ekin);

rho_pll0=sqrt(1-alphas_lambda).*alphas_v/Omega_ci;
v_pll0=sqrt(1-alphas_lambda).*alphas_v;
sign_v=zeros(Nalphas_simulated,1);
% Omega_phi=zeros(Nalphas_simulated,1);
omega_ik_avg=zeros(Nalphas_simulated,1);
cos_theta_avg=zeros(Nalphas_simulated,1);
vpll_avg=zeros(Nalphas_simulated,1);
omega_ae_avg=zeros(Nalphas_simulated,1);
omega_ntm32_avg=zeros(Nalphas_simulated,1);
omega_ntm43_avg=zeros(Nalphas_simulated,1);
omega_psi_avg=zeros(Nalphas_simulated,1);
r_vpll_minus_avg=zeros(Nalphas_simulated,1);
r_vpll_plus_avg=zeros(Nalphas_simulated,1);
omega_r_avg=zeros(Nalphas_simulated,1);
delta_r_avg=zeros(Nalphas_simulated,1);
omega_precess_avg=zeros(Nalphas_simulated,1);
omega_bounce=zeros(Nalphas_simulated,1);
max_vpll_avg=zeros(Nalphas_simulated,1);
sigma_avg=zeros(Nalphas_simulated,1);
theta_bounce=zeros(Nalphas_simulated,1);
phidot_avg=zeros(Nalphas_simulated,1);

for number=1:Nalphas_simulated
    polynomial_fit=polyfit(time_scale_corr,phi_output_pop(1:end_ts,number),1);
    phidot_avg(number)=polynomial_fit(1);
    polynomial_fit=polyfit(time_scale_corr,6.5*theta_output_pop(:,number)-6*phi_output_pop(:,number),1);
	omega_ae_avg(number)=polynomial_fit(1);
    polynomial_fit=polyfit(time_scale_corr,3*theta_output_pop(:,number)-2*phi_output_pop(:,number),1);
	omega_ntm32_avg(number)=polynomial_fit(1);
    polynomial_fit=polyfit(time_scale_corr,4*theta_output_pop(:,number)-3*phi_output_pop(:,number),1);
	omega_ntm43_avg(number)=polynomial_fit(1);
	polynomial_fit=polyfit(time_scale_corr,theta_output_pop(:,number)-phi_output_pop(:,number),1);
	omega_ik_avg(number)=polynomial_fit(1);
end


for number=1:Nalphas_simulated
    [max_vpll max_index]=max(abs(vpll_output_pop(:,number)));
    max_vpll_avg(number)=vpll_output_pop(max_index,number);
%     sign_v(number)=sign(vpll_output_pop(max_index,number));
end
sign_v=sign(max_vpll_avg);

rho_pll0=sign_v.*rho_pll0;
v_pll0=sign_v.*v_pll0;

disp('time_scale(end_ts)=');
disp(time_scale_corr(end_ts));


vDZ_avg=mean(vDZ_output_corr(:,:),1);
vDphi_avg=mean(vDphi_output_corr(:,:),1);
vpll_avg=mean(vpll_output_pop(:,:),1);
cos_theta_avg=mean(cos(theta_output_pop(:,:)),1);
r_avg=mean(r_output_pop(:,:),1)';
q_avg=mean(q_output_pop(:,:),1)';


alphas_Eperp=B_avg.*alphas_mm;
alphas_vperp=sqrt(2*alphas_Eperp*eV/mHe);
alphas_rhoL=(mHe/eV)*(1/ZHe)*alphas_vperp./B_avg;

r_avg=r_avg';
Ravg=R0+mean(X_output_pop(:,:),1)';
time=time_scale(end)

for number=1:Nalphas_simulated
    fft_theta_pop=abs(fft(theta_output_pop(:,number)));
    [fft_r_max fft_theta_freq]=max(fft_theta_pop(2:end));
    omega_bounce(number)=2*pi*((fft_theta_freq))/time;
    [min_value theta_pos]=min(abs(theta_output_pop(2:end,number)-pi));
    theta_bounce(number)=theta_output_pop(theta_pos,number);
end

clear theta_output_pop q_output_pop



Nfft=size(r_output_pop,1);
for number=1:Nalphas_simulated
    fft_r_pop=abs(fft(r_output_pop(:,number)));
    [fft_r_max fft_r_freq]=max(fft_r_pop(2:end));
    omega_r_avg(number)=2*pi*((fft_r_freq))/time;
    r_max=max(r_output_pop(:,number));
    r_min=min(r_output_pop(:,number));
    if ~isempty(find(abs(POTATOES-number)==0))
        delta_r_avg(number)=r_max+r_min;
    else
        delta_r_avg(number)=r_max-r_min;
    end
end


clear r_output_pop

bX_output_pop=interp2(scale_X,scale_Z,bX_XZ_map',X_output_pop,Z_output_pop,'*linear');
bZ_output_pop=interp2(scale_X,scale_Z,bZ_XZ_map',X_output_pop,Z_output_pop,'*linear');
bphi_output_pop=sqrt(1-(bX_output_pop.^2+bZ_output_pop.^2));


for number=1:Nalphas_simulated
%     number=ALL_TRAPPED(n);
    omega_precess_avg(number)=mean(GX_theta_output_pop(:,number).*vDX_output_corr(:,number)+...
	GZ_theta_output_pop(:,number).*vDZ_output_corr(:,number)-...
       vDphi_output_corr(:,number)./(R0+X_output_pop(:,number)));
    omega_psi_avg(number)=mean(GX_theta_output_pop(:,number).*bX_output_pop(:,number).*vpll_output_pop(:,number)+...
	   GZ_theta_output_pop(:,number).*bZ_output_pop(:,number).*vpll_output_pop(:,number)-...
       bphi_output_pop(:,number).*vpll_output_pop(:,number)./(R0+X_output_pop(:,number)));
end

clear bphi_output_pop vDphi_output_corr vDZ_output_corr


rX_output_pop=interp2(scale_X,scale_Z,drX_map',X_output_pop,Z_output_pop,'*linear');
rZ_output_pop=interp2(scale_X,scale_Z,drZ_map',X_output_pop,Z_output_pop,'*linear');

% omega_psi_avg(ALL_PASSING)=mean((q_output_pop(:,ALL_PASSING)-1).*vpll_output_pop(:,ALL_PASSING)./(R0+X_output_pop(:,ALL_PASSING)),1);

for number=1:Nalphas_simulated
    r_vpll_plus_avg(number)=mean(squeeze((sign(vpll_output_pop(:,number))==1).*sqrt((rX_output_pop(:,number).*bX_output_pop(:,number).*vpll_output_pop(:,number)).^2+...
       (rZ_output_pop(:,number).*bZ_output_pop(:,number).*vpll_output_pop(:,number)).^2)));
	r_vpll_minus_avg(number)=mean(squeeze((sign(vpll_output_pop(:,number))==-1).*sqrt((rX_output_pop(:,number).*bX_output_pop(:,number).*vpll_output_pop(:,number)).^2+...
       (rZ_output_pop(:,number).*bZ_output_pop(:,number).*vpll_output_pop(:,number)).^2)));
end

clear bX_output_pop bZ_output_pop rX_output_pop rZ_output_pop




% time=time_scale(end)





tau_precess=2*pi./abs(omega_precess_avg);

NO_REDIST=find((tau_precess<0.7e-4).*abs((omega_precess_avg./omega_psi_avg)>0.5));


omega_crash=0.5*pi/(0.5*tau_cr)

ratio_longitudinal_limit=1/3

