




q0=q_initial_profile(1)
Omega_ci=(ZHe*eV)*Bavg/mHe;
alphas_lambda=Bavg*alphas_mm./alphas_Ekin;
alphas_v=sqrt(2*(eV/mHe)*alphas_Ekin);

rho_pll0=sqrt(1-alphas_lambda).*alphas_v/Omega_ci;
v_pll0=sqrt(1-alphas_lambda).*alphas_v;
sign_v=zeros(Nalphas_simulated,1);
Omega_phi=zeros(Nalphas_simulated,1);
Omega_omega=zeros(Nalphas_simulated,1);
cos_theta_avg=zeros(Nalphas_simulated,1);
vpll_avg=zeros(Nalphas_simulated,1);
omega_psi_avg=zeros(Nalphas_simulated,1);
omega_r_avg=zeros(Nalphas_simulated,1);
delta_r_avg=zeros(Nalphas_simulated,1);
omega_precess_avg=zeros(Nalphas_simulated,1);
omega_bounce=zeros(Nalphas_simulated,1);
max_vpll_avg=zeros(Nalphas_simulated,1);
sigma_avg=zeros(Nalphas_simulated,1);
theta_bounce=zeros(Nalphas_simulated,1);



% load('initialG_alphas_precession.mat')

% for number=1:size(theta_output,2)
%     theta_timeline=theta_output(1:end_ts,number);
%     straighten_theta_timeline;
%     theta_timeline=theta_timeline_corr;
%     straighten_theta_timeline;
%     theta_output_corr(1:end_ts,number)=theta_timeline_corr;
% end

for number=1:Nalphas_simulated
    [max_vpll max_index]=max(abs(vpll_output_pop(:,number)));
    max_vpll_avg(number)=vpll_output_pop(max_index,number);
%     sign_v(number)=sign(vpll_output_pop(max_index,number));
end
sign_v=sign(max_vpll_avg);

% for number=1:Nalphas_simulated
%     [max_value pos_theta]=max(abs(r_output_pop(:,number)));
%     sign_v(number)=sign(vpll_output_pop(pos_theta,number));
% end
rho_pll0=sign_v.*rho_pll0;
v_pll0=sign_v.*v_pll0;

disp('time_scale(end_ts)=');
disp(time_scale_corr(end_ts));

% for number=1:Nalphas_simulated
%     polynomial_fit=polyfit(time_scale_corr,phipos_output_corr(1:end_ts,number),1);
%     Omega_phi(number)=polynomial_fit(1);
% end

% for number=1:Nalphas_simulated
%     polynomial_fit=polyfit(time_scale_corr,omega_output_corr(1:end_ts,number),1);
%     Omega_omega(number)=polynomial_fit(1);
% end

vDZ_avg=mean(vDZ_output_corr(:,:),1);
vDphi_avg=mean(vDphi_output_corr(:,:),1);
vpll_avg=mean(vpll_output_pop(:,:),1);
cos_theta_avg=mean(cos(theta_output_pop(:,:)),1);
r_avg=mean(r_output_pop(:,:),1)';
q_avg=mean(q_output_pop(:,:),1)';
B_avg=mean(B_output_pop(:,:),1)';

alphas_Eperp=B_avg.*alphas_mm;
alphas_vperp=sqrt(2*alphas_Eperp*eV/mHe);
alphas_rhoL=(mHe/eV)*(1/ZHe)*alphas_vperp./B_avg;

% alphas_kappa(ALL_TRAPPED)=sin();

% alphas_kappa=sqrt((alphas_Ekin*Raxis+Bavg*alphas_mm.*(r_avg-Raxis))./(2*alphas_mm.*r_avg*Bavg));
% omega_precess_avg=mean((-vDphi_output_corr/R0),1);
r_avg=r_avg';
% omega_precess_avg=(1./r_avg).*(-mean(vDZ_output_corr.*cos(theta_output_pop(:,:)),1))-vDphi_avg/R0;
Ravg=R0+mean(X_output_pop(:,:),1)';

% for n=1:length(ALL_PASSING)
%     number=ALL_PASSING(n);
%     
%     omega_precess_avg(number)=-mean(vDphi_output_corr(:,number)./(R0+X_output_pop(:,number)),1);
% end
% 
% for n=1:length(ALL_TRAPPED)
%     number=ALL_TRAPPED(n);
%     omega_precess_avg(number)=mean(-vDZ_output_corr(1:last_pos_pop(number),number).*cos(theta_output_pop(1:last_pos_pop(number),number))-...
%         vDphi_output_corr(1:last_pos_pop(number),number)./(R0+X_output_pop(1:last_pos_pop(number),number)));
% end

for number=1:Nalphas_simulated
%     number=ALL_TRAPPED(n);
    omega_precess_avg(number)=mean(GZ_theta_output_pop(:,number).*vDZ_output_corr(:,number)-...
       vDphi_output_corr(:,number)./(R0+X_output_pop(:,number)));
    omega_psi_avg(number)=mean(GZ_theta_output_pop(:,number).*bZ_output_pop(:,number).*vpll_output_pop(:,number)-...
       bphi_output_pop(:,number).*vpll_output_pop(:,number)./(R0+X_output_pop(:,number)));
end




% L=length(time_scale_corr);
% NFFT = 2^nextpow2(L)


time=time_scale(end)
Nfft=size(r_output_pop,1);
for number=1:Nalphas_simulated
    fft_r_pop=abs(fft(r_output_pop(:,number)));
    [fft_r_max fft_r_freq]=max(fft_r_pop(2:end));
    omega_r_avg(number)=2*pi*((fft_r_freq))/time;
    delta_r_avg(number)=4*fft_r_max/Nfft;
end

% time=time_scale(end)

for number=1:Nalphas_simulated
    fft_theta_pop=abs(fft(theta_output_pop(:,number)));
    [fft_r_max fft_theta_freq]=max(fft_theta_pop(2:end));
    omega_bounce(number)=2*pi*((fft_theta_freq))/time;
    [min_value theta_pos]=min(abs(theta_output_pop(2:end,number)-pi));
    theta_bounce(number)=theta_output_pop(theta_pos,number);
end

% omega_psi_avg(ALL_PASSING)=omega_bounce(ALL_PASSING).*mean((q_output_pop(ALL_PASSING)-1),1);
% omega_psi_avg(ALL_TRAPPED)=omega_precess_avg(ALL_TRAPPED).*mean((q_output_pop(ALL_PASSING)-1),1);

% omega_psi_avg=omega_bounce'.*abs(q_avg-1);

% omega_bounce_avg=((alphas_mm.*Bavg)./(q_avg.*Ravg)).*sqrt(r_avg./(2*Ravg));


% for number=1:Nalphas_simulated
%     omega_psi_avg(number)=-mean((1./(R0+X_output_pop(1:end_ts,number))).*abs(q_output_pop(1:end_ts,number)-1).*vpll_output_pop(1:end_ts,number));
%     omega_precess_avg(number)=mean((1./(1e-6+r_output_pop(1:end_ts,number))).*(-vDZ_output_corr(1:end_ts,number).*cos(theta_output_pop(1:end_ts,number)))-vDphi_output_corr(1:end_ts,number)/R0);
% %     omega_precess_avg(number)=mean((-vDphi_output_corr(1:end_ts,number)/R0));
% end

% alphas_vD_Z=interp2(scale_X,scale_Z,vD_Z_XZ_map',alphas_pos_x,alphas_pos_z,'*linear');
% alphas_vD_Z=((2*alphas_Ekin-Bavg*alphas_mm)/(ZHe*Bavg)).*alphas_vD_Z;
% Omega_precess_th=(alphas_vD_Z.*cos_theta_avg)./(radial_pos);

% TRAPPED=real(alphas_kappa)<=1;
% [K E]=ellipke(real(alphas_kappa(TRAPPED)));
% G=2*E./K-1;
% r_avg=r_avg';
% Omega_precess_trapped_th=(vDZ_avg(ALL_TRAPPED))./(r_avg(ALL_TRAPPED));
% 
% PASSING=find((real(alphas_kappa)>1).*(alphas_lambda<1));
% Omega_phi_passing_th=(v_pll0(PASSING)/R0).*(1-rho_pll0(PASSING)*(q0/(elongation*R0)).*(1-2*alphas_lambda(PASSING)+0.75*alphas_lambda(PASSING).^2)./((1-alphas_lambda(PASSING)).^2));

% Omega_precess_passing_th=Omega_phi(PASSING)-omega_psi_avg(PASSING);
% Omega_precess_passing_th=-Omega_precess_passing_th;
% Omega_precess_th=Omega_omega-omega_psi_avg';
% omega_psi=0.05*vpll_avg/R0;
% omega_pre=Omega_omega-omega_psi;

% OTHER=find((real(alphas_kappa)<1).*(alphas_lambda<1));

tau_precess=2*pi./abs(omega_precess_avg);

NO_REDIST=find((tau_precess<0.7e-4).*abs((omega_precess_avg./omega_psi_avg)>0.5));




figure(1);
xlabel('\lambda_0');
subplot(2,1,1)
ylabel('\omega_{pre}');
hold on
plot(alphas_lambda(ALL_TRAPPED),abs(omega_precess_avg(ALL_TRAPPED)),'r.');
plot(alphas_lambda(ALL_PASSING),abs(omega_precess_avg(ALL_PASSING)),'b.');
subplot(2,1,2)
hold on
plot(alphas_lambda(ALL_PASSING),abs(omega_precess_avg(ALL_PASSING))./abs(omega_psi_avg(ALL_PASSING)),'b.');
plot(alphas_lambda(ALL_TRAPPED),abs(omega_precess_avg(ALL_TRAPPED))./abs(omega_psi_avg(ALL_TRAPPED)),'r.');
ylim([0 10])






% plot(alphas_kappa(TRAPPED),log10(1e-6+abs(omega_precess_avg(TRAPPED)./omega_psi_avg(TRAPPED))),'r.');
% plot(alphas_kappa(PASSING),log10(1e-6+abs(omega_precess_avg(PASSING)./omega_psi_avg(PASSING))),'b.');

figure(4);
xlabel('Ekin_0');
subplot(2,1,1)
ylabel('\omega_{pre}');
hold on
plot(alphas_Ekin(ALL_PASSING),abs(omega_precess_avg(ALL_PASSING)),'b.');
plot(alphas_Ekin(ALL_TRAPPED),abs(omega_precess_avg(ALL_TRAPPED)),'r.');
subplot(2,1,2)
hold on
plot(alphas_Ekin(ALL_PASSING),abs(omega_psi_avg(ALL_PASSING)),'b.');
plot(alphas_Ekin(ALL_TRAPPED),abs(omega_psi_avg(ALL_TRAPPED)),'r.');


PRECESS_LAMBDA_BIN_SIZE=0.1;
precession_lambda_bins=(0:PRECESS_LAMBDA_BIN_SIZE:12*PRECESS_LAMBDA_BIN_SIZE);
lambda_values=(0.5*PRECESS_LAMBDA_BIN_SIZE:PRECESS_LAMBDA_BIN_SIZE:11.5*PRECESS_LAMBDA_BIN_SIZE);

precession_Ekin_bins=[0 400 1000 2000 3600]*1e3;
Ekin_values=[200 700 1500 2800]*1e3;

lambda_precess_value=zeros(length(precession_lambda_bins)-1,length(precession_Ekin_bins)-1);

for Ebin=1:(length(precession_Ekin_bins)-1)
    PRECESS_EKIN_POP=((alphas_Ekin>precession_Ekin_bins(Ebin)).*(alphas_Ekin<=precession_Ekin_bins(Ebin+1)));
    for bin=1:(length(precession_lambda_bins)-1)
        PRECESS_LAMBDA_POP=find((alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*PRECESS_EKIN_POP);
        lambda_precess_value(bin,Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP)));
    end
end

omega_crash=0.5*pi/(0.5*tau_cr)

figure(5);
set(gca,'FontSize',16);

hold on
grid on
hold on
grid on
plot(lambda_values,lambda_precess_value(:,1),'g','LineWidth',1);
plot(lambda_values,lambda_precess_value(:,2),'b','LineWidth',2);
plot(lambda_values,lambda_precess_value(:,3),'k','LineWidth',3);
plot(lambda_values,lambda_precess_value(:,4),'r','LineWidth',4);
plot(lambda_values,lambda_values.*0+omega_crash,'y--','LineWidth',2);
xlim([lambda_values(1) lambda_values(end)])

Ekin_values=Ekin_values*1e-3;
legend(strcat('Ekin=',num2str(Ekin_values(1)),'keV'),strcat('Ekin=',num2str(Ekin_values(2)),'keV'),strcat('Ekin=',num2str(Ekin_values(3)),'keV'),strcat('Ekin=',num2str(Ekin_values(4)),'keV'));

xlabel('\lambda_0')
ylabel('\omega_{pre}')


lambda_precess_value=zeros(length(precession_lambda_bins)-1,length(precession_Ekin_bins)-1);

for Ebin=1:(length(precession_Ekin_bins)-1)
    PRECESS_EKIN_POP=((alphas_Ekin>precession_Ekin_bins(Ebin)).*(alphas_Ekin<=precession_Ekin_bins(Ebin+1)));
    for bin=1:(length(precession_lambda_bins)-1)
        PRECESS_LAMBDA_POP=find((alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*PRECESS_EKIN_POP);
        lambda_precess_value(bin,Ebin)=mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP)));
    end
end

omega_crash=0.5*pi/(0.5*tau_cr)

figure(5);
set(gca,'FontSize',16);

hold on
grid on
hold on
grid on
plot(lambda_values,lambda_precess_value(:,1),'g-.','LineWidth',1);
plot(lambda_values,lambda_precess_value(:,2),'b-.','LineWidth',2);
plot(lambda_values,lambda_precess_value(:,3),'k-.','LineWidth',3);
plot(lambda_values,lambda_precess_value(:,4),'r-.','LineWidth',4);
plot(lambda_values,lambda_values.*0+omega_crash,'y--','LineWidth',2);
xlim([lambda_values(1) lambda_values(end)])

Ekin_values=Ekin_values*1e-3;
legend(strcat('Ekin=',num2str(Ekin_values(1)),'keV'),strcat('Ekin=',num2str(Ekin_values(2)),'keV'),strcat('Ekin=',num2str(Ekin_values(3)),'keV'),strcat('Ekin=',num2str(Ekin_values(4)),'keV'));

xlabel('\lambda_0')
ylabel('\omega_{pre}')






PRECESS_KAPPA_BIN_SIZE=0.05;
precession_kappa_bins=(8*PRECESS_KAPPA_BIN_SIZE:PRECESS_KAPPA_BIN_SIZE:82*PRECESS_KAPPA_BIN_SIZE);
kappa_values=(8.5*PRECESS_KAPPA_BIN_SIZE:PRECESS_KAPPA_BIN_SIZE:81.5*PRECESS_KAPPA_BIN_SIZE);

PRECESS_EKIN_BIN_SIZE=400;
precession_Ekin_bins=[0 400 1000 2000 3600]*1e3;
Ekin_values=[200 700 1500 2800]*1e3;

kappa_precess_value=zeros(length(precession_kappa_bins)-1,length(precession_Ekin_bins)-1);

for Ebin=1:(length(precession_Ekin_bins)-1)
    PRECESS_EKIN_POP=((alphas_Ekin>precession_Ekin_bins(Ebin)).*(alphas_Ekin<=precession_Ekin_bins(Ebin+1)));
    for bin=1:(length(precession_kappa_bins)-1)
        PRECESS_KAPPA_POP=find((alphas_kappa>precession_kappa_bins(bin)).*(alphas_kappa<=precession_kappa_bins(bin+1)).*PRECESS_EKIN_POP);
        kappa_precess_value(bin,Ebin)=mean(abs(omega_precess_avg(PRECESS_KAPPA_POP)));
    end
end


figure(6);
hold on
grid on
plot(kappa_values,kappa_precess_value(:,1),'g','LineWidth',2);
plot(kappa_values,kappa_precess_value(:,2),'b','LineWidth',2);
plot(kappa_values,kappa_precess_value(:,3),'k','LineWidth',2);
plot(kappa_values,kappa_precess_value(:,4),'r','LineWidth',2);
plot(kappa_values,kappa_values.*0+omega_crash,'y--','LineWidth',2);
xlim([kappa_values(1) kappa_values(end)])







PRECESS_KAPPA_BIN_SIZE=0.05;
precession_kappa_bins=(8*PRECESS_KAPPA_BIN_SIZE:PRECESS_KAPPA_BIN_SIZE:82*PRECESS_KAPPA_BIN_SIZE);
kappa_values=(8.5*PRECESS_KAPPA_BIN_SIZE:PRECESS_KAPPA_BIN_SIZE:81.5*PRECESS_KAPPA_BIN_SIZE);

PRECESS_EKIN_BIN_SIZE=400;
precession_Ekin_bins=[0 400 1000 2000 3600]*1e3;
Ekin_values=[200 700 1500 2800]*1e3;

kappa_precess_value=zeros(length(precession_kappa_bins)-1,length(precession_Ekin_bins)-1);

for Ebin=1:(length(precession_Ekin_bins)-1)
    PRECESS_EKIN_POP=((alphas_Ekin>precession_Ekin_bins(Ebin)).*(alphas_Ekin<=precession_Ekin_bins(Ebin+1)));
    for bin=1:(length(precession_kappa_bins)-1)
        PRECESS_KAPPA_POP=find((alphas_kappa>precession_kappa_bins(bin)).*(alphas_kappa<=precession_kappa_bins(bin+1)).*PRECESS_EKIN_POP);
        kappa_precess_value(bin,Ebin)=mean(abs(omega_psi_avg(PRECESS_KAPPA_POP)));
    end
end


% figure(7);
hold on
grid on
plot(kappa_values,kappa_precess_value(:,1),'g-.','LineWidth',2);
plot(kappa_values,kappa_precess_value(:,2),'b-.','LineWidth',2);
plot(kappa_values,kappa_precess_value(:,3),'k-.','LineWidth',2);
plot(kappa_values,kappa_precess_value(:,4),'r-.','LineWidth',2);
% plot(kappa_values,kappa_values.*0+omega_crash,'y--','LineWidth',2);
xlim([kappa_values(1) kappa_values(end)])






