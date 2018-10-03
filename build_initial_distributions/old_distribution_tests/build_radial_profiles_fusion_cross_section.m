clear all;
initialize_folder_names;

close all;

filename=strcat(DATA_FOLDER,'physics_constants.mat');
load(filename);
mr=mDT*mHe/(mDT+mHe);
mrDT=mD*mT/(mD+mT);
filename=strcat(DATA_FOLDER,'pressure_profile.mat');
load(filename);
filename=strcat(DATA_FOLDER,'motions_map_dimensions.mat');
load(filename);


Te_psi=Te_profile/eV;
Ti_psi=Te_profile/eV;
Ni0_psi=Ne_profile;
Ne0_psi=Ne_profile;

disp('-------------------------------------------------')
vthe=sqrt(2*eV*Te_psi/me);
vthi=sqrt(2*eV*Te_psi/mDT);


lambda_D=sqrt(epsilon0*Te_psi*eV/(Ne0_psi*eV^2));
Lambda_ei=9*(4*pi/3*Ne0_psi*lambda_D^3);

log_lambda=log(Lambda_ei);

tau_e=(3*(2*pi)^1.5)*(epsilon0^2*sqrt(me)*(eV*Te_psi).^1.5)./(Ne0_psi.*log_lambda*eV^4);
tau_eHe=(3*(2*pi)^1.5)*(epsilon0^2*sqrt(me)*(eV*Te_psi).^1.5)./(Ne0_psi.*log_lambda*4*eV^4);

tau_s=((Ni0_psi*mDT)./(Ne0_psi*me)).*tau_eHe;

% Number of radial points
NumberX=4000;

% Main important values for description of the distribution
w0=3.5*(10^6)*eV;
v0=sqrt(2*w0/mHe);

% approximative parameters
vc=sqrt(44*2*Te_psi*eV/mHe);


DeltaX=1/(NumberX-1);
%Ti0 in keV
averaged_cross_section=1.1*(1e-24)*(Ti_psi/1000).^2;
Salpha=0.25*Ni0_psi.*Ni0_psi.*averaged_cross_section;

barn=1e-28;
Tmax=296;
sigma_m=5.03;
vth_r=sqrt(2*eV*1000*Tmax/mrDT);

sigma_v_profile=(4*sigma_m/sqrt(3))*vth_r*((1000*Tmax./Ti_psi).^(2/3)).*exp(-3*(1000*Tmax./Ti_psi).^(1/3)+2);
sigma_v_max=16/(9*sqrt(3))*sigma_m*vth_r;


C1=1.17302e-9;
C2=1.51361e-2;
C3=7.51886e-2;
C4=4.60643e-3;
C5=1.35e-2;
C6=-1.0675e-4;
C7=1.366e-5;
Bg=34.3827;
mBH=1124656;

Ti_keV=Ti_psi/1000; 
Ti_range=Ti_keV;

theta_bh=Ti_range./(1-(Ti_range.*(C2+Ti_range.*(C4+Ti_range*C6)))./(1+Ti_range.*(C3+Ti_range.*(C5+Ti_range*C7))));
ksi_bh=((Bg^2)./(4*theta_bh)).^(1/3);
sigma_v_bosch_hale=C1*theta_bh.*sqrt(ksi_bh./(mBH*Ti_range.^3)).*exp(-3*ksi_bh);
sigma_v_bosch_hale=sigma_v_bosch_hale*1e-6; % SI units

% 
% for (x=1:2*NumberX-1)
%     speed_range(x)=(x-1)*DeltaX;
%     vb(x)=vc*speed_range(x);
% end
% 
% for (x=0.25*NumberX:2*NumberX-1)
%     nubi(x)=(1/(4*pi))*(4*Ni0_psi*log_lambda*eV^4)./((epsilon0^2)*mr*mHe)*(1/(vb(x)^3+1.3*vthi.^3));
%     nube(x)=(1/(4*pi))*(4*Ne0_psi*log_lambda*eV^4)./((epsilon0^2)*me*mHe)*(1/(vb(x)^3+1.3*vthe.^3));
% end
% 
% % need to saturate collision rate
% % because the low speed part of the curve is not
% % what we have modelled and anyway we only use these curves
% % to find the critical speed transition from electrons to ions collision
% % dominated slowing down
% 
% for (x=1:1000)
%     nubi(x)=nubi(1000);
%     nube(x)=nube(1000);
% end
% 
% [vc_eps vc_pos]=min(abs(nubi-nube));
% vc=speed_range(vc_pos)*vc;
% wcrit=(0.5*mHe*vc^2)/eV;
% 
% 
% % redo the collision rate curves with corrected speed range
% for (x=1:2*NumberX-1)
%     vb(x)=v0*speed_range(x);
%     vbc(x)=vc*speed_range(x);
% end
% for (x=1000:2*NumberX-1)
%     nubi(x)=(1/(4*pi))*(4*Ni0_psi*log_lambda*eV^4)/((epsilon0^2)*mr*mHe)*(1/(vbc(x)^3+1.3*vthi^3));
%     nube(x)=(1/(4*pi))*(4*Ne0_psi*log_lambda*eV^4)/((epsilon0^2)*me*mHe)*(1/(vbc(x)^3+1.3*vthe^3));
% end
% for (x=1:1000)
%     nubi(x)=nubi(1000);
%     nube(x)=nube(1000);
% end
% 
% 
