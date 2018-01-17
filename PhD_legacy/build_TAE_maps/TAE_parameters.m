NB_THETA=NP;

NB_PHI=48+1;

rho_scale=radial_r_value_flux/max(radial_r_value_flux);
Btot_PR_map=sqrt(Bpol_PR_map.^2+Btor_PR_map.^2);

Bavg_radial=mean(Btot_PR_map(1:end-1,:),1);
Btot_radial_profile=mean(Btot_PR_map,1);
volume_flux_diff=volume_flux*0;
volume_flux_diff(2:Nradial)=volume_flux(2:end)-volume_flux(1:end-1);
% plasma_beta_recalc=mean(2*mu0*volume_flux_diff.*P_initial_profile)/mean(volume_flux_diff.*Btot_radial_profile.^2)

Ne_PR_map=zeros(NP,Nradial);
for p=1:NP
    Ne_PR_map(p,:)=Ne_profile;
end
vA_PR_map=Btot_PR_map./sqrt(mu0*mH*Ne_PR_map);
vA_PR_map(isnan(vA_PR_map))=0;
vA_profile=mean(vA_PR_map(1:end-1,:),1);


% build n=18 ; m=19 and m=20 TAE
Phi_PR_map=zeros(NP,Nradial);

nTAE=6
mtheta0=nTAE+3
mtheta1=mtheta0+1;
mtheta2=mtheta1+1;
% mtheta3=mtheta2+1;
mTAE=0.5*(mtheta1+mtheta2);

qmin=q_initial_profile(1)

% q0=max(mtheta0/nTAE,min(q_initial_profile));
q1=max(mtheta1/nTAE,qmin);
q2=mtheta2/nTAE;
% q3=mtheta3/nTAE;
qTAE=(mTAE)/nTAE
% r0=interp1(q_initial_profile,radial_r_value_flux,q0)
r1=interp1(q_initial_profile,radial_r_value_flux,(0.05*q1+0.95*qTAE))
r2=interp1(q_initial_profile,radial_r_value_flux,(0.05*q2+0.95*qTAE))
% r3=interp1(q_initial_profile,radial_r_value_flux,q3)
rTAE=interp1(q_initial_profile,radial_r_value_flux,qTAE);
psiTAE=interp1(q_initial_profile,psi_scale,qTAE);
psiposTAE=interp1(q_initial_profile,1:Nradial,qTAE);
psi1=interp1(q_initial_profile,psi_scale,q1);
psi2=interp1(q_initial_profile,psi_scale,q2);

Phi0=100;
MODE_WIDTH=2.2
TAE_WIDTH=3.2

% from benchmark
Ni0=2e19



TAE_angle=2*pi/nTAE
DPHI=TAE_angle/(NB_PHI-1)

MIN_PSI=2

pTAE_inf=max(round(interp1(radial_r_value_flux,1:Nradial,r1-7.2*(r2-r1))),MIN_PSI)
pTAE_sup=min(round(interp1(radial_r_value_flux,1:Nradial,r2+7.0*(r2-r1))),Nradial)


posTAE=round(interp1(q_initial_profile,1:Nradial,qTAE))
for p=1:NP
    q_PR_map(p,:)=q_initial_profile;
end
flc_inc=sqrt((X_PR_map.^2+Z_PR_map.^2)+(Rpos_PR_map.*q_PR_map).^2);

flc_s=flc_inc*0;


DTHETA=2*pi*(1/(NP-1));

for r=1:Nradial
%     flc_s(1,r)=0.5*flc_inc(1,r)*DPHI*q_initial_profile(r);
    flc_s(1,r)=0;
    for sval=2:NP
        flc_s(sval,r)=flc_s(sval-1,r)+0.5*(flc_inc(sval-1,r)+flc_inc(sval,r))*DTHETA;
    end
end

length_fl_profile=flc_s(end,:);
length_fl_TAE=interp1(q_initial_profile,flc_s(end,:),qTAE);

kTAE=(2*pi)/(2*length_fl_TAE);

omega_TAE=4.1e5 % from benchmark
vA_TAE=omega_TAE/kTAE

%hydrogen
pphi_TAE=mH*R0*vA_TAE-psiTAE;

% vA_TAE=interp1(q_initial_profile,vA_profile,qTAE)
% omega_TAE=vA_TAE*kTAE

xTAE=omega_TAE/vA_TAE
kpllm=(nTAE-mtheta1./q_initial_profile)/R0;
r1_der=interp1(q_initial_profile,radial_r_value_flux,(0.5*q1+0.5*qTAE))
r2_der=interp1(q_initial_profile,radial_r_value_flux,(0.5*q2+0.5*qTAE))

derivative_ratio=((radial_r_value_flux-r2_der)./(radial_r_value_flux-r1_der)).*(2.5*(radial_r_value_flux/R0)*xTAE^2-xTAE)./(xTAE^2-xTAE-kpllm.^2);
coupling_m_mp1=abs(interp1(radial_r_value_flux,derivative_ratio,rTAE))
% coupling_m_mp1=1.0
% coupling_m_mp1=1.3
coupling_m_mp1=1.1


USE_GRIDFIT=0
UNIFORM_TAE_PROPAGATION=1
% co propagatings waves : 1
% opposite waves : -1
OPP_W=-1;

% dominant mode number sup : 2
% dominant mode number : 1
% both harmonices : 0
DOM_M=0;
