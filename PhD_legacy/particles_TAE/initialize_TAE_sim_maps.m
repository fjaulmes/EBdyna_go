filename=strcat(DATA_FOLDER,'q_profile.mat');
load(filename,'q_initial_profile');
%filename=strcat(DATA_FOLDER,'flux_geometry.mat');
%load(filename,'X_PR_map');
%load(filename,'Z_PR_map');

%for p=1:NB_THETA
%    q_PR_map(p,:)=q_initial_profile;
%end
%flc_inc=sqrt((X_PR_map.^2+Z_PR_map.^2)+((R0+X_PR_map).*q_PR_map).^2);
%
%flc_s=flc_inc*0;
%
%DTHETA=2*pi*(1/(NB_THETA-1));
%
%for r=1:Nradial
%    flc_s(1,r)=0;
%    for sval=2:NB_THETA
%        flc_s(sval,r)=flc_s(sval-1,r)+0.5*(flc_inc(sval-1,r)+flc_inc(sval,r))*DTHETA;
%    end
%end


	
gamma_TAE_min=-2*omega_TAE;
gamma_TAE_max=2*omega_TAE;

Btot_PR_map_ini=sqrt(Bpol_PR_map.^2+Btor_PR_map.^2);

dr_profile=radial_r_value_flux*0;
dr_profile(2:end)=radial_r_value_flux(2:end)-radial_r_value_flux(1:end-1);
dr_profile(2:end-1)=0.5*(-radial_r_value_flux(1:end-2)+radial_r_value_flux(3:end));
dr_PR_map=Rpos_PR_map*0;
r_PR_map=Rpos_PR_map*0;
for p=1:NB_THETA
    dr_PR_map(p,:)=dr_profile;
    r_PR_map(p,:)=radial_r_value_flux;
end

Btot_map_phi_ini=zeros(NB_PHI,NB_THETA,size_r);
BpolX_map_phi_ini=zeros(NB_PHI,NB_THETA,size_r);
BpolZ_map_phi_ini=zeros(NB_PHI,NB_THETA,size_r);
Btor_map_phi_ini=zeros(NB_PHI,NB_THETA,size_r);

for n=1:NB_PHI
    Btot_map_phi_ini(n,:,:)=Btot_PR_map_ini(:,pTAE_inf:pTAE_sup);
    BpolX_map_phi_ini(n,:,:)=BX_PR_map(:,pTAE_inf:pTAE_sup);
    BpolZ_map_phi_ini(n,:,:)=BZ_PR_map(:,pTAE_inf:pTAE_sup);
    Btor_map_phi_ini(n,:,:)=Btor_PR_map(:,pTAE_inf:pTAE_sup);
end

volume_circ_approx_PR_map=Rpos_PR_map.*r_PR_map.*dr_PR_map*DTHETA*DPHI;



% corrected theta maps for complete interpolation
QNB_THETA=round(0.25*NB_THETA);
HQNB_THETA=round(0.5*QNB_THETA);


% correction needed for old TAE_data ....
% need to fix this in map calaculations

kpllm1_profile(1:pTAE_inf)=kpllm1_profile(pTAE_inf);
kpllm1_profile(pTAE_sup:Nradial)=kpllm1_profile(pTAE_sup);

kpllm2_profile(1:pTAE_inf)=kpllm2_profile(pTAE_inf);
kpllm2_profile(pTAE_sup:Nradial)=kpllm2_profile(pTAE_sup);

m_adjust_profile1(1:pTAE_inf)=m_adjust_profile1(pTAE_inf);
m_adjust_profile1(pTAE_sup:Nradial)=m_adjust_profile1(pTAE_sup);

m_adjust_profile2(1:pTAE_inf)=m_adjust_profile2(pTAE_inf);
m_adjust_profile2(pTAE_sup:Nradial)=m_adjust_profile2(pTAE_sup);