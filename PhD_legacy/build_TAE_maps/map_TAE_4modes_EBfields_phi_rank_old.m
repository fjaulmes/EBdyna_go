
phi_nTAE=2*pi*(phi_rank-1)/(NB_PHI-1)

for r=1:Nradial
    vAvalues=vA_PR_map(:,r);
    calculate_temporal_mode_structures;
    r_value=radial_r_value_flux(r);
    sb_width1=exp(-((r_value-r0)^2)/(2*MODE_WIDTH*(r3-r0)));
    sb_width2=exp(-((r_value-r3)^2)/(2*MODE_WIDTH*(r3-r0)));
    m_width1=exp(-((r_value-r1)^2)/(2*MODE_WIDTH*(r3-r0)));
    m_width2=exp(-((r_value-r2)^2)/(2*MODE_WIDTH*(r3-r0)));
    global_width=exp(-((r_value-rTAE)^2)/(2*TAE_WIDTH*(r2-r1)));
    
    Phi_PR_map(:,r)=Phi0*global_width*(sb_width1*ksi_mm1_t+m_width1*ksi_m_t+m_width2*ksi_mp1_t+sb_width2*ksi_mp2_t);
    iPhi_PR_map(:,r)=Phi0*global_width*(sb_width1*iksi_mm1_t+m_width1*iksi_m_t+m_width2*iksi_mp1_t+sb_width2*iksi_mp2_t);
    A_PR_map(:,r)=(kTAE./omega_values').*Phi_PR_map(:,r);
    iA_PR_map(:,r)=(kTAE./omega_values').*iPhi_PR_map(:,r);
    omega_PR_map(:,r)=omega_values';
end
% imagesc(Phi_PR_map')

psi_star_PR_map(:,1:pTAE_sup)=-A_PR_map(:,1:pTAE_sup).*Rpos_PR_map(:,1:pTAE_sup).*Btor_PR_map(:,1:pTAE_sup)./Btot_PR_map(:,1:pTAE_sup);


Phi_data=reshape((Phi_PR_map(:,1:pTAE_sup)),NP*pTAE_sup,1);
Phi_XZ_map=griddata(finesse_data_X,finesse_data_Z,Phi_data,XX_small,ZZ_small,'cubic');
Phi_XZ_map(isnan(Phi_XZ_map)) = 0; 
Phi_XZ_map=Phi_XZ_map';

omega_data=reshape((omega_PR_map(:,1:pTAE_sup)),NP*pTAE_sup,1);
omega_XZ_map=griddata(finesse_data_X,finesse_data_Z,omega_data,XX_small,ZZ_small,'cubic');
omega_XZ_map(isnan(omega_XZ_map)) = 0; 
omega_XZ_map=omega_XZ_map';
omega_XZ_map(omega_XZ_map==0) = omega_TAE; 


A_XZ_map=(kTAE./omega_XZ_map).*Phi_XZ_map;

A_data=reshape((iA_PR_map(:,1:pTAE_sup)),NP*pTAE_sup,1);
iA_XZ_map=griddata(finesse_data_X,finesse_data_Z,A_data,XX_small,ZZ_small,'cubic');
iA_XZ_map(isnan(iA_XZ_map)) = 0; 
iA_XZ_map=iA_XZ_map';

AX_XZ_map=A_XZ_map.*bX_XZ_small_map;
AZ_XZ_map=A_XZ_map.*bZ_XZ_small_map;

gAZ_R=zeros(sizeX,sizeZ);
gAX_Z=zeros(sizeX,sizeZ);

for (x=3:sizeX-2)
    gAZ_R(x,3:sizeZ-2)=(1/12)*(-AZ_XZ_map(x+2,(3:sizeZ-2))+AZ_XZ_map(x-2,(3:sizeZ-2)))+(2/3)*(AZ_XZ_map(x+1,(3:sizeZ-2))-AZ_XZ_map(x-1,(3:sizeZ-2)));
    gAX_Z(x,3:sizeZ-2)=(1/12)*(-AX_XZ_map(x,(3:sizeZ-2)+2)+AX_XZ_map(x,(3:sizeZ-2)-2))+(2/3)*(AX_XZ_map(x,(3:sizeZ-2)+1)-AX_XZ_map(x,(3:sizeZ-2)-1));   
end
gAZ_R=gAZ_R/DX;
gAX_Z=gAX_Z/DX;
Bstarphi_XZ_map=(gAZ_R-gAX_Z);


psi_star_XZ_map=-A_XZ_map.*Rpos_XZsmall_map.*bphi_XZsmall_map;
% psi_star_XZsmall_map=psi_star_XZ_map(Xinf:Xsup,Zinf:Zsup);
% psi_data=reshape(psi_star_XZ_map(:,:)',sizeX*sizeZ,1);
% psi_star_PR_map=griddata(X_scale_data,Z_scale_data,psi_data,RR,Z_PR_map(:,1:pTAE_sup),'cubic');
% psi_star_PR_map(isnan(psi_star_PR_map)) = 0; 


psi2D=psi_star_XZ_map;

% psi2D=psi_star_XZ_map;

gpsi_R=zeros(sizeX,sizeZ);
gpsi_Z=zeros(sizeX,sizeZ);

for (x=3:sizeX-2)
    gpsi_R(x,3:sizeZ-2)=(1/12)*(-psi2D(x+2,(3:sizeZ-2))+psi2D(x-2,(3:sizeZ-2)))+(2/3)*(psi2D(x+1,(3:sizeZ-2))-psi2D(x-1,(3:sizeZ-2)));
    gpsi_Z(x,3:sizeZ-2)=(1/12)*(-psi2D(x,(3:sizeZ-2)+2)+psi2D(x,(3:sizeZ-2)-2))+(2/3)*(psi2D(x,(3:sizeZ-2)+1)-psi2D(x,(3:sizeZ-2)-1));
end
gpsi_R=gpsi_R/DX;
gpsi_Z=gpsi_Z/DX;

% BstarX_XZ_map=-gpsi_Z./Rpos_XZsmall_map;
% BstarZ_XZ_map=gpsi_R./Rpos_XZsmall_map;

BstarX_XZ_map=-gpsi_Z./Rpos_XZsmall_map-(iA_XZ_map.*bZ_XZ_small_map)./Rpos_XZsmall_map;
BstarZ_XZ_map=gpsi_R./Rpos_XZsmall_map+(iA_XZ_map.*bX_XZ_small_map)./Rpos_XZsmall_map;

% BstarX_XZ_map=-gpsi_Z./Rpos_XZsmall_map-iA_XZ_map./Rpos_XZsmall_map;
% BstarZ_XZ_map=gpsi_R./Rpos_XZsmall_map+iA_XZ_map./Rpos_XZsmall_map;


BpolX_XZ_map=BstarX_XZ_map+BpolX_initial_XZsmall_map;
BpolZ_XZ_map=BstarZ_XZ_map+BpolZ_initial_XZsmall_map;

Bstar_XZ_map_frame_rank=sqrt(BstarX_XZ_map.^2+BstarZ_XZ_map.^2);
Btot_XZ_map_frame_rank=sqrt(BpolX_XZ_map.^2+BpolZ_XZ_map.^2+Bphi_XZsmall_map.^2);


% B_data=reshape(Bstar_XZ_map_frame_rank(:,:)',sizeX*sizeZ,1);
% Btot_PR_map=griddata(X_scale_data,Z_scale_data,B_data,RR,Z_PR_map(:,1:pTAE_sup),'cubic');

B_data=reshape(BstarX_XZ_map(:,:)',sizeX*sizeZ,1);
bX_PR_map=griddata(X_scale_data,Z_scale_data,B_data,RR,Z_PR_map(:,1:pTAE_sup),'cubic');
% bX_PR_map=bX_PR_map./Btot_PR_map;
bX_PR_map(isnan(bX_PR_map)) = 0; 

B_data=reshape(BstarZ_XZ_map(:,:)',sizeX*sizeZ,1);
bZ_PR_map=griddata(X_scale_data,Z_scale_data,B_data,RR,Z_PR_map(:,1:pTAE_sup),'cubic');
% bZ_PR_map=bZ_PR_map./Btot_PR_map;
bZ_PR_map(isnan(bZ_PR_map)) = 0; 

B_data=reshape(Bstarphi_XZ_map(:,:)',sizeX*sizeZ,1);
bphi_PR_map=griddata(X_scale_data,Z_scale_data,B_data,RR,Z_PR_map(:,1:pTAE_sup),'cubic');
bphi_PR_map(isnan(bphi_PR_map)) = 0; 


gPhi_R=zeros(sizeX,sizeZ);
gPhi_Z=zeros(sizeX,sizeZ);

for (x=3:sizeX-2)
    gPhi_R(x,3:sizeZ-2)=(1/12)*(-Phi_XZ_map(x+2,(3:sizeZ-2))+Phi_XZ_map(x-2,(3:sizeZ-2)))+(2/3)*(Phi_XZ_map(x+1,(3:sizeZ-2))-Phi_XZ_map(x-1,(3:sizeZ-2)));
    gPhi_Z(x,3:sizeZ-2)=(1/12)*(-Phi_XZ_map(x,(3:sizeZ-2)+2)+Phi_XZ_map(x,(3:sizeZ-2)-2))+(2/3)*(Phi_XZ_map(x,(3:sizeZ-2)+1)-Phi_XZ_map(x,(3:sizeZ-2)-1));
end
gPhi_R=gPhi_R/DX;
gPhi_Z=gPhi_Z/DX;


%simplified Efield expression for TAE (weak A||)

EX_XZ_map=-gPhi_R;
EZ_XZ_map=-gPhi_Z;

Ephi_XZ_map=-(EX_XZ_map.*BpolX_XZ_map+EZ_XZ_map.*BpolZ_XZ_map)./(Bphi_XZsmall_map+Bstarphi_XZ_map);
Etot_XZ_map_frame_rank=sqrt(EX_XZ_map.^2+EZ_XZ_map.^2+Ephi_XZ_map.^2);

%estimate the displacement
ksi_dot_XZ_map_frame_rank=Etot_XZ_map_frame_rank./Btot_XZ_map_frame_rank;
exb_X=EZ_XZ_map.*Bphi_XZsmall_map-Ephi_XZ_map.*BpolZ_XZ_map;
exb_Z=-EX_XZ_map.*Bphi_XZsmall_map+Ephi_XZ_map.*BpolX_XZ_map;
exb_phi=EX_XZ_map.*BpolZ_XZ_map-EZ_XZ_map.*BpolX_XZ_map;
norm_exb=sqrt(exb_X.^2+exb_Z.^2+exb_phi.^2);
exb_X=exb_X./norm_exb;
exb_Z=exb_Z./norm_exb;
exb_phi=exb_phi./norm_exb;
exb_X(isnan(exb_X))=0;
exb_Z(isnan(exb_Z))=0;
exb_phi(isnan(exb_phi))=0;

ksi_X=exb_X.*ksi_dot_XZ_map_frame_rank;
ksi_Z=exb_Z.*ksi_dot_XZ_map_frame_rank;
ksi_phi=exb_phi.*ksi_dot_XZ_map_frame_rank;

gRksi_R=zeros(sizeX,sizeZ);
gksi_Z=zeros(sizeX,sizeZ);

Rksi_X=Rpos_XZsmall_map.*ksi_X;

for (x=3:sizeX-2)
    gRksi_R(x,3:sizeZ-2)=(1/12)*(-Rksi_X(x+2,(3:sizeZ-2))+Rksi_X(x-2,(3:sizeZ-2)))+(2/3)*(Rksi_X(x+1,(3:sizeZ-2))-Rksi_X(x-1,(3:sizeZ-2)));
    gksi_Z(x,3:sizeZ-2)=(1/12)*(-ksi_Z(x,(3:sizeZ-2)+2)+ksi_Z(x,(3:sizeZ-2)-2))+(2/3)*(ksi_Z(x,(3:sizeZ-2)+1)-ksi_Z(x,(3:sizeZ-2)-1));
end
gRksi_R=gRksi_R/DX;
gksi_Z=gksi_Z/DX;

% approximate divergence, without toroidal derivative
div_ksi=gRksi_R./Rpos_XZsmall_map+gksi_Z;
P1=-(ksi_X.*gP0_X+ksi_Z.*gP0_Z)-gamma_thermo*P0_XZ_map.*div_ksi;
P1=P1./omega_XZ_map;
div_ksi=div_ksi./omega_XZ_map;

E_data=reshape(EX_XZ_map(:,:)',sizeX*sizeZ,1);
EX_PR_map=griddata(X_scale_data,Z_scale_data,E_data,RR,Z_PR_map(:,1:pTAE_sup),'cubic');
EX_PR_map(isnan(EX_PR_map))=0;

E_data=reshape(EZ_XZ_map(:,:)',sizeX*sizeZ,1);
EZ_PR_map=griddata(X_scale_data,Z_scale_data,E_data,RR,Z_PR_map(:,1:pTAE_sup),'cubic');
EZ_PR_map(isnan(EZ_PR_map))=0;

