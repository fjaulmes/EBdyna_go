
phi_rank
m_adjust_profile=zeros(Nradial,1);
m_guess=q_initial_profile(pTAE_inf-1)*nTAE;
m_adjust_profile1(pTAE_inf-1)=0.9*(m_guess-mtheta1);
m_adjust_profile2(pTAE_inf-1)=0.9*(m_guess-mtheta2);
for r=pTAE_inf:pTAE_sup
    r_value=radial_r_value_flux(r);
    Rvalues=Rpos_PR_map(:,r);
%     vAvalues=vA_PR_map(:,r);
    vAvalues=vA_TAE;
    calculate_temporal_mode_structures;
    
    %     sb_width1=exp(-((r_value-r0)^2)/(2*MODE_WIDTH*(r3-r0)));
    %     sb_width2=exp(-((r_value-r3)^2)/(2*MODE_WIDTH*(r3-r0)));
    m_width1=exp(-((r_value-r1)^2)/(2*MODE_WIDTH*(r3-r0)));
    m_width2=exp(-((r_value-r2)^2)/(2*MODE_WIDTH*(r3-r0)));
    global_width=exp(-((r_value-rTAE)^2)/(2*TAE_WIDTH*(r2-r1)));
    
    Phi_PR_map(:,r)=Phi0*global_width*(m_width1*ksi_m_t+m_width2*ksi_mp1_t);
    iPhi_PR_map(:,r)=Phi0*global_width*(m_width1*iksi_m_t+m_width2*iksi_mp1_t);
%     dPhi_PR_map(:,r)=-(omega_TAE/nTAE)*Phi0*global_width*(iksi_m_t);
    if UNIFORM_TAE_PROPAGATION==0
        A_PR_map(:,r)=(kTAE./omega_values').*Phi_PR_map(:,r);
        iA_PR_map(:,r)=(kTAE./omega_values').*iPhi_PR_map(:,r);
        dA_PR_map(:,r)=-(kTAE/nTAE).*iPhi_PR_map(:,r);
        omega_PR_map(:,r)=omega_values';
    else
        A_PR_map(:,r)=(kTAE/omega_TAE)*Phi_PR_map(:,r);
        iA_PR_map(:,r)=(kTAE/omega_TAE)*iPhi_PR_map(:,r);
        dA_PR_map(:,r)=-(kTAE/nTAE)*iPhi_PR_map(:,r);
    end
    
end

psi_star_PR_map(:,1:pTAE_sup)=-A_PR_map(:,1:pTAE_sup).*Rpos_PR_map(:,1:pTAE_sup).*Btor_PR_map(:,1:pTAE_sup)./Btot_PR_map(:,1:pTAE_sup);

Phi_data=reshape((Phi_PR_map(:,1:pTAE_sup)),NP*pTAE_sup,1);
Phi_XZ_map=griddata(finesse_data_X,finesse_data_Z,Phi_data,XX_small,ZZ_small,'cubic');
Phi_XZ_map(isnan(Phi_XZ_map)) = 0; 
Phi_XZ_map=Phi_XZ_map';

if UNIFORM_TAE_PROPAGATION==0
    omega_data=reshape((omega_PR_map(:,1:pTAE_sup)),NP*pTAE_sup,1);
    omega_XZ_map=griddata(finesse_data_X,finesse_data_Z,omega_data,XX_small,ZZ_small,'cubic');
    omega_XZ_map(isnan(omega_XZ_map)) = 0;
    omega_XZ_map=omega_XZ_map';
    omega_XZ_map(omega_XZ_map==0) = omega_TAE;
end

if UNIFORM_TAE_PROPAGATION==0
    A_XZ_map=(kTAE./omega_XZ_map).*Phi_XZ_map;
else
    A_XZ_map=(kTAE/omega_TAE)*Phi_XZ_map;
end


A_data=reshape((iA_PR_map(:,1:pTAE_sup)),NP*pTAE_sup,1);
iA_XZ_map=griddata(finesse_data_X,finesse_data_Z,A_data,XX_small,ZZ_small,'cubic');
iA_XZ_map(isnan(iA_XZ_map)) = 0; 
iA_XZ_map=iA_XZ_map';

iPhi_XZ_map=(omega_TAE/kTAE)*iA_XZ_map;



% J0 current vector

RBphi=Rpos_XZsmall_map.*Bphi_XZsmall_map;

gBphi_Z=zeros(sizeX,sizeZ);
gRBphi_X=zeros(sizeX,sizeZ);

for (x=3:sizeX-2)
    gBphi_Z(x,3:sizeZ-2)=(1/12)*(-Bphi_XZsmall_map(x,(3:sizeZ-2)+2)+Bphi_XZsmall_map(x,(3:sizeZ-2)-2))+(2/3)*(Bphi_XZsmall_map(x,(3:sizeZ-2)+1)-Bphi_XZsmall_map(x,(3:sizeZ-2)-1));
    gRBphi_X(x,3:sizeZ-2)=(1/12)*(-RBphi(x+2,(3:sizeZ-2))+RBphi(x-2,(3:sizeZ-2)))+(2/3)*(RBphi(x+1,(3:sizeZ-2))-RBphi(x-1,(3:sizeZ-2)));
end
gBphi_Z=gBphi_Z/DX;
gRBphi_X=gRBphi_X/DX;

J0_X_XZ_map=(gBphi_Z)/mu0;
J0_Z_XZ_map=(-gRBphi_X./Rpos_XZsmall_map)/mu0;

gBZ_R=zeros(sizeX,sizeZ);
gBX_Z=zeros(sizeX,sizeZ);

for (x=3:sizeX-2)
    gBZ_R(x,3:sizeZ-2)=(1/12)*(-BpolZ_initial_XZsmall_map(x+2,(3:sizeZ-2))+BpolZ_initial_XZsmall_map(x-2,(3:sizeZ-2)))+(2/3)*(BpolZ_initial_XZsmall_map(x+1,(3:sizeZ-2))-BpolZ_initial_XZsmall_map(x-1,(3:sizeZ-2)));
    gBX_Z(x,3:sizeZ-2)=(1/12)*(-BpolX_initial_XZsmall_map(x,(3:sizeZ-2)+2)+BpolX_initial_XZsmall_map(x,(3:sizeZ-2)-2))+(2/3)*(BpolX_initial_XZsmall_map(x,(3:sizeZ-2)+1)-BpolX_initial_XZsmall_map(x,(3:sizeZ-2)-1));   
end
gBZ_R=gBZ_R/DX;
gBX_Z=gBX_Z/DX;
J0_phi_XZ_map=(gBZ_R-gBX_Z)/mu0;




if UNIFORM_TAE_PROPAGATION==0
    A_data=reshape((dA_PR_map(:,1:pTAE_sup)),NP*pTAE_sup,1);
    dA_XZ_map=griddata(finesse_data_X,finesse_data_Z,A_data,XX_small,ZZ_small,'cubic');
    dA_XZ_map(isnan(dA_XZ_map)) = 0;
    dA_XZ_map=dA_XZ_map';
else
    dA_XZ_map=-(omega_TAE/nTAE)*iA_XZ_map;
end

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
ipsi_star_XZ_map=-iA_XZ_map.*Rpos_XZsmall_map.*bphi_XZsmall_map;
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

% would like to have the picture with potentials
% for algorithm validation

% BstarX_XZ_map=-gpsi_Z./Rpos_XZsmall_map;
% BstarZ_XZ_map=gpsi_R./Rpos_XZsmall_map;

BstarX_XZ_map=-gpsi_Z./Rpos_XZsmall_map-(iA_XZ_map.*bZ_XZ_small_map)./Rpos_XZsmall_map;
BstarZ_XZ_map=gpsi_R./Rpos_XZsmall_map+(iA_XZ_map.*bX_XZ_small_map)./Rpos_XZsmall_map;


BpolX_XZ_map=BstarX_XZ_map+BpolX_initial_XZsmall_map;
BpolZ_XZ_map=BstarZ_XZ_map+BpolZ_initial_XZsmall_map;

Bstar_XZ_map_frame_rank=sqrt(BstarX_XZ_map.^2+BstarZ_XZ_map.^2+Bstarphi_XZ_map.^2);
Btot_XZ_map_rank=sqrt(BpolX_XZ_map.^2+BpolZ_XZ_map.^2+(Bphi_XZsmall_map+Bstarphi_XZ_map).^2);


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


%time derivative of Efield 

igPhi_R=zeros(sizeX,sizeZ);
igPhi_Z=zeros(sizeX,sizeZ);

for (x=3:sizeX-2)
    igPhi_R(x,3:sizeZ-2)=(1/12)*(-iPhi_XZ_map(x+2,(3:sizeZ-2))+iPhi_XZ_map(x-2,(3:sizeZ-2)))+(2/3)*(iPhi_XZ_map(x+1,(3:sizeZ-2))-iPhi_XZ_map(x-1,(3:sizeZ-2)));
    igPhi_Z(x,3:sizeZ-2)=(1/12)*(-iPhi_XZ_map(x,(3:sizeZ-2)+2)+iPhi_XZ_map(x,(3:sizeZ-2)-2))+(2/3)*(iPhi_XZ_map(x,(3:sizeZ-2)+1)-iPhi_XZ_map(x,(3:sizeZ-2)-1));
end
igPhi_R=igPhi_R/DX;
igPhi_Z=igPhi_Z/DX;

iEX_XZ_map=-igPhi_R+(omega_TAE^2)*AZ_XZ_map.*bX_XZ_small_map;
iEZ_XZ_map=-igPhi_Z+(omega_TAE^2)*AZ_XZ_map.*bZ_XZ_small_map;
iEphi_XZ_map=-(iEX_XZ_map.*BpolX_XZ_map+iEZ_XZ_map.*BpolZ_XZ_map)./(Bphi_XZsmall_map+Bstarphi_XZ_map);

% now find vector ksi

iexb_X=iEZ_XZ_map.*Bphi_XZsmall_map-iEphi_XZ_map.*BpolZ_XZ_map;
iexb_Z=-iEX_XZ_map.*Bphi_XZsmall_map+iEphi_XZ_map.*BpolX_XZ_map;
iexb_phi=iEX_XZ_map.*BpolZ_XZ_map-iEZ_XZ_map.*BpolX_XZ_map;
inorm_exb=sqrt(iexb_X.^2+iexb_Z.^2+iexb_phi.^2);
% inorm_exb=inorm_exb+1;
iexb_X=iexb_X./inorm_exb;
iexb_Z=iexb_Z./inorm_exb;
iexb_phi=iexb_phi./inorm_exb;
iexb_X(isnan(iexb_X))=0;
iexb_Z(isnan(iexb_Z))=0;
iexb_phi(isnan(iexb_phi))=0;


% Efield 

EX_XZ_map=-gPhi_R-dA_XZ_map.*bX_XZ_small_map;
EZ_XZ_map=-gPhi_Z-dA_XZ_map.*bZ_XZ_small_map;

Ephi_XZ_map=-(EX_XZ_map.*BpolX_XZ_map+EZ_XZ_map.*BpolZ_XZ_map)./(Bphi_XZsmall_map+Bstarphi_XZ_map);
Etot_XZ_map_rank=sqrt(EX_XZ_map.^2+EZ_XZ_map.^2+Ephi_XZ_map.^2);

%estimate the displacement
ksi_dot_XZ_map_frame_rank=Etot_XZ_map_rank./Btot_XZ_map_rank;
% if UNIFORM_TAE_PROPAGATION==0
%     ksi_dot_XZ_map_frame_rank=omega_XZ_map.*A_XZ_map./Btot_XZ_map_rank;
% else
%     ksi_dot_XZ_map_frame_rank=omega_TAE*A_XZ_map./Btot_XZ_map_rank;
% end
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


% ksi is unknown

% ksi_X=exb_X;
% ksi_Z=exb_Z;
% ksi_phi=exb_phi;
% % if UNIFORM_TAE_PROPAGATION==0
% %     ksi_X=ksi_X./omega_XZ_map;
% %     ksi_Z=ksi_Z./omega_XZ_map;
% %     ksi_phi=ksi_phi./omega_XZ_map;
% % else
% %     ksi_X=ksi_X/omega_TAE;
% %     ksi_Z=ksi_Z/omega_TAE;
% %     ksi_phi=ksi_phi/omega_TAE;
% % end
% 
% ksi_tot=sqrt(ksi_X.^2+ksi_Z.^2+ksi_phi.^2);
% 
% ksi_B0_X=ksi_Z.*Bphi_XZsmall_map-ksi_phi.*BpolZ_XZ_map;
% ksi_B0_Z=-ksi_X.*Bphi_XZsmall_map+ksi_phi.*BpolX_XZ_map;
% ksi_B0_phi=ksi_X.*BpolZ_XZ_map-ksi_Z.*BpolX_XZ_map;
% 
% psi2D=-Rpos_XZsmall_map.*ksi_B0_phi;
% 
% % psi2D=psi_star_XZ_map;
% 
% gpsi_R=zeros(sizeX,sizeZ);
% gpsi_Z=zeros(sizeX,sizeZ);
% 
% for (x=3:sizeX-2)
%     gpsi_R(x,3:sizeZ-2)=(1/12)*(-psi2D(x+2,(3:sizeZ-2))+psi2D(x-2,(3:sizeZ-2)))+(2/3)*(psi2D(x+1,(3:sizeZ-2))-psi2D(x-1,(3:sizeZ-2)));
%     gpsi_Z(x,3:sizeZ-2)=(1/12)*(-psi2D(x,(3:sizeZ-2)+2)+psi2D(x,(3:sizeZ-2)-2))+(2/3)*(psi2D(x,(3:sizeZ-2)+1)-psi2D(x,(3:sizeZ-2)-1));
% end
% gpsi_R=gpsi_R/DX;
% gpsi_Z=gpsi_Z/DX;
% 
% 
% % rot_ksi_B0_X=rot_ksi_B0_X/omega_TAE;
% % rot_ksi_B0_Z=rot_ksi_B0_Z/omega_TAE;
% % rot_ksi_B0_phi=rot_ksi_B0_phi/omega_TAE;
% 
% gRksi_R=zeros(sizeX,sizeZ);
% gksi_Z=zeros(sizeX,sizeZ);
% 
% Rksi_X=Rpos_XZsmall_map.*ksi_X;
% 
% for (x=3:sizeX-2)
%     gRksi_R(x,3:sizeZ-2)=(1/12)*(-Rksi_X(x+2,(3:sizeZ-2))+Rksi_X(x-2,(3:sizeZ-2)))+(2/3)*(Rksi_X(x+1,(3:sizeZ-2))-Rksi_X(x-1,(3:sizeZ-2)));
%     gksi_Z(x,3:sizeZ-2)=(1/12)*(-ksi_Z(x,(3:sizeZ-2)+2)+ksi_Z(x,(3:sizeZ-2)-2))+(2/3)*(ksi_Z(x,(3:sizeZ-2)+1)-ksi_Z(x,(3:sizeZ-2)-1));
% end
% gRksi_R=gRksi_R/DX;
% gksi_Z=gksi_Z/DX;
% 
% % approximate divergence, without toroidal derivative
% div_ksi=gRksi_R./Rpos_XZsmall_map+gksi_Z;
% P1=-(ksi_X.*gP0_X+ksi_Z.*gP0_Z)-gamma_thermo*P0_XZ_map.*div_ksi;
% % if UNIFORM_TAE_PROPAGATION==0
% %     P1=P1./omega_XZ_map;
% %     div_ksi=div_ksi./omega_XZ_map;
% % else
% %     P1=P1/omega_TAE;
% %     div_ksi=div_ksi/omega_TAE;
% % end

E_data=reshape(EX_XZ_map(:,:)',sizeX*sizeZ,1);
EX_PR_map=griddata(X_scale_data,Z_scale_data,E_data,RR,Z_PR_map(:,1:pTAE_sup),'cubic');
EX_PR_map(isnan(EX_PR_map))=0;

E_data=reshape(EZ_XZ_map(:,:)',sizeX*sizeZ,1);
EZ_PR_map=griddata(X_scale_data,Z_scale_data,E_data,RR,Z_PR_map(:,1:pTAE_sup),'cubic');
EZ_PR_map(isnan(EZ_PR_map))=0;


%%
% check div(B)=0

gRB_R=zeros(sizeX,sizeZ);
gB_Z=zeros(sizeX,sizeZ);

RB_X=BstarX_XZ_map.*Rpos_XZsmall_map;
B_Z=BstarZ_XZ_map;

for (x=3:sizeX-2)
    gRB_R(x,3:sizeZ-2)=(1/12)*(-RB_X(x+2,(3:sizeZ-2))+RB_X(x-2,(3:sizeZ-2)))+(2/3)*(RB_X(x+1,(3:sizeZ-2))-RB_X(x-1,(3:sizeZ-2)));
    gB_Z(x,3:sizeZ-2)=(1/12)*(-B_Z(x,(3:sizeZ-2)+2)+B_Z(x,(3:sizeZ-2)-2))+(2/3)*(B_Z(x,(3:sizeZ-2)+1)-B_Z(x,(3:sizeZ-2)-1));
end
gRB_R=gRB_R/DX;
gB_Z=gB_Z/DX;

% approximate divergence, without toroidal derivative
% div_B=gRB_R./Rpos_XZsmall_map+gB_Z+iBstarphi_XZ_map./Rpos_XZsmall_map;


ksi_XZ_map=Etot_XZ_map_rank./Btot_XZ_map_rank/omega_TAE;

ksi_X=ksi_XZ_map.*iexb_X;
ksi_Z=ksi_XZ_map.*iexb_Z;
ksi_phi=ksi_XZ_map.*iexb_phi;

J0_B1_X=J0_Z_XZ_map.*Bstarphi_XZ_map-J0_phi_XZ_map.*BstarZ_XZ_map;
J0_B1_Z=-J0_X_XZ_map.*Bstarphi_XZ_map+J0_phi_XZ_map.*BstarX_XZ_map;
J0_B1_phi=J0_X_XZ_map.*BstarZ_XZ_map-J0_Z_XZ_map.*BstarX_XZ_map;
