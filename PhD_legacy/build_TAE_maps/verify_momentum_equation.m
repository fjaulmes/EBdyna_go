% estimate perpendicular laplacian of A


RgA_X=zeros(sizeX,sizeZ);
gA_Z=zeros(sizeX,sizeZ);

for (x=3:sizeX-2)
    RgA_X(x,3:sizeZ-2)=(1/12)*(-A_XZ_map(x+2,(3:sizeZ-2))+A_XZ_map(x-2,(3:sizeZ-2)))+(2/3)*(A_XZ_map(x+1,(3:sizeZ-2))-A_XZ_map(x-1,(3:sizeZ-2)));
    gA_Z(x,3:sizeZ-2)=(1/12)*(-A_XZ_map(x,(3:sizeZ-2)+2)+A_XZ_map(x,(3:sizeZ-2)-2))+(2/3)*(A_XZ_map(x,(3:sizeZ-2)+1)-A_XZ_map(x,(3:sizeZ-2)-1));
end
RgA_X=RgA_X/DX;
RgA_X=RgA_X.*Rpos_XZsmall_map;
gA_Z=gA_Z/DX;


Rg2A_X=zeros(sizeX,sizeZ);
g2A_Z=zeros(sizeX,sizeZ);

for (x=3:sizeX-2)
    Rg2A_X(x,3:sizeZ-2)=(1/12)*(-RgA_X(x+2,(3:sizeZ-2))+RgA_X(x-2,(3:sizeZ-2)))+(2/3)*(RgA_X(x+1,(3:sizeZ-2))-RgA_X(x-1,(3:sizeZ-2)));
    g2A_Z(x,3:sizeZ-2)=(1/12)*(-gA_Z(x,(3:sizeZ-2)+2)+gA_Z(x,(3:sizeZ-2)-2))+(2/3)*(gA_Z(x,(3:sizeZ-2)+1)-gA_Z(x,(3:sizeZ-2)-1));
end
Rg2A_X=Rg2A_X/DX;
Rg2A_X=Rg2A_X./Rpos_XZsmall_map;
g2A_Z=g2A_Z/DX;

jppl=-(Rg2A_X+g2A_Z)/mu0;




ipsi_star_XZ_map=-iA_XZ_map.*Rpos_XZsmall_map.*bphi_XZsmall_map;
ipsi2D=ipsi_star_XZ_map;

% psi2D=psi_star_XZ_map;

igpsi_R=zeros(sizeX,sizeZ);
igpsi_Z=zeros(sizeX,sizeZ);

for (x=3:sizeX-2)
    igpsi_R(x,3:sizeZ-2)=(1/12)*(-ipsi2D(x+2,(3:sizeZ-2))+ipsi2D(x-2,(3:sizeZ-2)))+(2/3)*(ipsi2D(x+1,(3:sizeZ-2))-ipsi2D(x-1,(3:sizeZ-2)));
    igpsi_Z(x,3:sizeZ-2)=(1/12)*(-ipsi2D(x,(3:sizeZ-2)+2)+ipsi2D(x,(3:sizeZ-2)-2))+(2/3)*(ipsi2D(x,(3:sizeZ-2)+1)-ipsi2D(x,(3:sizeZ-2)-1));
end
igpsi_R=igpsi_R/DX;
igpsi_Z=igpsi_Z/DX;

iBstarX_XZ_map=-igpsi_Z./Rpos_XZsmall_map-(-(nTAE^2)*A_XZ_map.*bZ_XZ_small_map)./Rpos_XZsmall_map;
iBstarZ_XZ_map=igpsi_R./Rpos_XZsmall_map+(-(nTAE^2)*A_XZ_map.*bX_XZ_small_map)./Rpos_XZsmall_map;



RBsphi=Rpos_XZsmall_map.*Bstarphi_XZ_map;

gBphi_Z=zeros(sizeX,sizeZ);
gRBphi_X=zeros(sizeX,sizeZ);

for (x=3:sizeX-2)
    gBphi_Z(x,3:sizeZ-2)=(1/12)*(-Bstarphi_XZ_map(x,(3:sizeZ-2)+2)+Bstarphi_XZ_map(x,(3:sizeZ-2)-2))+(2/3)*(Bstarphi_XZ_map(x,(3:sizeZ-2)+1)-Bstarphi_XZ_map(x,(3:sizeZ-2)-1));
    gRBphi_X(x,3:sizeZ-2)=(1/12)*(-RBsphi(x+2,(3:sizeZ-2))+RBsphi(x-2,(3:sizeZ-2)))+(2/3)*(RBsphi(x+1,(3:sizeZ-2))-RBsphi(x-1,(3:sizeZ-2)));
end
gBphi_Z=gBphi_Z/DX;
gRBphi_X=gRBphi_X/DX;


J1_X_XZ_map=(-iBstarZ_XZ_map./Rpos_XZsmall_map+gBphi_Z)/mu0;
J1_Z_XZ_map=(-gRBphi_X./Rpos_XZsmall_map+iBstarX_XZ_map./Rpos_XZsmall_map)/mu0;



% 
% 
% RBphi=Rpos_XZsmall_map.*Bphi_XZsmall_map;
% 
% gBphi_Z=zeros(sizeX,sizeZ);
% gRBphi_X=zeros(sizeX,sizeZ);
% 
% for (x=3:sizeX-2)
%     gBphi_Z(x,3:sizeZ-2)=(1/12)*(-Bphi_XZsmall_map(x,(3:sizeZ-2)+2)+Bphi_XZsmall_map(x,(3:sizeZ-2)-2))+(2/3)*(Bphi_XZsmall_map(x,(3:sizeZ-2)+1)-Bphi_XZsmall_map(x,(3:sizeZ-2)-1));
%     gRBphi_X(x,3:sizeZ-2)=(1/12)*(-RBphi(x+2,(3:sizeZ-2))+RBphi(x-2,(3:sizeZ-2)))+(2/3)*(RBphi(x+1,(3:sizeZ-2))-RBphi(x-1,(3:sizeZ-2)));
% end
% gBphi_Z=gBphi_Z/DX;
% gRBphi_X=gRBphi_X/DX;
% 
% J0_X_XZ_map=(gBphi_Z)/mu0;
% J0_Z_XZ_map=(-gRBphi_X./Rpos_XZsmall_map)/mu0;


% gBZ_R=zeros(sizeX,sizeZ);
% gBX_Z=zeros(sizeX,sizeZ);
% 
% for (x=3:sizeX-2)
%     gBZ_R(x,3:sizeZ-2)=(1/12)*(-BpolZ_XZ_map(x+2,(3:sizeZ-2))+BpolZ_XZ_map(x-2,(3:sizeZ-2)))+(2/3)*(BpolZ_XZ_map(x+1,(3:sizeZ-2))-BpolZ_XZ_map(x-1,(3:sizeZ-2)));
%     gBX_Z(x,3:sizeZ-2)=(1/12)*(-BpolX_XZ_map(x,(3:sizeZ-2)+2)+BpolX_XZ_map(x,(3:sizeZ-2)-2))+(2/3)*(BpolX_XZ_map(x,(3:sizeZ-2)+1)-BpolX_XZ_map(x,(3:sizeZ-2)-1));   
% end
% gBZ_R=gBZ_R/DX;
% gBX_Z=gBX_Z/DX;
% J_phi_XZ_map=(gBZ_R-gBX_Z)/mu0;


% gBZ_R=zeros(sizeX,sizeZ);
% gBX_Z=zeros(sizeX,sizeZ);
% 
% for (x=3:sizeX-2)
%     gBZ_R(x,3:sizeZ-2)=(1/12)*(-BpolZ_initial_XZsmall_map(x+2,(3:sizeZ-2))+BpolZ_initial_XZsmall_map(x-2,(3:sizeZ-2)))+(2/3)*(BpolZ_initial_XZsmall_map(x+1,(3:sizeZ-2))-BpolZ_initial_XZsmall_map(x-1,(3:sizeZ-2)));
%     gBX_Z(x,3:sizeZ-2)=(1/12)*(-BpolX_initial_XZsmall_map(x,(3:sizeZ-2)+2)+BpolX_initial_XZsmall_map(x,(3:sizeZ-2)-2))+(2/3)*(BpolX_initial_XZsmall_map(x,(3:sizeZ-2)+1)-BpolX_initial_XZsmall_map(x,(3:sizeZ-2)-1));   
% end
% gBZ_R=gBZ_R/DX;
% gBX_Z=gBX_Z/DX;
% J0_phi_XZ_map=(gBZ_R-gBX_Z)/mu0;



gBZ_R=zeros(sizeX,sizeZ);
gBX_Z=zeros(sizeX,sizeZ);

for (x=3:sizeX-2)
    gBZ_R(x,3:sizeZ-2)=(1/12)*(-BstarZ_XZ_map(x+2,(3:sizeZ-2))+BstarZ_XZ_map(x-2,(3:sizeZ-2)))+(2/3)*(BstarZ_XZ_map(x+1,(3:sizeZ-2))-BstarZ_XZ_map(x-1,(3:sizeZ-2)));
    gBX_Z(x,3:sizeZ-2)=(1/12)*(-BstarX_XZ_map(x,(3:sizeZ-2)+2)+BstarX_XZ_map(x,(3:sizeZ-2)-2))+(2/3)*(BstarX_XZ_map(x,(3:sizeZ-2)+1)-BstarX_XZ_map(x,(3:sizeZ-2)-1));   
end
gBZ_R=gBZ_R/DX;
gBX_Z=gBX_Z/DX;
J1_phi_XZ_map=(gBZ_R-gBX_Z)/mu0;


gBZ_R=zeros(sizeX,sizeZ);
gBX_Z=zeros(sizeX,sizeZ);

for (x=3:sizeX-2)
    gBZ_R(x,3:sizeZ-2)=(1/12)*(-iBstarZ_XZ_map(x+2,(3:sizeZ-2))+iBstarZ_XZ_map(x-2,(3:sizeZ-2)))+(2/3)*(iBstarZ_XZ_map(x+1,(3:sizeZ-2))-iBstarZ_XZ_map(x-1,(3:sizeZ-2)));
    gBX_Z(x,3:sizeZ-2)=(1/12)*(-iBstarX_XZ_map(x,(3:sizeZ-2)+2)+iBstarX_XZ_map(x,(3:sizeZ-2)-2))+(2/3)*(iBstarX_XZ_map(x,(3:sizeZ-2)+1)-iBstarX_XZ_map(x,(3:sizeZ-2)-1));   
end
gBZ_R=gBZ_R/DX;
gBX_Z=gBX_Z/DX;
iJ1_phi_XZ_map=(gBZ_R-gBX_Z)/mu0;



J1_B0_X=J1_Z_XZ_map.*Bphi_tot_XZ_map-J1_phi_XZ_map.*BpolZ_XZ_map;
J1_B0_Z=-J1_X_XZ_map.*Bphi_tot_XZ_map+J1_phi_XZ_map.*BpolX_XZ_map;
J1_B0_phi=J1_X_XZ_map.*BpolZ_XZ_map-J1_Z_XZ_map.*BpolX_XZ_map;



J0_B1_X=J0_Z_XZ_map.*Bstarphi_XZ_map-J0_phi_XZ_map.*BstarZ_XZ_map;
J0_B1_Z=-J0_X_XZ_map.*Bstarphi_XZ_map+J0_phi_XZ_map.*BstarX_XZ_map;
J0_B1_phi=J0_X_XZ_map.*BstarZ_XZ_map-J0_Z_XZ_map.*BstarX_XZ_map;


% momentum map
mom_X_XZ_map=-(omega_TAE^2)*(mDT*ion_density_XZ_map).*(ksi_X);
mom_Z_XZ_map=-(omega_TAE^2)*(mDT*ion_density_XZ_map).*(ksi_Z);
mom_phi_XZ_map=-(omega_TAE^2)*(mDT*ion_density_XZ_map).*(ksi_phi);

% JxB force map
F_X_XZ_map=J0_B1_X+J1_B0_X;
F_Z_XZ_map=J0_B1_Z+J1_B0_Z;
F_phi_XZ_map=J0_B1_phi+J1_B0_phi;


% save momentum_data.mat mom_X_XZ_map mom_Z_XZ_map mom_phi_XZ_map F_X_XZ_map F_Z_XZ_map F_phi_XZ_map A_XZ_map 