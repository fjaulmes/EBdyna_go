
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
div_B=gRB_R./Rpos_XZsmall_map+gB_Z+iBstarphi_XZ_map./Rpos_XZsmall_map;




% check div(J0)=0

gRJ_R=zeros(sizeX,sizeZ);
gJ_Z=zeros(sizeX,sizeZ);
div_J0=zeros(sizeX,sizeZ);

RJ_X=J0_X_XZ_map.*Rpos_XZsmall_map;
J_Z=J0_Z_XZ_map;

for (x=3:sizeX-2)
    gRJ_R(x,3:sizeZ-2)=(1/12)*(-RJ_X(x+2,(3:sizeZ-2))+RJ_X(x-2,(3:sizeZ-2)))+(2/3)*(RJ_X(x+1,(3:sizeZ-2))-RJ_X(x-1,(3:sizeZ-2)));
    gJ_Z(x,3:sizeZ-2)=(1/12)*(-J_Z(x,(3:sizeZ-2)+2)+J_Z(x,(3:sizeZ-2)-2))+(2/3)*(J_Z(x,(3:sizeZ-2)+1)-J_Z(x,(3:sizeZ-2)-1));
end
gRJ_R=gRJ_R/DX;
gJ_Z=gJ_Z/DX;

% approximate divergence, without toroidal derivative
% div_J=gRJ_R./Rpos_XZsmall_map+gJ_Z;
div_J0(3:sizeX-2,3:sizeZ-2)=gRJ_R(3:sizeX-2,3:sizeZ-2)./Rpos_XZsmall_map(3:sizeX-2,3:sizeZ-2)+gJ_Z(3:sizeX-2,3:sizeZ-2);




% check div(J1)=0

gRJ_R=zeros(sizeX,sizeZ);
gJ_Z=zeros(sizeX,sizeZ);

RJ_X=J1_X_XZ_map.*Rpos_XZsmall_map;
J_Z=J1_Z_XZ_map;

for (x=3:sizeX-2)
    gRJ_R(x,3:sizeZ-2)=(1/12)*(-RJ_X(x+2,(3:sizeZ-2))+RJ_X(x-2,(3:sizeZ-2)))+(2/3)*(RJ_X(x+1,(3:sizeZ-2))-RJ_X(x-1,(3:sizeZ-2)));
    gJ_Z(x,3:sizeZ-2)=(1/12)*(-J_Z(x,(3:sizeZ-2)+2)+J_Z(x,(3:sizeZ-2)-2))+(2/3)*(J_Z(x,(3:sizeZ-2)+1)-J_Z(x,(3:sizeZ-2)-1));
end
gRJ_R=gRJ_R/DX;
gJ_Z=gJ_Z/DX;

% approximate divergence, without toroidal derivative
% div_J=gRJ_R./Rpos_XZsmall_map+gJ_Z;
div_J1=gRJ_R./Rpos_XZsmall_map+gJ_Z+iJ1_phi_XZ_map./Rpos_XZsmall_map;
