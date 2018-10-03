
NX=550;
NZ=4*NX;
FILENAME=strcat(FINESSE_FOLDER,'finesse_data.mat');
load(FILENAME);
number_of_data_points=size(finesse_data,1);


% SIGN_TOROIDAL_FIELD=sign(mean(finesse_data(1:end,end-13)))

% need to rescale the dimensions
%flipud of the Z map if we need to consider opposite toroidal field and current
Z_axis=SIGN_TOROIDAL_FIELD*Z_axis;

finesse_data(:,1)=a*finesse_data(:,1);
if SIGN_TOROIDAL_FIELD>0
    finesse_data(:,2)=a*finesse_data(:,2);
else
    finesse_data(:,2)=-a*finesse_data(:,2);
end
finesse_data(:,end)=SIGN_TOROIDAL_FIELD*finesse_data(:,end);

NR=4*NX;

inf_X=min(finesse_data(:,1));
sup_X=max(finesse_data(:,1));

inf_Z=min(finesse_data(:,2));
sup_Z=max(finesse_data(:,2));
DX=(sup_Z-inf_Z)/(4*NX-1);

% max_X=max(abs(inf_X),sup_X);
% margin necessary becaus of poor centering
max_Z=max(abs(inf_Z),sup_Z)+150*DX;

inf_X=-max_Z;
sup_X=max_Z;
inf_Z=-max_Z;
sup_Z=max_Z;

mid_X=(2*NX)+1;
mid_Z=(2*NX)+1;

DX=(sup_X-inf_X)/(4*NX-1);
%DZ=(2*max_Z)/(NZ-1);
DZ=DX;


X_scale=(0:NR-1)'*DX-(mid_X-1)*DX;
Z_scale=(0:NZ-1)'*DZ-(mid_Z-1)*DX;


X_axis_pos=round(X_axis/DX)+mid_X;
Z_axis_pos=round(Z_axis/DX)+mid_Z;

[XX,ZZ] = meshgrid(X_scale,Z_scale);


%radial geometry
for x=1:2*NX
    radial_r_value(x)=a*(x-1)/(NX-1);
end


% prepare axis
% for (R,Z) coordinates

R_inf_lim=R0-2*a-2*DX;

for (z=1:NZ)
    pos_z(z)=Z_scale(z);
    pos_x(z)=X_scale(z);
end

Rpos=R0+pos_x;



%Shafranov

psi2D=zeros(NR,NZ);
P_map=zeros(NR,NZ);
F_2_map=zeros(NR,NZ);
gPsi_X=zeros(NR,NZ);
gPsi_Z=zeros(NR,NZ);
dfgPsi_X=zeros(NR,NZ);
dgPsi_Z=zeros(NR,NZ);

%geometry

Raxis_map=zeros(NR,NZ);
r_radial_mask_map=zeros(NR,NZ);
r_value_map=zeros(NR,NZ);
gradX_r_map=zeros(NR,NZ);
gradZ_r_map=zeros(NR,NZ);
grad_r_square_map=zeros(NR,NZ);
xi_map=zeros(NR,NZ);
ki_map=zeros(NR,NZ);
q_map=zeros(NR,NZ);
q_psi_map=zeros(NR,NZ);
q_map_recalc=zeros(NR,NZ);
theta_map=zeros(NR,NZ);
gtheta_X=zeros(NR,NZ);
gtheta_Z=zeros(NR,NZ);

%magnetic field

BHpol_X=zeros(NR,NZ);
BHpol_Z=zeros(NR,NZ);


