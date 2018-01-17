

NR=4*NX;


inf_X=min(finesse_data(:,1));
sup_X=max(finesse_data(:,1));

inf_Z=min(finesse_data(:,2));
sup_Z=max(finesse_data(:,2));

max_X=max(abs(inf_X),sup_X);
max_Z=max(abs(inf_Z),sup_Z);

inf_X=-max_Z;
sup_X=max_Z;
inf_Z=-max_Z-Z_axis;
sup_Z=max_Z-Z_axis;

mid_X=(2*NX)+1;
mid_Z=(2*NX)+1;

DX=(sup_X-inf_X)/(4*NX-1);
%DZ=(2*max_Z)/(NZ-1);
DZ=DX;


X_scale=(0:NR-1)'*DX-(mid_X-1)*DX;
Z_scale=(0:NZ-1)'*DZ-(mid_Z-1)*DX-2*Z_axis;


X_axis_pos=round(X_axis/DX)+mid_X;

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



