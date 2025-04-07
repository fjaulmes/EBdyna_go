
NG=200;
NZ=6*NG;

NR=3*NG;

FRACTION_XDATA=0.56
FRACTION_ZDATA=0.7

inf_X=(min(FIESTA.scale_R)-R0)*FRACTION_XDATA+R0;
sup_X=(max(FIESTA.scale_R)-R0)*FRACTION_XDATA+R0;

inf_Z=min(FIESTA.scale_Z)*FRACTION_ZDATA;
sup_Z=max(FIESTA.scale_Z)*FRACTION_ZDATA;
DZ=(sup_Z-inf_Z)/(NZ-1);

% max_X=max(abs(inf_X),sup_X);
% margin necessary becaus of poor centering


mid_X=(0.5*NR)+1;
mid_Z=(0.5*NZ)+1;

DX=DZ;
%DZ=(2*max_Z)/(NZ-1);
DZ=DX;


X_scale=(0:NR-1)'*DX-(mid_X-1)*DX;
Z_scale=(0:NZ-1)'*DZ-(mid_Z-1)*DX;


X_axis_pos=round(X_axis/DX)+mid_X;
Z_axis_pos=round(Z_axis/DX)+mid_Z;

[XX,ZZ] = meshgrid(X_scale,Z_scale);


Rpos=X_scale+R0;
Rpos_XZ_map=repmat(Rpos,1,NZ);
