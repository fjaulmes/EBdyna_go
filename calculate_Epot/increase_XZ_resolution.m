
% NX=NX*2;
% NZ=4*NX;
% mid_X=2*NX+2;
% mid_Z=2*NX+2;
% 
% inf_Z=min(min(Z_PR_map));
% sup_Z=max(max(Z_PR_map));
% 
% max_Z=max(abs(inf_Z),sup_Z);
% 
% inf_X=-max_Z;
% sup_X=max_Z;
% inf_Z=-max_Z;
% sup_Z=max_Z

INCREASE_RESOLUTION_FACTOR=3;

% X_scale=((0:4*NX-1)-mid_Z)'*DX+X_axis;
% Z_scale=((0:4*NX-1)-mid_Z)'*DX+Z_axis;
[XX,ZZ] = meshgrid(X_scale,Z_scale);

X_scale_data=reshape(XX(:,:),NZ*NZ,1);
Z_scale_data=reshape(ZZ(:,:),NZ*NZ,1);

RR=Rpos_PR_map(:,1:size_r)-R0;


DZ_zoom=DX/INCREASE_RESOLUTION_FACTOR;
DX_zoom=DZ_zoom;


ZOOM_RESIZE=1.15

mid_X_zoom=round(INCREASE_RESOLUTION_FACTOR*ZOOM_RESIZE*NX)+1;
mid_Z_zoom=round(INCREASE_RESOLUTION_FACTOR*ZOOM_RESIZE*NX)+1;

NZ_zoom=round(INCREASE_RESOLUTION_FACTOR*ZOOM_RESIZE*2*NX);

X_scale_zoom=((0:NZ_zoom-1)-mid_X_zoom+1)'*DX_zoom;
Z_scale_zoom=((0:NZ_zoom-1)-mid_Z_zoom+1)'*DX_zoom;



[XX_zoom,ZZ_zoom] = meshgrid(X_scale_zoom,Z_scale_zoom);

