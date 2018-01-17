
close all;

mid_X=INCREASE_RESOLUTION_FACTOR*NX+INCREASE_RESOLUTION_FACTOR;
mid_Z=INCREASE_RESOLUTION_FACTOR*NX+INCREASE_RESOLUTION_FACTOR;

NZ=2*INCREASE_RESOLUTION_FACTOR*NX;


ZOOM_FACTOR=3.5;

figure(1);
imagesc(X_scale_zoom(mid_Z-ZOOM_FACTOR*NX:mid_Z+ZOOM_FACTOR*NX),Z_scale_zoom(mid_Z-ZOOM_FACTOR*NX:mid_Z+ZOOM_FACTOR*NX),integ_element_XZ_map(mid_Z-ZOOM_FACTOR*NX:mid_Z+ZOOM_FACTOR*NX,mid_Z-ZOOM_FACTOR*NX:mid_Z+ZOOM_FACTOR*NX)');
axis xy;
hold on

n=1;
Nlength=Nlength_data(n);
plot(x_value_map(n,1:Nlength),z_value_map(n,1:Nlength),'k')

Nlength=Nlength_data(n+1);
plot(x_value_map(n+1,1:Nlength),z_value_map(n+1,1:Nlength),'k')
Nlength=Nlength_data(n+2);
plot(x_value_map(n+2,1:Nlength),z_value_map(n+2,1:Nlength),'k')

n=Ncontours-1;
Nlength=Nlength_data(n-1);
plot(x_value_map(n-1,1:Nlength),z_value_map(n-1,1:Nlength),'k')
Nlength=Nlength_data(n-2);
plot(x_value_map(n-2,1:Nlength),z_value_map(n-2,1:Nlength),'k')

n=Ncontours;
Nlength=Nlength_data(n);
plot(x_value_map(n,1:Nlength),z_value_map(n,1:Nlength),'k')


figure(2);
imagesc(X_scale_zoom(mid_Z-ZOOM_FACTOR*NX:mid_Z+ZOOM_FACTOR*NX),Z_scale_zoom(mid_Z-ZOOM_FACTOR*NX:mid_Z+ZOOM_FACTOR*NX),psi_star_dot_XZ_zoom_map(mid_Z-ZOOM_FACTOR*NX:mid_Z+ZOOM_FACTOR*NX,mid_Z-ZOOM_FACTOR*NX:mid_Z+ZOOM_FACTOR*NX)');
axis xy;
hold on;

n=Ncontours;
Nlength=Nlength_data(n);
plot(x_value_map(n,1:Nlength),z_value_map(n,1:Nlength),'k')




figure(4);
imagesc(X_scale_zoom(mid_Z-ZOOM_FACTOR*NX:mid_Z+ZOOM_FACTOR*NX),Z_scale_zoom(mid_Z-ZOOM_FACTOR*NX:mid_Z+ZOOM_FACTOR*NX),Bstar_XZ_zoom_map(mid_Z-ZOOM_FACTOR*NX:mid_Z+ZOOM_FACTOR*NX,mid_Z-ZOOM_FACTOR*NX:mid_Z+ZOOM_FACTOR*NX)');
axis xy;
hold on;


figure(5);
imagesc(X_scale_zoom(mid_Z-ZOOM_FACTOR*NX:mid_Z+ZOOM_FACTOR*NX),Z_scale_zoom(mid_Z-ZOOM_FACTOR*NX:mid_Z+ZOOM_FACTOR*NX),psi_star_XZ_zoom_map(mid_Z-ZOOM_FACTOR*NX:mid_Z+ZOOM_FACTOR*NX,mid_Z-ZOOM_FACTOR*NX:mid_Z+ZOOM_FACTOR*NX)');
axis xy;
hold on;



ZOOM_FACTOR=2.7;

figure(3);
imagesc(X_scale_zoom(mid_Z-ZOOM_FACTOR*NX:mid_Z+ZOOM_FACTOR*NX),Z_scale_zoom(mid_Z-ZOOM_FACTOR*NX:mid_Z+ZOOM_FACTOR*NX),E_potential_XZ_zoom_map(mid_Z-ZOOM_FACTOR*NX:mid_Z+ZOOM_FACTOR*NX,mid_Z-ZOOM_FACTOR*NX:mid_Z+ZOOM_FACTOR*NX)');
axis xy;
hold on;

n=Ncontours;
Nlength=Nlength_data(n);
plot(x_value_map(n,1:Nlength),z_value_map(n,1:Nlength),'k')


% figure(6);
% imagesc(X_scale_zoom(mid_Z-ZOOM_FACTOR*NX:mid_Z+ZOOM_FACTOR*NX),Z_scale_zoom(mid_Z-ZOOM_FACTOR*NX:mid_Z+ZOOM_FACTOR*NX),Bstar_2_XZ_zoom_map(mid_Z-ZOOM_FACTOR*NX:mid_Z+ZOOM_FACTOR*NX,mid_Z-ZOOM_FACTOR*NX:mid_Z+ZOOM_FACTOR*NX)');
% axis xy;
% hold on;

% mid_X=INCREASE_RESOLUTION_FACTOR*NX+INCREASE_RESOLUTION_FACTOR;
% mid_Z=INCREASE_RESOLUTION_FACTOR*NX+INCREASE_RESOLUTION_FACTOR;
% 
% NZ=2*INCREASE_RESOLUTION_FACTOR*NX;