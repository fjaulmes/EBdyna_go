close all;

figure(2);
set(gca,'FontSize',22);

ZF=0.9*INCREASE_RESOLUTION_FACTOR;
ZFO=round(ZF*ZOOM_RESIZE*NX);
%imagesc(X_scale_zoom(2*NX:9*NX),Z_scale_zoom(1*NX:9*NX),  psi_star_XZ_zoom_map(2*NX:9*NX,1*NX:9*NX)',[0.8*psi_star_max psi_star_max]);
%imagesc(X_scale_zoom(3*NX:8*NX),Z_scale_zoom(2*NX:8*NX),  Bstar_XZ_zoom_map(3*NX:8*NX,2*NX:8*NX)',[0.8*psi_star_max psi_star_max]);
% imagesc(X_scale_zoom(NX:ZF*2*NX),Z_scale_zoom(NX:ZF*2*NX),  integ_element_XZ_map(NX:ZF*2*NX,NX:ZF*2*NX)',[-8000 8000]);
imagesc(X_scale_zoom(mid_Z_zoom-ZFO:mid_Z_zoom+ZFO),Z_scale_zoom(mid_Z_zoom-ZFO:mid_Z_zoom+ZFO),integ_element_XZ_map(mid_Z_zoom-ZFO:mid_Z_zoom+ZFO,mid_Z_zoom-ZFO:mid_Z_zoom+ZFO)',[-20000 20000]);
colormap summer
brighten(0.5);

axis xy
hold on
xlim([-1.52 1.52]);
ylim([-1.5 1.5]);

xlabel('X (m)');
ylabel('Z (m)');


for(n=1:Ncontours-1)
%for(n=24:38)
    Nlength=Nlength_data(n);
    Nlength_half=Nhalf_data(n);
    %Nlength_end=Nend_data(n);

    %plot(x_value_map(n,:),z_value_map(n,:),'k')
    plot(x_value_map(n,1:Nlength_half),z_value_map(n,1:Nlength_half),'k')
%    plot(x_value_map(n,Nlength_half:Nlength),z_value_map(n,Nlength_half:Nlength),'k')
end
Nlength=Nlength_data(1);
Nlength_half=Nhalf_data(1);
plot(x_value_map(1,1:Nlength),z_value_map(1,1:Nlength),'r--')
Nlength=Nlength_data(Ncontours);
plot(x_value_map(Ncontours,1:Nlength),z_value_map(Ncontours,1:Nlength),'b')
plot(x_value_map(Ncontours,1:Nlength),z_value_map(Ncontours,1:Nlength),'k--')


plot(Line_ref_X,Line_ref_Z,'r')