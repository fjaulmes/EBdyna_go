MIN_BSTAR_VALUE=2e-6;
MIN_BSTAR_COORD_VALUE=1e-6;



gPsi_R=zeros(NZ_zoom,NZ_zoom);
gPsi_Z=zeros(NZ_zoom,NZ_zoom);
Bstar_X_XZ_zoom_map=zeros(NZ_zoom,NZ_zoom);
Bstar_Z_XZ_zoom_map=zeros(NZ_zoom,NZ_zoom);
Bstar_X_PR_map=zeros(NP,Nradial);
Bstar_Z_PR_map=zeros(NP,Nradial);


gPsi_R((minX:maxX),(minZ:maxZ))=(1/12)*(-psi_star_XZ_zoom_map((minX:maxX)+2,(minZ:maxZ))+psi_star_XZ_zoom_map((minX:maxX)-2,(minZ:maxZ)))+(2/3)*(psi_star_XZ_zoom_map((minX:maxX)+1,(minZ:maxZ))-psi_star_XZ_zoom_map((minX:maxX)-1,(minZ:maxZ)));
gPsi_Z((minX:maxX),(minZ:maxZ))=(1/12)*(-psi_star_XZ_zoom_map((minX:maxX),(minZ:maxZ)+2)+psi_star_XZ_zoom_map((minX:maxX),(minZ:maxZ)-2))+(2/3)*(psi_star_XZ_zoom_map((minX:maxX),(minZ:maxZ)+1)-psi_star_XZ_zoom_map((minX:maxX),(minZ:maxZ)-1));
gPsi_R=gPsi_R.*mask_XZ_zoom_map_reconnection;
gPsi_Z=gPsi_Z.*mask_XZ_zoom_map_reconnection;

gPsi_R=gPsi_R/DX_zoom;
gPsi_Z=gPsi_Z/DZ_zoom;

Bstar_X_XZ_zoom_map=-gPsi_Z./Rpos_XZ_zoom_map;
Bstar_Z_XZ_zoom_map=gPsi_R./Rpos_XZ_zoom_map;
Bstar_X_XZ_zoom_map=Bstar_X_XZ_zoom_map.*mask_XZ_zoom_map_reconnection;
Bstar_Z_XZ_zoom_map=Bstar_Z_XZ_zoom_map.*mask_XZ_zoom_map_reconnection;





Bstar_X=max(abs(Bstar_X_XZ_zoom_map),MIN_BSTAR_COORD_VALUE);
Bstar_Z=max(abs(Bstar_Z_XZ_zoom_map),MIN_BSTAR_COORD_VALUE);
Bstar_X_XZ_zoom_map=sign(Bstar_X_XZ_zoom_map).*Bstar_X;
Bstar_Z_XZ_zoom_map=sign(Bstar_Z_XZ_zoom_map).*Bstar_Z;


Bstar_XZ_zoom_map=sqrt(Bstar_X_XZ_zoom_map.^2+Bstar_Z_XZ_zoom_map.^2);
Bstar_XZ_zoom_map=max(Bstar_XZ_zoom_map,MIN_BSTAR_VALUE);

Bstar_XZ_zoom_map=Bstar_XZ_zoom_map.*mask_XZ_zoom_map_reconnection;
Bstar_XZ_zoom_map=Bstar_XZ_zoom_map+abs(Bstar_initial_map).*(1-mask_XZ_zoom_map_reconnection);

Bstar_X_RZ_map=interp2(X_scale_zoom,Z_scale_zoom,Bstar_X_XZ_zoom_map,XX,ZZ,'cubic');
Bstar_Z_RZ_map=interp2(X_scale_zoom,Z_scale_zoom,Bstar_Z_XZ_zoom_map,XX,ZZ,'cubic');
Bstar_X_RZ_map(isnan(Bstar_X_RZ_map))=0;
Bstar_Z_RZ_map(isnan(Bstar_Z_RZ_map))=0;



if DISPLAY_OUTPUTS==1
    ZOOM_FACTOR=0.8*INCREASE_RESOLUTION_FACTOR;
    figure(1)
    Xdisp_inf=round(mid_Z_zoom-ZOOM_FACTOR*ZOOM_RESIZE*NX);
    Zdisp_inf=round(mid_Z_zoom-ZOOM_FACTOR*ZOOM_RESIZE*NX);
    Xdisp_sup=round(mid_Z_zoom+ZOOM_FACTOR*ZOOM_RESIZE*NX);
    Zdisp_sup=round(mid_Z_zoom+ZOOM_FACTOR*ZOOM_RESIZE*NX);
    imagesc(X_scale_zoom(Xdisp_inf:Xdisp_sup),Z_scale_zoom(Zdisp_inf:Zdisp_sup),Bstar_XZ_zoom_map(Xdisp_inf:Xdisp_sup,Zdisp_inf:Zdisp_sup)',[-0.04 0.04]);
    axis square xy;
    title('B_* tot')
    xlabel('X');
    ylabel('Z');
    colorbar;
    pause(0.1);
end