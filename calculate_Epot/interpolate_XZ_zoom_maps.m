


% equilibrium maps

%Bphi_XZ_zoom_map=interp2(X_scale,Z_scale,sqrt(Bphi_XZ_map.^2+BHpol_X.^2+BHpol_Z.^2),XX,ZZ,'cubic');
Bphi_XZ_zoom_map=interp2(X_scale,Z_scale,Bphi_XZ_map,XX_zoom,ZZ_zoom,'cubic');

% BHpol_X_XZ_zoom_map=interp2(X_scale,Z_scale,BHpol_X_map,XX,ZZ,'cubic');
% BHpol_Z_XZ_zoom_map=interp2(X_scale,Z_scale,BHpol_Z_map,XX,ZZ,'cubic');


Rpos_XZ_zoom_map=interp2(X_scale,Z_scale,Rpos_map,XX_zoom,ZZ_zoom,'cubic');

%clear Rpos_map;

theta_XZ_zoom_map=interp2(X_scale,Z_scale,theta_XZ_map,XX_zoom,ZZ_zoom,'cubic');

%clear theta_XZ_map;

radial_XZ_zoom_map=interp2(X_scale,Z_scale,radial_XZ_map,XX_zoom,ZZ_zoom,'cubic');

% X_data=reshape(XX,NZ*NZ,1);
% Z_data=reshape(ZZ,NZ*NZ,1);
% initial_data=reshape(radial_XZ_map,NZ*NZ,1);
% 
% radial_XZ_zoom_map=griddata(X_data,Z_data,initial_data,XX_zoom,ZZ_zoom,'cubic');
% radial_XZ_zoom_map=radial_XZ_zoom_map';
% radial_XZ_zoom_map(isnan(radial_XZ_zoom_map))=0;

%clear radial_XZ_map;



% psi_star_XZ_zoom_map=interp2(X_scale,Z_scale,psi_star_RZ_map,XX,ZZ,'cubic');
% 
% clear psi_star_RZ_map;


% Bstar fields

% Bstar_XZ_zoom_map=interp2(X_scale,Z_scale,Bstar_RZ_map,XX,ZZ,'cubic');
% Bstar_XZ_zoom_map(isnan(Bstar_XZ_zoom_map))=0;
% 
% clear Bstar_RZ_map;
% 
% Bstar_X_XZ_zoom_map=interp2(X_scale,Z_scale,Bstar_X_RZ_map,XX,ZZ,'cubic');
% Bstar_X_XZ_zoom_map(isnan(Bstar_X_XZ_zoom_map))=0;
% 
% clear Bstar_X_RZ_map;
% 
% Bstar_Z_XZ_zoom_map=interp2(X_scale,Z_scale,Bstar_Z_RZ_map,XX,ZZ,'cubic');
% Bstar_Z_XZ_zoom_map(isnan(Bstar_Z_XZ_zoom_map))=0;
% 
% clear Bstar_Z_RZ_map;
% 

