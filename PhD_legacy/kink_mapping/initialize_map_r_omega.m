
% Dimensions of the region inside the mixing radius
% in helical coordinates (r,omega)

% For the equation of psi* over time
% we need r to be between 0 and 1
PRECISE_MESH=4.0

size_r=xMixing_radius+8;
NR=size_r*PRECISE_MESH;


% Dr=mean(mean(dist_surf_PR_map(:,1:size_r)));
% Dr=1/(xPsih_zero-1);
% scale_r=Dr*[0:size_r-1];
scale_r=interp1((0:size_r-1)/(size_r-1),radial_r_value_flux(1:size_r),(0:NR-1)/(NR-1));
Psih_final=interp1((0:size_r-1)/(size_r-1),Psih_final(1:size_r),(0:NR-1)/(NR-1));
Psih=interp1((0:size_r-1)/(size_r-1),Psih(1:size_r),(0:NR-1)/(NR-1));

surf_flux_precise=interp1((0:size_r-1)/(size_r-1),surf_flux(1:size_r),(0:NR-1)/(NR-1));
P_ini_precise=interp1((0:size_r-1)/(size_r-1),P_initial_profile(1:size_r),(0:NR-1)/(NR-1));
scale_X_precise=interp1((0:size_r-1)/(size_r-1),scale_X(1:size_r),(0:NR-1)/(NR-1));

Dr_avg(1)=0;
for (r=2:NR)
    Dr_avg(r)=scale_r(r)-scale_r(r-1);
end

% Dr=mean(Dr_avg(1:size_r));

% we need w to be between 0 and pi

size_omega=pi;
Nomega=round((NP+1)/2);
Domega=size_omega/(Nomega-1);

half_Nomega=round(Nomega/2);

psi_star_2D=zeros(NR,Nomega);
scale_omega=[0:Nomega-1]*Domega;

% r_PR_map=sqrt((Z_PR_map-Z_axis).^2+(Rpos_PR_map-R0-X_axis).^2);
%dl_avg=2*pi*scale_r/(Nomega-1);
dl_avg_raw=mean(dl_PR_map(1:Nomega-1,:),1);
dl_avg_raw=interp1((0:size_r-1)/(size_r-1),dl_avg_raw(1:size_r),(0:NR-1)/(NR-1));

dl_avg(1)=0.5*dl_avg_raw(1);
for (r=2:NR)
    dl_avg(r)=0.5*(dl_avg_raw(r)+dl_avg_raw(r-1));
end

% sqrtg_S_avg=mean((sqrtg_B_PR_map(1:Nomega-1,:)).*Rpos_PR_map(1:Nomega-1,:).*BHpol_PR_map(1:Nomega-1,:)./(Btor_PR_map(1:Nomega-1,:)),1);


% for (p=1:NP)
% scale_X_PR_map(p,:)=scale_X;
% end

%r_value_r_omega=griddata(scale_X_PR_map(:,1:xPsih_zero)',theta_PR_map(:,1:xPsih_zero)',r_PR_map(:,1:xPsih_zero)',scale_r',scale_omega);
