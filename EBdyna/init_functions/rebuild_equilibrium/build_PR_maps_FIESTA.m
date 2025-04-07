
disp('reversing sign of psi profile due to change of clockwise convention !')
psi_FIESTA=FIESTA.prof.psi;
psi_FIESTA=[FIESTA.psi_axis   psi_FIESTA FIESTA.psi_boundary];
psi_n_FIESTA=psi_FIESTA-psi_FIESTA(1);
psi_n_FIESTA=(psi_n_FIESTA)/DPSI;

% surf_r_values=FIESTA.prof.volume./(2*pi*FIESTA.params.r0);
% surf_r_values=surf_r_values*FIESTA.params.area/surf_r_values(end);
surf_FIESTA=[0 FIESTA.prof.volume/(2*pi*0.894)   FIESTA.params.volume/(2*pi*0.894) ];

radial_r_values=sqrt(surf_FIESTA/pi);

r_scale=((1:Nradial)-1)/(Nradial-1);
% strange distribution of points in FIESTA when extremities are considered
% fiesta_scale=[0 (1:length(psi_FIESTA))/(length(psi_FIESTA)+1) 1];


psi_scale=interp1(psi_n_FIESTA,psi_FIESTA,r_scale,'makima');
psi_n_scale=psi_scale-(psi_scale(1));
psi_n_scale=(psi_n_scale)/DPSI;

q_profile_METIS=interp1(post.profil0d.temps,post.profil0d.qjli,TMETIS);
q_FIESTA=[FIESTA.params.q0 FIESTA.prof.q   q_profile_METIS(end) ];    
q_initial_profile=interp1(psi_n_FIESTA,q_FIESTA,r_scale,'makima');
tor_flux_scale=cumtrapz(psi_scale,q_initial_profile);

tor_flux_n_profile=tor_flux_scale/tor_flux_scale(end);
rho_tor_scale=sqrt(tor_flux_n_profile);

F_FIESTA=[FIESTA.params.b0*FIESTA.params.r0  FIESTA.prof.f   FVAC];
Fdia_profile=interp1(psi_n_FIESTA,F_FIESTA,r_scale,'makima');
F2_profile=Fdia_profile.^2;

radial_r_value_flux=interp1(psi_n_FIESTA,radial_r_values,r_scale,'makima');

vol_FIESTA=[0  FIESTA.prof.volume   FIESTA.params.volume];

volume_flux=interp1(psi_n_FIESTA,vol_FIESTA,r_scale,'makima');
surf_flux=interp1(psi_n_FIESTA,surf_FIESTA,r_scale,'makima');

edge_pressure=2.*eV*interp1(post.zerod.temps,post.zerod.nebord.*post.zerod.tebord,TMETIS);
P_FIESTA=[FIESTA.params.P0 FIESTA.prof.pressure   edge_pressure];
P_initial_profile=interp1(psi_n_FIESTA,P_FIESTA,r_scale,'makima');

% reversed toroidal direction in tokamak coordinates
Iaxis=FIESTA.ip;
P0=P_initial_profile(1);
X_axis=FIESTA.params.r0-R0;
Z_axis=FIESTA.params.z0;
elongation=FIESTA.params.kappa;

psi_global_recalc=abs(psi_scale(1)-psi_scale(end));

if abs(psi_global_recalc-psi_global)>0.0001
    warning('THERE IS aN ISSUE WITH PSI_GLOBAL')
    warning('recheck the psi map for double nulls : they can cause FIESTA to give incorrect axis values')
end