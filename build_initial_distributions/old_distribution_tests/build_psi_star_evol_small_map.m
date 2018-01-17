clear all
load psi_star_evol.mat

psi_star_map_phi_evol=interp1(1:1001,psi_star_2D_evol_interp,1:101,'cubic');

save psi_star_map_phi_evol.mat psi_star_map_phi_evol