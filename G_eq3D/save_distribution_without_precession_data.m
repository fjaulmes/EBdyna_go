



%PART_POP=find(~alphas_ejected);
PART_POP=find((~alphas_ejected).*(~isnan(alphas_vpll)));

alphas_pos_x=alphas_pos_x(PART_POP);
alphas_pos_z=alphas_pos_z(PART_POP);
alphas_pos_phi=alphas_pos_phi(PART_POP);
v_X=v_X(PART_POP);
v_Z=v_Z(PART_POP);
v_phi=v_phi(PART_POP);
alphas_vpll=alphas_vpll(PART_POP);
alphas_mm=alphas_mm(PART_POP);
alphas_psi=alphas_psi(PART_POP);
alphas_Ekin=alphas_Ekin(PART_POP);
alphas_pphi0=alphas_pphi0(PART_POP);
alphas_momentum=alphas_momentum(PART_POP);

Nalphas_simulated=length(alphas_Ekin);

save(SAVENAME,'mHe','ZHe','alphas_momentum','alphas_pos_x','alphas_pos_z','alphas_pos_phi','v_X','v_Z','v_phi','alphas_vpll','alphas_mm','alphas_psi','alphas_Ekin','alphas_pphi0','Nalphas_simulated');
