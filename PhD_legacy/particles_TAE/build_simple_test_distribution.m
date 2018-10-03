DPOP=100

load('initial_MBD_600_pre_collapse_all.mat')
PART_POP=find(alphas_psi<200);
alphas_Ekin=alphas_Ekin(PART_POP(1:DPOP:end));
alphas_mm=alphas_mm(PART_POP(1:DPOP:end));
alphas_weight=alphas_weight(PART_POP(1:DPOP:end));
alphas_pphi0=alphas_pphi0(PART_POP(1:DPOP:end));
alphas_pos_x=alphas_pos_x(PART_POP(1:DPOP:end));
alphas_pos_z=alphas_pos_z(PART_POP(1:DPOP:end));
alphas_pos_phi=alphas_pos_phi(PART_POP(1:DPOP:end));
alphas_vpll=alphas_vpll(PART_POP(1:DPOP:end));
alphas_psi=alphas_psi(PART_POP(1:DPOP:end));
v_X=v_X(PART_POP(1:DPOP:end));
v_Z=v_Z(PART_POP(1:DPOP:end));
v_phi=v_phi(PART_POP(1:DPOP:end));
pos_X_gc=pos_X_gc(PART_POP(1:DPOP:end));
pos_Z_gc=pos_Z_gc(PART_POP(1:DPOP:end));

Nalphas_simulated=length(alphas_Ekin)

save('initial_MBD_600_pre_collapse_super_light.mat','mHe', 'ZHe','Nalphas_simulated','alphas_weight','alphas_Ekin','alphas_pphi0','alphas_mm', 'alphas_pphi0', 'pos_X_gc','pos_Z_gc','alphas_pos_x', 'alphas_pos_z', 'alphas_pos_phi', 'alphas_vpll', 'alphas_psi', 'v_X', 'v_Z', 'v_phi');


load('initial_MBD_600_precession_stats_all.mat')

vpll_avg=vpll_avg(PART_POP(1:DPOP:end));
r_avg=r_avg(PART_POP(1:DPOP:end));
q_avg=q_avg(PART_POP(1:DPOP:end));
alphas_lambda=alphas_lambda(PART_POP(1:DPOP:end));
ALL_TRAPPED_POP=ALL_TRAPPED_POP(PART_POP(1:DPOP:end));
ALL_PASSING_POP=ALL_PASSING_POP(PART_POP(1:DPOP:end));
TRAPPED_MINUS_POP=TRAPPED_MINUS_POP(PART_POP(1:DPOP:end));
TRAPPED_PLUS_POP=TRAPPED_PLUS_POP(PART_POP(1:DPOP:end));
POTATOES_POP=POTATOES_POP(PART_POP(1:DPOP:end));
STAGNATION_POP=STAGNATION_POP(PART_POP(1:DPOP:end));
CO_PASSING_POP=CO_PASSING_POP(PART_POP(1:DPOP:end));
COUNTER_PASSING_POP=COUNTER_PASSING_POP(PART_POP(1:DPOP:end));
omega_psi_avg=omega_psi_avg(PART_POP(1:DPOP:end));
r_vpll_minus_avg=r_vpll_minus_avg(PART_POP(1:DPOP:end));
r_vpll_plus_avg=r_vpll_plus_avg(PART_POP(1:DPOP:end));
omega_r_avg=omega_r_avg(PART_POP(1:DPOP:end));
delta_r_avg=delta_r_avg(PART_POP(1:DPOP:end));
omega_precess_avg=omega_precess_avg(PART_POP(1:DPOP:end));
omega_bounce=omega_bounce(PART_POP(1:DPOP:end));
phidot_avg=phidot_avg(PART_POP(1:DPOP:end));


ALL_TRAPPED=find(ALL_TRAPPED_POP);
ALL_PASSING=find(ALL_PASSING_POP);
TRAPPED_PLUS=find(TRAPPED_PLUS_POP);
TRAPPED_MINUS=find(TRAPPED_MINUS_POP);
POTATOES=find(POTATOES_POP);
STAGNATION=find(STAGNATION_POP);
CO_PASSING=find(CO_PASSING_POP);
COUNTER_PASSING=find(COUNTER_PASSING_POP);


save('initial_MBD_600_precession_stats_super_light.mat','mHe', 'ZHe', 'vpll_avg', 'r_avg', 'q_avg',  'alphas_lambda', 'ALL_TRAPPED_POP', 'ALL_PASSING_POP', 'TRAPPED_MINUS_POP', 'TRAPPED_PLUS_POP', 'POTATOES_POP', 'STAGNATION_POP', 'CO_PASSING_POP', 'COUNTER_PASSING_POP', 'ALL_TRAPPED', 'ALL_PASSING', 'TRAPPED_PLUS', 'TRAPPED_MINUS', 'POTATOES', 'STAGNATION', 'CO_PASSING', 'COUNTER_PASSING', 'omega_psi_avg', 'r_vpll_minus_avg', 'r_vpll_plus_avg', 'omega_r_avg', 'delta_r_avg', 'omega_precess_avg', 'omega_bounce', 'phidot_avg')