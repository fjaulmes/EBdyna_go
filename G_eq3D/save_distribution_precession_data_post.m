% SAVENAME='initialG_flat2800_precession_stats.mat'
r_avg=r_avg';

SAVENAME

save(SAVENAME,'alphas_lambda','alphas_ejected','r_avg','q_avg','phidot_avg','omega_ae_avg','omega_ik_avg',...
    'omega_psi_avg','omega_precess_avg','omega_r_avg','delta_r_avg','r_vpll_plus_avg','r_vpll_minus_avg','omega_bounce',...
    'ALL_PASSING','ALL_TRAPPED','TRAPPED_PLUS','TRAPPED_MINUS','POTATOES','STAGNATION','CO_PASSING','COUNTER_PASSING',...
    'ALL_TRAPPED_POP','ALL_PASSING_POP','TRAPPED_PLUS_POP','TRAPPED_MINUS_POP','POTATOES_POP','STAGNATION_POP','CO_PASSING_POP','COUNTER_PASSING_POP');


