%%
NB_PROCESS=16
load('initial_NBI_1MEV_D_distribution1.mat')
% precession stats
SAVENAME='../final_NBI_1MEV_D_precession_stats_all.mat'

n=1;
INPUTNAME=strcat('final_NBI_1MEV_D_precession_stats',num2str(n),'.mat')
load(INPUTNAME)
% r_avg=r_avg';

for n=2:NB_PROCESS
    r_avg_tot=r_avg;
    q_avg_tot=q_avg;
%     alphas_kappa_tot=alphas_kappa;
    alphas_lambda_tot=alphas_lambda;
    ALL_TRAPPED_POP_tot=ALL_TRAPPED_POP;
    ALL_PASSING_POP_tot=ALL_PASSING_POP;
    TRAPPED_MINUS_POP_tot=TRAPPED_MINUS_POP;
    TRAPPED_PLUS_POP_tot=TRAPPED_PLUS_POP;
    POTATOES_POP_tot=POTATOES_POP;
    STAGNATION_POP_tot=STAGNATION_POP;
    CO_PASSING_POP_tot=CO_PASSING_POP;
    COUNTER_PASSING_POP_tot=COUNTER_PASSING_POP;
    omega_ik_avg_tot=omega_ik_avg;
    omega_psi_avg_tot=omega_psi_avg;
    r_vpll_minus_avg_tot=r_vpll_minus_avg;
    r_vpll_plus_avg_tot=r_vpll_plus_avg;
    omega_r_avg_tot=omega_r_avg;
    delta_r_avg_tot=delta_r_avg;
    omega_precess_avg_tot=omega_precess_avg;
    omega_bounce_tot=omega_bounce;
    phidot_avg_tot=phidot_avg;
    alphas_ejected_tot=alphas_ejected;
    
    INPUTNAME=strcat('final_NBI_1MEV_D_precession_stats',num2str(n),'.mat')
    load(INPUTNAME)
    
        
%     r_avg=r_avg';
    
    r_avg_p2=r_avg;
    q_avg_p2=q_avg;
%     alphas_kappa_p2=alphas_kappa;
    alphas_lambda_p2=alphas_lambda;
    ALL_TRAPPED_POP_p2=ALL_TRAPPED_POP;
    ALL_PASSING_POP_p2=ALL_PASSING_POP;
    TRAPPED_MINUS_POP_p2=TRAPPED_MINUS_POP;
    TRAPPED_PLUS_POP_p2=TRAPPED_PLUS_POP;
    POTATOES_POP_p2=POTATOES_POP;
    STAGNATION_POP_p2=STAGNATION_POP;
    CO_PASSING_POP_p2=CO_PASSING_POP;
    COUNTER_PASSING_POP_p2=COUNTER_PASSING_POP;
    omega_ik_avg_p2=omega_ik_avg;
    omega_psi_avg_p2=omega_psi_avg;
    r_vpll_minus_avg_p2=r_vpll_minus_avg;
    r_vpll_plus_avg_p2=r_vpll_plus_avg;
    omega_r_avg_p2=omega_r_avg;
    delta_r_avg_p2=delta_r_avg;
    omega_precess_avg_p2=omega_precess_avg;
    omega_bounce_p2=omega_bounce;
    phidot_avg_p2=phidot_avg;
    alphas_ejected_p2=alphas_ejected;
    
    
    
    r_avg=[r_avg_tot ; r_avg_p2 ];
    q_avg=[q_avg_tot ; q_avg_p2 ];
%     alphas_kappa=[alphas_kappa_tot ; alphas_kappa_p2 ];
    alphas_lambda=[alphas_lambda_tot ; alphas_lambda_p2];
    ALL_TRAPPED_POP=[ALL_TRAPPED_POP_tot ; ALL_TRAPPED_POP_p2 ];
    ALL_PASSING_POP=[ALL_PASSING_POP_tot ; ALL_PASSING_POP_p2 ];
    TRAPPED_MINUS_POP=[TRAPPED_MINUS_POP_tot ; TRAPPED_MINUS_POP_p2];
    TRAPPED_PLUS_POP=[TRAPPED_PLUS_POP_tot ; TRAPPED_PLUS_POP_p2 ];
    POTATOES_POP=[POTATOES_POP_tot ; POTATOES_POP_p2 ];
    STAGNATION_POP=[STAGNATION_POP_tot ;STAGNATION_POP_p2 ];
    CO_PASSING_POP=[CO_PASSING_POP_tot ; CO_PASSING_POP_p2 ];
    COUNTER_PASSING_POP=[COUNTER_PASSING_POP_tot ; COUNTER_PASSING_POP_p2 ];
    omega_ik_avg_p3=[omega_ik_avg_tot ; omega_ik_avg_p2 ];
    omega_psi_avg=[omega_psi_avg_tot ; omega_psi_avg_p2 ];
    r_vpll_minus_avg=[r_vpll_minus_avg_tot ; r_vpll_minus_avg_p2 ];
    r_vpll_plus_avg=[r_vpll_plus_avg_tot ; r_vpll_plus_avg_p2 ];
    omega_r_avg=[omega_r_avg_tot ; omega_r_avg_p2 ];
    delta_r_avg=[delta_r_avg_tot ; delta_r_avg_p2 ];
    omega_precess_avg=[omega_precess_avg_tot ; omega_precess_avg_p2 ];
    omega_bounce=[omega_bounce_tot ; omega_bounce_p2 ];
    phidot_avg=[phidot_avg_tot ; phidot_avg_p2 ];
    alphas_ejected=[alphas_ejected_tot ; alphas_ejected_p2 ];
end

ALL_TRAPPED=find(ALL_TRAPPED_POP);
ALL_PASSING=find(ALL_PASSING_POP);
TRAPPED_PLUS=find(TRAPPED_PLUS_POP);
TRAPPED_MINUS=find(TRAPPED_MINUS_POP);
POTATOES=find(POTATOES_POP);
STAGNATION=find(STAGNATION_POP);
CO_PASSING=find(CO_PASSING_POP);
COUNTER_PASSING=find(COUNTER_PASSING_POP);

save(SAVENAME, 'mHe', 'ZHe', 'r_avg', 'q_avg','alphas_lambda', 'ALL_TRAPPED_POP', 'ALL_PASSING_POP', 'TRAPPED_MINUS_POP', 'TRAPPED_PLUS_POP', 'POTATOES_POP', 'STAGNATION_POP', 'CO_PASSING_POP', 'COUNTER_PASSING_POP', 'ALL_TRAPPED', 'ALL_PASSING', 'TRAPPED_PLUS', 'TRAPPED_MINUS', 'POTATOES', 'STAGNATION', 'CO_PASSING', 'COUNTER_PASSING', 'omega_ik_avg','omega_psi_avg', 'r_vpll_minus_avg', 'r_vpll_plus_avg', 'omega_r_avg', 'delta_r_avg', 'omega_precess_avg', 'omega_bounce', 'phidot_avg','alphas_ejected')

