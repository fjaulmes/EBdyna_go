%%
% initial distributions
SAVENAME='initialG_DT_MB_all_pre_collapse.mat'

load('initialG_DT_MB5_pre_collapse.mat')

alphas_Ekin_p1=alphas_Ekin;
alphas_mm_p1=alphas_mm;
alphas_pphi0_p1=alphas_pphi0;
alphas_pos_x_p1=alphas_pos_x;
alphas_pos_z_p1=alphas_pos_z;
alphas_pos_phi_p1=alphas_pos_phi;
alphas_vpll_p1=alphas_vpll;
alphas_psi_p1=alphas_psi;
v_X_p1=v_X;
v_Z_p1=v_Z;
v_phi_p1=v_phi;
% alphas_ejected_p1=alphas_ejected;
% alphas_vE_sq_p1=alphas_vE_sq;
% alphas_psi_value_corr_p1=alphas_psi_value_corr;
% alphas_grad_Phi_p1=alphas_grad_Phi;
% alphas_omega_p1=alphas_omega;
% alphas_Etot_p1=alphas_Etot;
% pos_X_gc_p1=pos_X_gc;
% pos_Z_gc_p1=pos_Z_gc;

load('initialG_DT_MB6_pre_collapse.mat')

alphas_Ekin_p2=alphas_Ekin;
alphas_mm_p2=alphas_mm;
alphas_pphi0_p2=alphas_pphi0;
alphas_pos_x_p2=alphas_pos_x;
alphas_pos_z_p2=alphas_pos_z;
alphas_pos_phi_p2=alphas_pos_phi;
alphas_vpll_p2=alphas_vpll;
alphas_psi_p2=alphas_psi;
v_X_p2=v_X;
v_Z_p2=v_Z;
v_phi_p2=v_phi;
% alphas_ejected_p2=alphas_ejected;
% alphas_vE_sq_p2=alphas_vE_sq;
% alphas_psi_value_corr_p2=alphas_psi_value_corr;
% alphas_grad_Phi_p2=alphas_grad_Phi;
% alphas_omega_p2=alphas_omega;
% alphas_Etot_p2=alphas_Etot;
% pos_X_gc_p2=pos_X_gc;
% pos_Z_gc_p2=pos_Z_gc;




load('initialG_DT_MB7_pre_collapse.mat')

alphas_Ekin_p3=alphas_Ekin;
alphas_mm_p3=alphas_mm;
alphas_pphi0_p3=alphas_pphi0;
alphas_pos_x_p3=alphas_pos_x;
alphas_pos_z_p3=alphas_pos_z;
alphas_pos_phi_p3=alphas_pos_phi;
alphas_vpll_p3=alphas_vpll;
alphas_psi_p3=alphas_psi;
v_X_p3=v_X;
v_Z_p3=v_Z;
v_phi_p3=v_phi;
% alphas_ejected_p3=alphas_ejected;
% alphas_vE_sq_p3=alphas_vE_sq;
% alphas_psi_value_corr_p3=alphas_psi_value_corr;
% alphas_grad_Phi_p3=alphas_grad_Phi;
% alphas_omega_p3=alphas_omega;
% alphas_Etot_p3=alphas_Etot;
% pos_X_gc_p3=pos_X_gc;
% pos_Z_gc_p3=pos_Z_gc;

load('initialG_DT_MB8_pre_collapse.mat')

alphas_Ekin_p4=alphas_Ekin;
alphas_mm_p4=alphas_mm;
alphas_pphi0_p4=alphas_pphi0;
alphas_pos_x_p4=alphas_pos_x;
alphas_pos_z_p4=alphas_pos_z;
alphas_pos_phi_p4=alphas_pos_phi;
alphas_vpll_p4=alphas_vpll;
alphas_psi_p4=alphas_psi;
v_X_p4=v_X;
v_Z_p4=v_Z;
v_phi_p4=v_phi;

alphas_Ekin=[alphas_Ekin_p1 ; alphas_Ekin_p2 ; alphas_Ekin_p3; alphas_Ekin_p4];
alphas_mm=[alphas_mm_p1 ; alphas_mm_p2 ; alphas_mm_p3 ; alphas_mm_p4];
alphas_pphi0=[alphas_pphi0_p1 ; alphas_pphi0_p2 ; alphas_pphi0_p3 ; alphas_pphi0_p4];
alphas_pos_x=[alphas_pos_x_p1 ; alphas_pos_x_p2 ; alphas_pos_x_p3 ; alphas_pos_x_p4];
alphas_pos_z=[alphas_pos_z_p1 ; alphas_pos_z_p2 ; alphas_pos_z_p3 ; alphas_pos_z_p4];
alphas_pos_phi=[alphas_pos_phi_p1 ; alphas_pos_phi_p2 ; alphas_pos_phi_p3 ; alphas_pos_phi_p4];
alphas_vpll=[alphas_vpll_p1 ; alphas_vpll_p2 ; alphas_vpll_p3 ; alphas_vpll_p4];
alphas_psi=[alphas_psi_p1 ; alphas_psi_p2 ; alphas_psi_p3 ; alphas_psi_p4];
v_X=[v_X_p1 ; v_X_p2 ; v_X_p3 ; v_X_p4];
v_Z=[v_Z_p1 ;v_Z_p2 ; v_Z_p3 ; v_Z_p4];
v_phi=[v_phi_p1 ; v_phi_p2 ; v_phi_p3 ; v_phi_p4];
% alphas_ejected=[alphas_ejected_p1 ; alphas_ejected_p2 ; alphas_ejected_p3];
% alphas_vE_sq=[alphas_vE_sq_p1 ; alphas_vE_sq_p2 ; alphas_vE_sq_p3];
% alphas_psi_value_corr=[alphas_psi_value_corr_p1 ; alphas_psi_value_corr_p2 ; alphas_psi_value_corr_p3];
% alphas_grad_Phi=[alphas_grad_Phi_p1 ; alphas_grad_Phi_p2 ; alphas_grad_Phi_p3];
% alphas_omega=[alphas_omega_p1 ; alphas_omega_p2 ; alphas_omega_p3];
% alphas_Etot=[alphas_Etot_p1 ; alphas_Etot_p2 ; alphas_Etot_p3];
% pos_X_gc=[pos_X_gc_p1 ; pos_X_gc_p2 ; pos_X_gc_p3];
% pos_Z_gc=[pos_Z_gc_p1 ; pos_Z_gc_p2 ; pos_Z_gc_p3];

save(SAVENAME,'alphas_Ekin', 'alphas_mm', 'alphas_pphi0', 'alphas_pos_x', 'alphas_pos_z', 'alphas_pos_phi', 'alphas_vpll', 'alphas_psi', 'v_X', 'v_Z', 'v_phi');

%%
% precession stats
SAVENAME='initialG_DT_MB_all_precession_stats.mat'

load('initialG_DT_MB5_precession_stats.mat')

r_avg=r_avg';
 
r_avg_p1=r_avg;
q_avg_p1=q_avg;
alphas_kappa_p1=alphas_kappa;
alphas_lambda_p1=alphas_lambda;
ALL_TRAPPED_POP_p1=ALL_TRAPPED_POP;
ALL_PASSING_POP_p1=ALL_PASSING_POP;
TRAPPED_MINUS_POP_p1=TRAPPED_MINUS_POP;
TRAPPED_PLUS_POP_p1=TRAPPED_PLUS_POP;
POTATOES_POP_p1=POTATOES_POP;
STAGNATION_POP_p1=STAGNATION_POP;
CO_PASSING_POP_p1=CO_PASSING_POP;
COUNTER_PASSING_POP_p1=COUNTER_PASSING_POP;
omega_ik_avg_p1=omega_ik_avg;
omega_psi_avg_p1=omega_psi_avg;
r_vpll_minus_avg_p1=r_vpll_minus_avg;
r_vpll_plus_avg_p1=r_vpll_plus_avg;
omega_r_avg_p1=omega_r_avg;
delta_r_avg_p1=delta_r_avg;
omega_precess_avg_p1=omega_precess_avg;
omega_bounce_p1=omega_bounce;
phidot_avg_p1=phidot_avg;

load('initialG_DT_MB6_precession_stats.mat')

r_avg=r_avg';
 
r_avg_p2=r_avg;
q_avg_p2=q_avg;
alphas_kappa_p2=alphas_kappa;
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

load('initialG_DT_MB7_precession_stats.mat')

r_avg=r_avg';
 
r_avg_p3=r_avg;
q_avg_p3=q_avg;
alphas_kappa_p3=alphas_kappa;
alphas_lambda_p3=alphas_lambda;
ALL_TRAPPED_POP_p3=ALL_TRAPPED_POP;
ALL_PASSING_POP_p3=ALL_PASSING_POP;
TRAPPED_MINUS_POP_p3=TRAPPED_MINUS_POP;
TRAPPED_PLUS_POP_p3=TRAPPED_PLUS_POP;
POTATOES_POP_p3=POTATOES_POP;
STAGNATION_POP_p3=STAGNATION_POP;
CO_PASSING_POP_p3=CO_PASSING_POP;
COUNTER_PASSING_POP_p3=COUNTER_PASSING_POP;
omega_ik_avg_p3=omega_ik_avg;
omega_psi_avg_p3=omega_psi_avg;
r_vpll_minus_avg_p3=r_vpll_minus_avg;
r_vpll_plus_avg_p3=r_vpll_plus_avg;
omega_r_avg_p3=omega_r_avg;
delta_r_avg_p3=delta_r_avg;
omega_precess_avg_p3=omega_precess_avg;
omega_bounce_p3=omega_bounce;
phidot_avg_p3=phidot_avg;

load('initialG_DT_MB8_precession_stats.mat')

r_avg=r_avg';
 
r_avg_p4=r_avg;
q_avg_p4=q_avg;
alphas_kappa_p4=alphas_kappa;
alphas_lambda_p4=alphas_lambda;
ALL_TRAPPED_POP_p4=ALL_TRAPPED_POP;
ALL_PASSING_POP_p4=ALL_PASSING_POP;
TRAPPED_MINUS_POP_p4=TRAPPED_MINUS_POP;
TRAPPED_PLUS_POP_p4=TRAPPED_PLUS_POP;
POTATOES_POP_p4=POTATOES_POP;
STAGNATION_POP_p4=STAGNATION_POP;
CO_PASSING_POP_p4=CO_PASSING_POP;
COUNTER_PASSING_POP_p4=COUNTER_PASSING_POP;
omega_ik_avg_p4=omega_ik_avg;
omega_psi_avg_p4=omega_psi_avg;
r_vpll_minus_avg_p4=r_vpll_minus_avg;
r_vpll_plus_avg_p4=r_vpll_plus_avg;
omega_r_avg_p4=omega_r_avg;
delta_r_avg_p4=delta_r_avg;
omega_precess_avg_p4=omega_precess_avg;
omega_bounce_p4=omega_bounce;
phidot_avg_p4=phidot_avg;

r_avg=[r_avg_p1 ; r_avg_p2 ; r_avg_p3 ; r_avg_p4];
q_avg=[q_avg_p1 ; q_avg_p2 ; q_avg_p3 ; q_avg_p4];
alphas_kappa=[alphas_kappa_p1 ; alphas_kappa_p2 ; alphas_kappa_p3 ; alphas_kappa_p4];
alphas_lambda=[alphas_lambda_p1 ; alphas_lambda_p2 ; alphas_lambda_p3 ; alphas_lambda_p4];
ALL_TRAPPED_POP=[ALL_TRAPPED_POP_p1 ; ALL_TRAPPED_POP_p2 ; ALL_TRAPPED_POP_p3 ; ALL_TRAPPED_POP_p4];
ALL_PASSING_POP=[ALL_PASSING_POP_p1 ; ALL_PASSING_POP_p2 ; ALL_PASSING_POP_p3 ; ALL_PASSING_POP_p4];
TRAPPED_MINUS_POP=[TRAPPED_MINUS_POP_p1 ; TRAPPED_MINUS_POP_p2 ; TRAPPED_MINUS_POP_p3 ; TRAPPED_MINUS_POP_p4];
TRAPPED_PLUS_POP=[TRAPPED_PLUS_POP_p1 ; TRAPPED_PLUS_POP_p2 ; TRAPPED_PLUS_POP_p3 ; TRAPPED_PLUS_POP_p4];
POTATOES_POP=[POTATOES_POP_p1 ; POTATOES_POP_p2 ; POTATOES_POP_p3 ; POTATOES_POP_p4];
STAGNATION_POP=[STAGNATION_POP_p1 ;STAGNATION_POP_p2 ; STAGNATION_POP_p3 ; STAGNATION_POP_p4];
CO_PASSING_POP=[CO_PASSING_POP_p1 ; CO_PASSING_POP_p2 ; CO_PASSING_POP_p3 ; CO_PASSING_POP_p4];
COUNTER_PASSING_POP=[COUNTER_PASSING_POP_p1 ; COUNTER_PASSING_POP_p2 ; COUNTER_PASSING_POP_p3 ; COUNTER_PASSING_POP_p4];
omega_ik_avg_p3=[omega_ik_avg_p1 ; omega_ik_avg_p2 ; omega_ik_avg_p3 ; omega_ik_avg_p4];
omega_psi_avg=[omega_psi_avg_p1 ; omega_psi_avg_p2 ; omega_psi_avg_p3 ; omega_psi_avg_p4];
r_vpll_minus_avg=[r_vpll_minus_avg_p1 ; r_vpll_minus_avg_p2 ; r_vpll_minus_avg_p3 ; r_vpll_minus_avg_p4];
r_vpll_plus_avg=[r_vpll_plus_avg_p1 ; r_vpll_plus_avg_p2 ; r_vpll_plus_avg_p3 ; r_vpll_plus_avg_p4];
omega_r_avg=[omega_r_avg_p1 ; omega_r_avg_p2 ; omega_r_avg_p3 ; omega_r_avg_p4];
delta_r_avg=[delta_r_avg_p1 ; delta_r_avg_p2 ; delta_r_avg_p3 ; delta_r_avg_p4];
omega_precess_avg=[omega_precess_avg_p1 ; omega_precess_avg_p2 ; omega_precess_avg_p3 ; omega_precess_avg_p4];
omega_bounce=[omega_bounce_p1 ; omega_bounce_p2 ; omega_bounce_p3 ; omega_bounce_p4];
phidot_avg=[phidot_avg_p1 ; phidot_avg_p2 ; phidot_avg_p3 ; phidot_avg_p4];


ALL_TRAPPED=find(ALL_TRAPPED_POP);
ALL_PASSING=find(ALL_PASSING_POP);
TRAPPED_PLUS=find(TRAPPED_PLUS_POP);
TRAPPED_MINUS=find(TRAPPED_MINUS_POP);
POTATOES=find(POTATOES_POP);
STAGNATION=find(STAGNATION_POP);
CO_PASSING=find(CO_PASSING_POP);
COUNTER_PASSING=find(COUNTER_PASSING_POP);

save(SAVENAME, 'mHe', 'ZHe', 'r_avg', 'q_avg', 'alphas_kappa', 'alphas_lambda', 'ALL_TRAPPED_POP', 'ALL_PASSING_POP', 'TRAPPED_MINUS_POP', 'TRAPPED_PLUS_POP', 'POTATOES_POP', 'STAGNATION_POP', 'CO_PASSING_POP', 'COUNTER_PASSING_POP', 'ALL_TRAPPED', 'ALL_PASSING', 'TRAPPED_PLUS', 'TRAPPED_MINUS', 'POTATOES', 'STAGNATION', 'CO_PASSING', 'COUNTER_PASSING', 'omega_ik_avg','omega_psi_avg', 'r_vpll_minus_avg', 'r_vpll_plus_avg', 'omega_r_avg', 'delta_r_avg', 'omega_precess_avg', 'omega_bounce', 'phidot_avg', 'omega_crash')

%%
% post simulation results
SAVENAME='DT_MBall_collapse_Glisa_fc0p8h2_G160414.mat'

load('DT_MB5_collapse_Glisa_fc0p8h2_G160414.mat');
alphas_Ekin_p1=alphas_Ekin;
alphas_mm_p1=alphas_mm;
alphas_pphi0_p1=alphas_pphi0;
alphas_pos_x_p1=alphas_pos_x;
alphas_pos_z_p1=alphas_pos_z;
alphas_pos_phi_p1=alphas_pos_phi;
alphas_vpll_p1=alphas_vpll;
alphas_psi_p1=alphas_psi;
v_X_p1=v_X;
v_Z_p1=v_Z;
v_phi_p1=v_phi;
alphas_ejected_p1=alphas_ejected;
alphas_vE_sq_p1=alphas_vE_sq;
alphas_psi_value_corr_p1=alphas_psi_value_corr;
alphas_grad_Phi_p1=alphas_grad_Phi;
alphas_omega_p1=alphas_omega;
alphas_Etot_p1=alphas_Etot;
pos_X_gc_p1=pos_X_gc;
pos_Z_gc_p1=pos_Z_gc;


load('DT_MB6_collapse_Glisa_fc0p8h2_G160414.mat')
alphas_Ekin_p2=alphas_Ekin;
alphas_mm_p2=alphas_mm;
alphas_pphi0_p2=alphas_pphi0;
alphas_pos_x_p2=alphas_pos_x;
alphas_pos_z_p2=alphas_pos_z;
alphas_pos_phi_p2=alphas_pos_phi;
alphas_vpll_p2=alphas_vpll;
alphas_psi_p2=alphas_psi;
v_X_p2=v_X;
v_Z_p2=v_Z;
v_phi_p2=v_phi;
alphas_ejected_p2=alphas_ejected;
alphas_vE_sq_p2=alphas_vE_sq;
alphas_psi_value_corr_p2=alphas_psi_value_corr;
alphas_grad_Phi_p2=alphas_grad_Phi;
alphas_omega_p2=alphas_omega;
alphas_Etot_p2=alphas_Etot;
pos_X_gc_p2=pos_X_gc;
pos_Z_gc_p2=pos_Z_gc;


load('DT_MB7_collapse_Glisa_fc0p8h2_G160414.mat')

alphas_Ekin_p3=alphas_Ekin;
alphas_mm_p3=alphas_mm;
alphas_pphi0_p3=alphas_pphi0;
alphas_pos_x_p3=alphas_pos_x;
alphas_pos_z_p3=alphas_pos_z;
alphas_pos_phi_p3=alphas_pos_phi;
alphas_vpll_p3=alphas_vpll;
alphas_psi_p3=alphas_psi;
v_X_p3=v_X;
v_Z_p3=v_Z;
v_phi_p3=v_phi;
alphas_ejected_p3=alphas_ejected;
alphas_vE_sq_p3=alphas_vE_sq;
alphas_psi_value_corr_p3=alphas_psi_value_corr;
alphas_grad_Phi_p3=alphas_grad_Phi;
alphas_omega_p3=alphas_omega;
alphas_Etot_p3=alphas_Etot;
pos_X_gc_p3=pos_X_gc;
pos_Z_gc_p3=pos_Z_gc;

load('DT_MB8_collapse_Glisa_fc0p8h2_G160414.mat')

alphas_Ekin_p4=alphas_Ekin;
alphas_mm_p4=alphas_mm;
alphas_pphi0_p4=alphas_pphi0;
alphas_pos_x_p4=alphas_pos_x;
alphas_pos_z_p4=alphas_pos_z;
alphas_pos_phi_p4=alphas_pos_phi;
alphas_vpll_p4=alphas_vpll;
alphas_psi_p4=alphas_psi;
v_X_p4=v_X;
v_Z_p4=v_Z;
v_phi_p4=v_phi;
alphas_ejected_p4=alphas_ejected;
alphas_vE_sq_p4=alphas_vE_sq;
alphas_psi_value_corr_p4=alphas_psi_value_corr;
alphas_grad_Phi_p4=alphas_grad_Phi;
alphas_omega_p4=alphas_omega;
alphas_Etot_p4=alphas_Etot;
pos_X_gc_p4=pos_X_gc;
pos_Z_gc_p4=pos_Z_gc;

alphas_Ekin=[alphas_Ekin_p1 ; alphas_Ekin_p2 ; alphas_Ekin_p3 ; alphas_Ekin_p4];
alphas_mm=[alphas_mm_p1 ; alphas_mm_p2 ; alphas_mm_p3 ; alphas_mm_p4];
alphas_pphi0=[alphas_pphi0_p1 ; alphas_pphi0_p2 ; alphas_pphi0_p3 ; alphas_pphi0_p4];
alphas_pos_x=[alphas_pos_x_p1 ; alphas_pos_x_p2 ; alphas_pos_x_p3 ; alphas_pos_x_p4];
alphas_pos_z=[alphas_pos_z_p1 ; alphas_pos_z_p2 ; alphas_pos_z_p3 ; alphas_pos_z_p4];
alphas_pos_phi=[alphas_pos_phi_p1 ; alphas_pos_phi_p2 ; alphas_pos_phi_p3 ; alphas_pos_phi_p4];
alphas_vpll=[alphas_vpll_p1 ; alphas_vpll_p2 ; alphas_vpll_p3 ; alphas_vpll_p4];
alphas_psi=[alphas_psi_p1 ; alphas_psi_p2 ; alphas_psi_p3 ; alphas_psi_p4];
v_X=[v_X_p1 ; v_X_p2 ; v_X_p3 ; v_X_p4];
v_Z=[v_Z_p1 ;v_Z_p2 ; v_Z_p3 ; v_Z_p4];
v_phi=[v_phi_p1 ; v_phi_p2 ; v_phi_p3 ; v_phi_p4];
alphas_ejected=[alphas_ejected_p1 ; alphas_ejected_p2 ; alphas_ejected_p3; alphas_ejected_p4];
alphas_vE_sq=[alphas_vE_sq_p1 ; alphas_vE_sq_p2 ; alphas_vE_sq_p3 ; alphas_vE_sq_p4];
alphas_psi_value_corr=[alphas_psi_value_corr_p1 ; alphas_psi_value_corr_p2 ; alphas_psi_value_corr_p3 ; alphas_psi_value_corr_p4];
alphas_grad_Phi=[alphas_grad_Phi_p1 ; alphas_grad_Phi_p2 ; alphas_grad_Phi_p3 ; alphas_grad_Phi_p4];
alphas_omega=[alphas_omega_p1 ; alphas_omega_p2 ; alphas_omega_p3 ; alphas_omega_p4];
alphas_Etot=[alphas_Etot_p1 ; alphas_Etot_p2 ; alphas_Etot_p3 ; alphas_Etot_p4];
pos_X_gc=[pos_X_gc_p1 ; pos_X_gc_p2 ; pos_X_gc_p3 ; pos_X_gc_p4];
pos_Z_gc=[pos_Z_gc_p1 ; pos_Z_gc_p2 ; pos_Z_gc_p3 ; pos_Z_gc_p4];

save(SAVENAME,'alphas_Ekin', 'alphas_mm', 'alphas_pphi0', 'alphas_pos_x', 'alphas_pos_z', 'alphas_pos_phi', 'alphas_vpll', 'alphas_psi', 'v_X', 'v_Z', 'v_phi', 'alphas_ejected', 'alphas_vE_sq', 'alphas_psi_value_corr', 'alphas_grad_Phi', 'time', 'frame_rank_precise', 'alphas_omega', 'alphas_Etot', 'pos_X_gc', 'pos_Z_gc');

disp('done....')

