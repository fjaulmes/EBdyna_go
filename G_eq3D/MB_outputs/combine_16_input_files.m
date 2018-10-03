reset_data_analysis_environment
format compact
%%
NB_PROCESS=16
% initial distributions
SAVENAME='./initial_MBD_100_pre_collapse_all.mat'

n=1;
%INPUTNAME=strcat('initial_MBD_100_TAE_distribution',num2str(n),'.mat')
%load(INPUTNAME,'alphas_weight');
INPUTNAME=strcat('initial_MBD_100_pre_collapse',num2str(n),'.mat')
load(INPUTNAME)



for n=2:NB_PROCESS
    alphas_weight_tot=[alphas_weight];
    alphas_Ekin_tot=[alphas_Ekin ];
    alphas_mm_tot=[alphas_mm  ];
    alphas_pphi0_tot=[alphas_pphi0  ];
    alphas_pos_x_tot=[alphas_pos_x];
    alphas_pos_z_tot=[alphas_pos_z  ];
    alphas_pos_phi_tot=[alphas_pos_phi  ];
    alphas_vpll_tot=[alphas_vpll ];
    alphas_psi_tot=[alphas_psi ];
    v_X_tot=[v_X  ];
    v_Z_tot=[v_Z];
    v_phi_tot=[v_phi ];

    %INPUTNAME=strcat('initial_MBD_100_TAE_distribution',num2str(n),'.mat')
    %load(INPUTNAME,'alphas_weight');
    INPUTNAME=strcat('initial_MBD_100_pre_collapse',num2str(n),'.mat')
    load(INPUTNAME)
	
    alphas_weight_p2=alphas_weight;
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
    
    
    alphas_weight=[alphas_weight_tot ; alphas_weight_p2 ];
    alphas_Ekin=[alphas_Ekin_tot ; alphas_Ekin_p2 ];
    alphas_mm=[alphas_mm_tot ; alphas_mm_p2 ];
    alphas_pphi0=[alphas_pphi0_tot ; alphas_pphi0_p2 ];
    alphas_pos_x=[alphas_pos_x_tot ; alphas_pos_x_p2 ];
    alphas_pos_z=[alphas_pos_z_tot ; alphas_pos_z_p2 ];
    alphas_pos_phi=[alphas_pos_phi_tot ; alphas_pos_phi_p2 ];
    alphas_vpll=[alphas_vpll_tot ; alphas_vpll_p2 ];
    alphas_psi=[alphas_psi_tot ; alphas_psi_p2 ];
    v_X=[v_X_tot ; v_X_p2 ];
    v_Z=[v_Z_tot ;v_Z_p2 ];
    v_phi=[v_phi_tot ; v_phi_p2 ];
    % alphas_ejected=[alphas_ejected_tot ; alphas_ejected_p2 ; alphas_ejected_p3];
    % alphas_vE_sq=[alphas_vE_sq_p1 ; alphas_vE_sq_p2 ; alphas_vE_sq_p3];
    % alphas_psi_value_corr=[alphas_psi_value_corr_p1 ; alphas_psi_value_corr_p2 ; alphas_psi_value_corr_p3];
    % alphas_grad_Phi=[alphas_grad_Phi_p1 ; alphas_grad_Phi_p2 ; alphas_grad_Phi_p3];
    % alphas_omega=[alphas_omega_p1 ; alphas_omega_p2 ; alphas_omega_p3];
    % alphas_Etot=[alphas_Etot_p1 ; alphas_Etot_p2 ; alphas_Etot_p3];
    % pos_X_gc=[pos_X_gc_p1 ; pos_X_gc_p2 ; pos_X_gc_p3];
    % pos_Z_gc=[pos_Z_gc_p1 ; pos_Z_gc_p2 ; pos_Z_gc_p3];
end
    bX=interp2(scale_X,scale_Z,bX_XZ_map',alphas_pos_x,alphas_pos_z);
    bZ=interp2(scale_X,scale_Z,bZ_XZ_map',alphas_pos_x,alphas_pos_z);
    alphas_Bfield=interp2(scale_X,scale_Z,Btot_XZ_map',alphas_pos_x,alphas_pos_z);
    bphi=sqrt(1-bX.^2-bZ.^2);
    pos_X_gc=alphas_pos_x+(mHe/eV)*(1/ZHe)*(v_Z.*bphi-v_phi.*bZ)./alphas_Bfield;
    pos_Z_gc=alphas_pos_z+(mHe/eV)*(1/ZHe)*(v_phi.*bX-v_X.*bphi)./alphas_Bfield;

Nalphas_simulated=length(alphas_vpll)


save(SAVENAME,'mHe', 'ZHe','Nalphas_simulated','alphas_weight','alphas_Ekin','alphas_mm', 'alphas_pphi0', 'pos_X_gc','pos_Z_gc','alphas_pos_x', 'alphas_pos_z', 'alphas_pos_phi', 'alphas_vpll', 'alphas_psi', 'v_X', 'v_Z', 'v_phi');




%%
% precession stats
SAVENAME='./initial_MBD_100_precession_stats_all.mat'

n=1;
INPUTNAME=strcat('initial_MBD_100_precession_stats',num2str(n),'.mat')
load(INPUTNAME)
% r_avg=r_avg';
vpll_avg=vpll_avg';

for n=2:NB_PROCESS
    vpll_avg_tot=vpll_avg;
    r_avg_tot=r_avg;
    q_avg_tot=q_avg;
    alphas_lambda_tot=alphas_lambda;
    ALL_TRAPPED_POP_tot=ALL_TRAPPED_POP;
    ALL_PASSING_POP_tot=ALL_PASSING_POP;
    TRAPPED_MINUS_POP_tot=TRAPPED_MINUS_POP;
    TRAPPED_PLUS_POP_tot=TRAPPED_PLUS_POP;
    POTATOES_POP_tot=POTATOES_POP;
    STAGNATION_POP_tot=STAGNATION_POP;
    CO_PASSING_POP_tot=CO_PASSING_POP;
    COUNTER_PASSING_POP_tot=COUNTER_PASSING_POP;
    omega_ae_avg_tot=omega_ae_avg;
    omega_psi_avg_tot=omega_psi_avg;
    r_vpll_minus_avg_tot=r_vpll_minus_avg;
    r_vpll_plus_avg_tot=r_vpll_plus_avg;
    omega_r_avg_tot=omega_r_avg;
    delta_r_avg_tot=delta_r_avg;
    omega_precess_avg_tot=omega_precess_avg;
    omega_bounce_tot=omega_bounce;
    phidot_avg_tot=phidot_avg;
    
    INPUTNAME=strcat('initial_MBD_100_precession_stats',num2str(n),'.mat')
    load(INPUTNAME)
    vpll_avg=vpll_avg';   
        
%     r_avg=r_avg';
    
    vpll_avg_p2=vpll_avg;
    r_avg_p2=r_avg;
    q_avg_p2=q_avg;
    alphas_lambda_p2=alphas_lambda;
    ALL_TRAPPED_POP_p2=ALL_TRAPPED_POP;
    ALL_PASSING_POP_p2=ALL_PASSING_POP;
    TRAPPED_MINUS_POP_p2=TRAPPED_MINUS_POP;
    TRAPPED_PLUS_POP_p2=TRAPPED_PLUS_POP;
    POTATOES_POP_p2=POTATOES_POP;
    STAGNATION_POP_p2=STAGNATION_POP;
    CO_PASSING_POP_p2=CO_PASSING_POP;
    COUNTER_PASSING_POP_p2=COUNTER_PASSING_POP;
    omega_ae_avg_p2=omega_ae_avg;
    omega_psi_avg_p2=omega_psi_avg;
    r_vpll_minus_avg_p2=r_vpll_minus_avg;
    r_vpll_plus_avg_p2=r_vpll_plus_avg;
    omega_r_avg_p2=omega_r_avg;
    delta_r_avg_p2=delta_r_avg;
    omega_precess_avg_p2=omega_precess_avg;
    omega_bounce_p2=omega_bounce;
    phidot_avg_p2=phidot_avg;
    
    
    
    
    vpll_avg=[vpll_avg_tot ; vpll_avg_p2 ];
    r_avg=[r_avg_tot ; r_avg_p2 ];
    q_avg=[q_avg_tot ; q_avg_p2 ];
    alphas_lambda=[alphas_lambda_tot ; alphas_lambda_p2];
    ALL_TRAPPED_POP=[ALL_TRAPPED_POP_tot ; ALL_TRAPPED_POP_p2 ];
    ALL_PASSING_POP=[ALL_PASSING_POP_tot ; ALL_PASSING_POP_p2 ];
    TRAPPED_MINUS_POP=[TRAPPED_MINUS_POP_tot ; TRAPPED_MINUS_POP_p2];
    TRAPPED_PLUS_POP=[TRAPPED_PLUS_POP_tot ; TRAPPED_PLUS_POP_p2 ];
    POTATOES_POP=[POTATOES_POP_tot ; POTATOES_POP_p2 ];
    STAGNATION_POP=[STAGNATION_POP_tot ;STAGNATION_POP_p2 ];
    CO_PASSING_POP=[CO_PASSING_POP_tot ; CO_PASSING_POP_p2 ];
    COUNTER_PASSING_POP=[COUNTER_PASSING_POP_tot ; COUNTER_PASSING_POP_p2 ];
    omega_ae_avg=[omega_ae_avg_tot ; omega_ae_avg_p2 ];
    omega_psi_avg=[omega_psi_avg_tot ; omega_psi_avg_p2 ];
    r_vpll_minus_avg=[r_vpll_minus_avg_tot ; r_vpll_minus_avg_p2 ];
    r_vpll_plus_avg=[r_vpll_plus_avg_tot ; r_vpll_plus_avg_p2 ];
    omega_r_avg=[omega_r_avg_tot ; omega_r_avg_p2 ];
    delta_r_avg=[delta_r_avg_tot ; delta_r_avg_p2 ];
    omega_precess_avg=[omega_precess_avg_tot ; omega_precess_avg_p2 ];
    omega_bounce=[omega_bounce_tot ; omega_bounce_p2 ];
    phidot_avg=[phidot_avg_tot ; phidot_avg_p2 ];
    
end

ALL_TRAPPED=find(ALL_TRAPPED_POP);
ALL_PASSING=find(ALL_PASSING_POP);
TRAPPED_PLUS=find(TRAPPED_PLUS_POP);
TRAPPED_MINUS=find(TRAPPED_MINUS_POP);
POTATOES=find(POTATOES_POP);
STAGNATION=find(STAGNATION_POP);
CO_PASSING=find(CO_PASSING_POP);
COUNTER_PASSING=find(COUNTER_PASSING_POP);

save(SAVENAME, 'mHe', 'ZHe', 'vpll_avg', 'r_avg', 'q_avg',  'alphas_lambda', 'ALL_TRAPPED_POP', 'ALL_PASSING_POP', 'TRAPPED_MINUS_POP', 'TRAPPED_PLUS_POP', 'POTATOES_POP', 'STAGNATION_POP', 'CO_PASSING_POP', 'COUNTER_PASSING_POP', 'ALL_TRAPPED', 'ALL_PASSING', 'TRAPPED_PLUS', 'TRAPPED_MINUS', 'POTATOES', 'STAGNATION', 'CO_PASSING', 'COUNTER_PASSING', 'omega_ae_avg','omega_psi_avg', 'r_vpll_minus_avg', 'r_vpll_plus_avg', 'omega_r_avg', 'delta_r_avg', 'omega_precess_avg', 'omega_bounce', 'phidot_avg')
