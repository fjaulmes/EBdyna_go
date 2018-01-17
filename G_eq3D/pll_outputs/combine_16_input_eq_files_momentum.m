%%
NB_PROCESS=16
load('initial_alphas_MB_D_distribution1.mat')

% initial distributions
SAVENAME='./initial_MB_D_pre_collapse_all.mat'

n=1;
INPUTNAME=strcat('initial_MB_D_pre_collapse',num2str(n),'.mat')
load(INPUTNAME)
load(strcat('initial_MB_D_precession',num2str(n),'.mat'),'alphas_ejected');
load(strcat('initial_MB_D_precession_stats',num2str(n),'.mat'),'r_avg');
alphas_momentum=alphas_momentum(~alphas_ejected);
alphas_momentum_interp=interp1(radial_r_value_flux,rotation_profile,r_avg);
alphas_momentum=alphas_momentum_interp;

for n=2:NB_PROCESS
    alphas_Ekin_tot=[alphas_Ekin ];
    alphas_mm_tot=[alphas_mm  ];
    alphas_pphi0_tot=[alphas_pphi0  ];
    alphas_pos_x_tot=[alphas_pos_x];
    alphas_pos_z_tot=[alphas_pos_z  ];
    alphas_pos_phi_tot=[alphas_pos_phi  ];
    alphas_vpll_tot=[alphas_vpll ];
    alphas_psi_tot=[alphas_psi ];
    alphas_momentum_tot=[alphas_momentum ];
    v_X_tot=[v_X  ];
    v_Z_tot=[v_Z];
    v_phi_tot=[v_phi ];

    INPUTNAME=strcat('initial_MB_D_pre_collapse',num2str(n),'.mat')
    load(INPUTNAME)
    load(strcat('initial_MB_D_precession',num2str(n),'.mat'),'alphas_ejected');
    load(strcat('initial_MB_D_precession_stats',num2str(n),'.mat'),'r_avg');
    alphas_momentum=alphas_momentum(~alphas_ejected);
    alphas_momentum_interp=interp1(radial_r_value_flux,rotation_profile,r_avg);
    alphas_momentum=alphas_momentum_interp;
    
    alphas_Ekin_p2=alphas_Ekin;
    alphas_mm_p2=alphas_mm;
    alphas_pphi0_p2=alphas_pphi0;
    alphas_pos_x_p2=alphas_pos_x;
    alphas_pos_z_p2=alphas_pos_z;
    alphas_pos_phi_p2=alphas_pos_phi;
    alphas_vpll_p2=alphas_vpll;
    alphas_psi_p2=alphas_psi;
    alphas_momentum_p2=alphas_momentum;
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
    
    
    alphas_Ekin=[alphas_Ekin_tot ; alphas_Ekin_p2 ];
    alphas_mm=[alphas_mm_tot ; alphas_mm_p2 ];
    alphas_pphi0=[alphas_pphi0_tot ; alphas_pphi0_p2 ];
    alphas_pos_x=[alphas_pos_x_tot ; alphas_pos_x_p2 ];
    alphas_pos_z=[alphas_pos_z_tot ; alphas_pos_z_p2 ];
    alphas_pos_phi=[alphas_pos_phi_tot ; alphas_pos_phi_p2 ];
    alphas_vpll=[alphas_vpll_tot ; alphas_vpll_p2 ];
    alphas_psi=[alphas_psi_tot ; alphas_psi_p2 ];
    alphas_momentum=[alphas_momentum_tot ; alphas_momentum_p2 ];
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

Nalphas_simulated=length(alphas_vpll)


save(SAVENAME,'Nalphas_simulated','alphas_Ekin', 'alphas_momentum','alphas_mm', 'alphas_pphi0', 'alphas_pos_x', 'alphas_pos_z', 'alphas_pos_phi', 'alphas_vpll', 'alphas_psi', 'v_X', 'v_Z', 'v_phi');

    alphas_Ekin_tot=[alphas_Ekin ];
    alphas_mm_tot=[alphas_mm  ];
    alphas_pphi0_tot=[alphas_pphi0  ];
    alphas_pos_x_tot=[alphas_pos_x];
    alphas_pos_z_tot=[alphas_pos_z  ];
    alphas_pos_phi_tot=[alphas_pos_phi  ];
    alphas_vpll_tot=[alphas_vpll ];
    alphas_psi_tot=[alphas_psi ];
    alphas_momentum_tot=[alphas_momentum ];
    v_X_tot=[v_X  ];
    v_Z_tot=[v_Z];
    v_phi_tot=[v_phi ];

%%
% precession stats
SAVENAME='./initial_MB_D_precession_stats_all.mat'

n=1;
INPUTNAME=strcat('initial_MB_D_precession_stats',num2str(n),'.mat')
load(INPUTNAME)
% r_avg=r_avg';

for n=2:NB_PROCESS
    r_avg_tot=r_avg;
    q_avg_tot=q_avg;
    alphas_kappa_tot=alphas_kappa;
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
    
    INPUTNAME=strcat('initial_MB_D_precession_stats',num2str(n),'.mat')
    load(INPUTNAME)
    
        
%     r_avg=r_avg';
    
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
    
    
    
    
    r_avg=[r_avg_tot ; r_avg_p2 ];
    q_avg=[q_avg_tot ; q_avg_p2 ];
    alphas_kappa=[alphas_kappa_tot ; alphas_kappa_p2 ];
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
    
end

ALL_TRAPPED=find(ALL_TRAPPED_POP);
ALL_PASSING=find(ALL_PASSING_POP);
TRAPPED_PLUS=find(TRAPPED_PLUS_POP);
TRAPPED_MINUS=find(TRAPPED_MINUS_POP);
POTATOES=find(POTATOES_POP);
STAGNATION=find(STAGNATION_POP);
CO_PASSING=find(CO_PASSING_POP);
COUNTER_PASSING=find(COUNTER_PASSING_POP);

save(SAVENAME, 'mHe', 'ZHe', 'r_avg', 'q_avg', 'alphas_kappa', 'alphas_lambda', 'ALL_TRAPPED_POP', 'ALL_PASSING_POP', 'TRAPPED_MINUS_POP', 'TRAPPED_PLUS_POP', 'POTATOES_POP', 'STAGNATION_POP', 'CO_PASSING_POP', 'COUNTER_PASSING_POP', 'ALL_TRAPPED', 'ALL_PASSING', 'TRAPPED_PLUS', 'TRAPPED_MINUS', 'POTATOES', 'STAGNATION', 'CO_PASSING', 'COUNTER_PASSING', 'omega_ik_avg','omega_psi_avg', 'r_vpll_minus_avg', 'r_vpll_plus_avg', 'omega_r_avg', 'delta_r_avg', 'omega_precess_avg', 'omega_bounce', 'phidot_avg', 'omega_crash')


disp('done....')



%%

b_phi=interp2(scale_X,scale_Z,bphi_XZ_map',alphas_pos_x,alphas_pos_z);
alphas_Bfield=interp2(scale_X,scale_Z,Btot_XZ_map',alphas_pos_x,alphas_pos_z);

%radial bins for toroidal momentum
psipos_avg=interp1(radial_r_value_flux,1:Nradial,r_avg);
radial_bin_size=16;
radial_bin_half_size=0.5*(radial_bin_size);

radial_bins=[radial_bin_half_size+1:radial_bin_size:simulation_size_r+3*radial_bin_half_size+1]';
N_radial_bins=size(radial_bins,1)
radial_bins_lims=[1:radial_bin_size:radial_bin_size*N_radial_bins+1]';

psipos_profile=zeros(N_radial_bins,1);
phidot_profile=zeros(N_radial_bins,1);
trapped_fraction_profile=zeros(N_radial_bins,1);

phidot_global=mean(phidot_avg(ALL_PASSING_POP|STAGNATION_POP));

for n=1:N_radial_bins
    BINPOP=find((psipos_avg>=radial_bins_lims(n)).*(psipos_avg<radial_bins_lims(n+1)).*(ALL_PASSING_POP|STAGNATION_POP));
    psipos_profile(n)=mean(psipos_avg(BINPOP));
    phidot_profile(n)=mean(phidot_avg(BINPOP));
    BINPOPTOT=find((psipos_avg>=radial_bins_lims(n)).*(psipos_avg<radial_bins_lims(n+1)));
    BINPOPTR=find((psipos_avg>=radial_bins_lims(n)).*(psipos_avg<radial_bins_lims(n+1)).*ALL_TRAPPED_POP);
    trapped_fraction_profile(n)=length(BINPOPTR)/length(BINPOPTOT);
end
v_pll_tot=alphas_vpll+((alphas_pos_x_tot+R0).*(alphas_momentum_tot/R0-phidot_global))./b_phi;
v_phi_tot=v_pll_tot.*b_phi;
vtot2=v_X_tot.^2+v_Z_tot.^2+v_phi_tot.^2;
alphas_Ekin_tot=0.5*(mHe/eV)*vtot2;
alphas_vpll_tot=v_pll_tot;
alphas_mm_tot=(alphas_Ekin_tot-0.5*((mHe/eV)*alphas_vpll_tot.^2))./alphas_Bfield;

[psipos_dist_profile hbins]=histc(psipos_avg,radial_bins_lims);

psipos_profile_mom=zeros(N_radial_bins,1);
phidot_profile_mom=zeros(N_radial_bins,1);

for n=1:N_radial_bins
    BINPOP=find((psipos_avg>=radial_bins_lims(n)).*(psipos_avg<radial_bins_lims(n+1)).*(ALL_PASSING_POP|STAGNATION_POP));
    psipos_profile_mom(n)=mean(psipos_avg(BINPOP));
    phidot_profile_mom(n)=mean(v_phi_tot(BINPOP)./(alphas_pos_x_tot(BINPOP)+R0));
end

%%
% new initial distributions

init_rank=1;

for n=1:NB_PROCESS
    INPUTNAME=strcat('initial_alphas_MB_D_distribution',num2str(n),'.mat');
    load(INPUTNAME);
    load(strcat('initial_MB_D_pre_collapse',num2str(n),'.mat'),'Nalphas_simulated')
    end_rank=init_rank+Nalphas_simulated-1;
    alphas_pos_x=alphas_pos_x_tot(init_rank:end_rank);
    alphas_pos_z=alphas_pos_z_tot(init_rank:end_rank);
    alphas_pos_phi=alphas_pos_phi_tot(init_rank:end_rank);
    alphas_Ekin=alphas_Ekin_tot(init_rank:end_rank);
    alphas_vpll=alphas_vpll_tot(init_rank:end_rank);
    alphas_mm=alphas_mm_tot(init_rank:end_rank);
    alphas_momentum=alphas_momentum_tot(init_rank:end_rank);
    SAVENAME=strcat('initial_alphas_MB_D_distribution_m',num2str(n),'.mat')
    save(SAVENAME, 'mHe', 'ZHe', 'Nalphas_simulated', 'alphas_pos_x','alphas_pos_z', 'alphas_pos_phi', 'alphas_Ekin', 'alphas_vpll', 'alphas_mm', 'density_part_ratio', 'frac_fusion', 'alphas_momentum', 'PROCESS_NUMBER')
    init_rank=init_rank+Nalphas_simulated;
end
