%%
reset_data_analysis_environment;
NB_PROCESS=128
% initial distributions
SAVENAME='./postRMP_NBI60keV_all.mat'

n=1;
INPUTNAME=strcat('post_RMP_NBI60keV_',num2str(n),'.mat')
load(INPUTNAME)



for n=2:NB_PROCESS
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
    alphas_ejected_tot=[alphas_ejected ];
    alphas_momentum_tot=[alphas_momentum ];

    INPUTNAME=strcat('post_RMP_NBI60keV_',num2str(n),'.mat')
    load(INPUTNAME)
    
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
    alphas_momentum_p2=alphas_momentum;
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
    v_X=[v_X_tot ; v_X_p2 ];
    v_Z=[v_Z_tot ;v_Z_p2 ];
    v_phi=[v_phi_tot ; v_phi_p2 ];
    alphas_ejected=[alphas_ejected_tot ; alphas_ejected_p2 ];
    alphas_momentum=[alphas_momentum_tot ; alphas_momentum_p2];
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


save(SAVENAME,'mHe', 'ZHe','Nalphas_simulated','alphas_ejected','alphas_momentum','alphas_Ekin', 'alphas_mm', 'alphas_pphi0', 'pos_X_gc','pos_Z_gc','alphas_pos_x', 'alphas_pos_z', 'alphas_pos_phi', 'alphas_vpll', 'alphas_psi', 'v_X', 'v_Z', 'v_phi');


