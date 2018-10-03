NB_PROCESS=16

for frame_number=200:200:19800
    
    % post simulation results
    SAVENAME=strcat('./W40_fc1p6h1p6_',num2str(frame_number),'_all.mat')
    
    n=1;
    INPUTNAME=strcat('W40_fc1p6h1p6_',num2str(frame_number),'_',num2str(n),'.mat')
    load(INPUTNAME)
    
    
    for n=2:NB_PROCESS
%         if n~=10
        alphas_Ekin_tot=alphas_Ekin;
        alphas_mm_tot=alphas_mm;
        alphas_pphi0_tot=alphas_pphi0;
        alphas_pos_x_tot=alphas_pos_x;
        alphas_pos_z_tot=alphas_pos_z;
        alphas_pos_phi_tot=alphas_pos_phi;
        alphas_vpll_tot=alphas_vpll;
%         alphas_psi_tot=alphas_psi;
%         v_X_tot=v_X;
%         v_Z_tot=v_Z;
%         v_phi_tot=v_phi;
        alphas_ejected_tot=alphas_ejected;
        alphas_psi_value_corr_tot=alphas_psi_value_corr;
        alphas_omega_tot=alphas_omega;
        alphas_Etot_tot=alphas_Etot;
        pos_X_gc_tot=pos_X_gc;
        pos_Z_gc_tot=pos_Z_gc;
        
        INPUTNAME=strcat('W40_fc1p6h1p6_',num2str(frame_number),'_',num2str(n),'.mat')
        load(INPUTNAME)
        
        alphas_Ekin_p2=alphas_Ekin;
        alphas_mm_p2=alphas_mm;
        alphas_pphi0_p2=alphas_pphi0;
        alphas_pos_x_p2=alphas_pos_x;
        alphas_pos_z_p2=alphas_pos_z;
        alphas_pos_phi_p2=alphas_pos_phi;
        alphas_vpll_p2=alphas_vpll;
%         alphas_psi_p2=alphas_psi;
%         v_X_p2=v_X;
%         v_Z_p2=v_Z;
%         v_phi_p2=v_phi;
        alphas_ejected_p2=alphas_ejected;
        alphas_psi_value_corr_p2=alphas_psi_value_corr;
        alphas_omega_p2=alphas_omega;
        alphas_Etot_p2=alphas_Etot;
        pos_X_gc_p2=pos_X_gc;
        pos_Z_gc_p2=pos_Z_gc;
        
        
        
        
        alphas_Ekin=[alphas_Ekin_tot ; alphas_Ekin_p2 ];
        alphas_mm=[alphas_mm_tot ; alphas_mm_p2 ];
        alphas_pphi0=[alphas_pphi0_tot ; alphas_pphi0_p2 ];
        alphas_pos_x=[alphas_pos_x_tot ; alphas_pos_x_p2 ];
        alphas_pos_z=[alphas_pos_z_tot ; alphas_pos_z_p2 ];
        alphas_pos_phi=[alphas_pos_phi_tot ; alphas_pos_phi_p2 ];
        alphas_vpll=[alphas_vpll_tot ; alphas_vpll_p2 ];
%         alphas_psi=[alphas_psi_tot ; alphas_psi_p2 ];
%         v_X=[v_X_tot ; v_X_p2 ];
%         v_Z=[v_Z_tot ;v_Z_p2 ];
%         v_phi=[v_phi_tot ; v_phi_p2 ];
        alphas_ejected=[alphas_ejected_tot ; alphas_ejected_p2];
        alphas_psi_value_corr=[alphas_psi_value_corr_tot ; alphas_psi_value_corr_p2 ];
        alphas_omega=[alphas_omega_tot ; alphas_omega_p2 ];
        alphas_Etot=[alphas_Etot_tot ; alphas_Etot_p2 ];
        pos_X_gc=[pos_X_gc_tot ; pos_X_gc_p2 ];
        pos_Z_gc=[pos_Z_gc_tot ; pos_Z_gc_p2 ];
%         end
    end
    
    save(SAVENAME,'alphas_Ekin', 'alphas_mm', 'alphas_pphi0', 'alphas_pos_x', 'alphas_pos_z', 'alphas_pos_phi', 'alphas_vpll',  'alphas_ejected', 'alphas_psi_value_corr',  'time', 'frame_rank_precise', 'alphas_omega', 'alphas_Etot', 'pos_X_gc', 'pos_Z_gc');
    
    
end
disp('done....')

