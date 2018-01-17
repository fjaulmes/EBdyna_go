function [ output_args ] = correct_values_struct( input_args )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
error('Not finished writing')
%% Correct values
    %calculating integer time step velocity values
    if CALCULATE_TRUE_PPHI==0
        update_GT_3D_eq;
        v_X_step=(0.5*v_X+0.5*v_X_prev);
        v_Z_step=(0.5*v_Z+0.5*v_Z_prev);
        v_phi_step=(0.5*v_phi+0.5*v_phi_prev);
        v_X=v_X_prev;
        v_Z=v_Z_prev;
        v_phi=v_phi_prev;
    else
        alphas_pphi0=alphas_pphi0_prev+alphas_Delta_pphi;
        alphas_pphi0_prev=alphas_pphi0;
        alphas_Delta_pphi=zeros(Nalphas_simulated,1);
    end
    %applying vphi correction
    if PPHI_CORR_FACTOR~=0
        v_phi_step_recalc=(alphas_pphi0+(ZHe)*alphas_psi_value)./(alphas_Rpos);
        v_phi_step_recalc=(eV/mHe)*v_phi_step_recalc;
        adapt_speed_pphi_G;
    end
    %applying Ekin correction
    if ECOEF~=0
        adapt_speed_Ekin_G;
    end
    %applying vphi correction
    if PPHI_CORR_FACTOR~=0
        update_GT_3D_eq;
        v_X_step=(0.5*v_X+0.5*v_X_prev);
        v_Z_step=(0.5*v_Z+0.5*v_Z_prev);
        v_phi_step=(0.5*v_phi+0.5*v_phi_prev);
        v_X=v_X_prev;
        v_Z=v_Z_prev;
        v_phi=v_phi_prev;
        adapt_speed_pphi_G;
    end

end

