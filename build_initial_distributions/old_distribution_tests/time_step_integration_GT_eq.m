    
    % -------------------------------------------------------------------
    % ----------------------- half time step ----------------------------
    % -------------------------------------------------------------------

    
    time=time+0.5*DELTA_TIME;

    
    update_GT_3D_eq;

%     estimate_half_time_step_positions_eq_G;

%     adapt_speed_pphi_G;
%     v_phi=v_phi_tilde;
%     adapt_speed_Ekin_G;
%     v_phi=v_phi_tilde;
    
    alphas_pphi0_prev=alphas_pphi0;

% 
    % -------------------------------------------------------------------
    % ----------------------- next time step ----------------------------
    % -------------------------------------------------------------------
    time=time+0.5*DELTA_TIME;
    
%     tic
%     update_Gposition_tokamak_coordinates;
    alphas_Rpos=alphas_pos_x+R0;
    alphas_pos_x=alphas_pos_x+h*v_X;
    alphas_Rpos=0.5*(alphas_Rpos+alphas_pos_x+R0);
    alphas_pos_phi=alphas_pos_phi+h*(v_phi./alphas_Rpos);
    alphas_pos_z=alphas_pos_z+h*v_Z;
    alphas_Rpos=alphas_pos_x+R0;
    alphas_pos_phi_wrapped=wrap2pi(alphas_pos_phi);
    
    % limited map range
    alphas_pos_x=min(alphas_pos_x,scale_X(end-3));
    alphas_pos_z=min(alphas_pos_z,scale_Z(end-3));
    alphas_pos_x=max(alphas_pos_x,scale_X(4));
    alphas_pos_z=max(alphas_pos_z,scale_Z(4));


    %     toc
    % estimating psi and theta values
    interpolate_theta_psi_fromXZ;


    % amplitude of the B field and potential at local positions
    alphas_Bfield=interp2_XZ(interp_x,interp_z,Btot_XZ_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
%     
% 
    %evaluating the time step value
    v_X_step=(1.5*v_X-0.5*v_X_prev);
    v_Z_step=(1.5*v_Z-0.5*v_Z_prev);
    v_phi_step=(1.5*v_phi-0.5*v_phi_prev);
    
    v_X_prev=v_X;
    v_Z_prev=v_Z;
    v_phi_prev_prev=v_phi_prev;
    v_phi_prev=v_phi;

    % Kinetic energy evolution
    alphas_Ekin_prev=alphas_Ekin;
%     alphas_Ekin=max(alphas_Etot-ZHe*alphas_Epot,0);