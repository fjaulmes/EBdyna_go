    
    % -------------------------------------------------------------------
    % ----------------------- half time step ----------------------------
    % -------------------------------------------------------------------
    v_X_prev=v_X;
    v_Z_prev=v_Z;
    v_phi_prev_prev=v_phi_prev;
    v_phi_prev=v_phi;

    
    time=time+0.5*DELTA_TIME;

    
    update_GT_3D_eq_centrifug;

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
    alphas_Rpos_int=0.5*(alphas_Rpos+alphas_pos_x+R0);
    alphas_pos_phi=alphas_pos_phi+h*(v_phi./alphas_Rpos_int);
    alphas_pos_z=alphas_pos_z+h*v_Z;
    alphas_Rpos=alphas_pos_x+R0;
	alphas_pos_phi=double(alphas_pos_phi);
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
    if APPLY_RMP_FIELD==1
        [IL3D_1 IL3D_2 IL3D_3 IL3D_4 IL3D_5 IL3D_6 IL3D_7 IL3D_8 slopex slopey slopez] = ...
            build_3Dinterp_indexarrays(scale_phi, scale_theta, scale_psi, DPHI,DTHETA,DPSIVAL,alphas_pos_phi_wrapped, alphas_theta, alphas_psi);
        BR_tilde=lininterp3( BR_RMP,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
        BZ_tilde=lininterp3( BZ_RMP,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
        Bphi_tilde=lininterp3( Bphi_RMP,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
        if CALCULATE_TRUE_PPHI==1
            alphas_dAphi_dphi=lininterp3( dAphi_dphi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
        end
    end


    %evaluating the time step value

    if CALCULATE_TRUE_PPHI==1
        v_X_prev=v_X;
        v_Z_prev=v_Z;
        v_phi_prev_prev=v_phi_prev;
        v_phi_prev=v_phi;
        update_GT_3D_eq;
        v_X_step=(0.5*v_X+0.5*v_X_prev);
        v_Z_step=(0.5*v_Z+0.5*v_Z_prev);
        v_phi_step=(0.5*v_phi+0.5*v_phi_prev);
        v_X=v_X_prev;
        v_Z=v_Z_prev;
        v_phi=v_phi_prev;
    else
        v_X_step=(1.5*v_X-0.5*v_X_prev);
        v_Z_step=(1.5*v_Z-0.5*v_Z_prev);
        v_phi_step=(1.5*v_phi-0.5*v_phi_prev);
    end
    
    v_X_prev=v_X;
    v_Z_prev=v_Z;
    v_phi_prev_prev=v_phi_prev;
    v_phi_prev=v_phi;

    % pphi and energy evolution
    alphas_Ekin_prev=alphas_Ekin;
    if CALCULATE_TRUE_PPHI==1
        % this calculates the Delta pphi between half time step values
        alphas_Delta_pphi=alphas_Delta_pphi+(ZHe*h)*(v_phi_step.*alphas_dAphi_dphi);
    end
    

    % Kinetic energy evolution
    alphas_Ekin_prev=alphas_Ekin;
 	% work of the centrifugal force
    alphas_Ekin=max(alphas_Ekin+((mHe/eV)*h)*v_X.*alphas_Rpos_int.*alphas_momentum.^2,MIN_ALPHAS_EKIN);

%     alphas_Ekin=max(alphas_Etot-ZHe*alphas_Epot,0);