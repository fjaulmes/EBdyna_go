    
    % -------------------------------------------------------------------
    % ----------------------- half time step ----------------------------
    % -------------------------------------------------------------------


    time=time+0.5*DELTA_TIME;

    
    update_GT_3D_collapse;

%    alphas_pphi0_prev=alphas_pphi0;
%       alphas_pphi0(INNER_PART)=alphas_pphi0(INNER_PART)-...
%         (h*ZHe).*(alphas_Rpos(INNER_PART).*(alphas_grad_Phi(INNER_PART))+2*(v_phi_step(INNER_PART)).*alphas_grad_psi_star(INNER_PART));
%     v_X_step=0.5*(v_X+v_X_prev);
%     v_Z_step=0.5*(v_Z+v_Z_prev);
%     alphas_psi_star_half_prev=0.5*(alphas_psi_star+alphas_psi_star_prev);
%     alphas_vphi_grad_psi_star=(1/h)*((alphas_psi_star-alphas_psi_star_part_prev))-alphas_Rpos.*alphas_Bfield.*(-v_X_step.*bZ+v_Z_step.*bX);
%     alphas_vphi_grad_psi_star=(1/h)*((alphas_psi_star-alphas_psi_star_part_prev)-delta_ps)+alphas_Rpos.*alphas_Bfield.*(-v_X_step.*bZ+v_Z_step.*bX);
    
%     psi_star_omega_map_half_step=psi_star_omega_map+(0.5)*D_psi_star_omega_map;

   % Canonical angular momentum evolution
    if CALCULATE_TRUE_PPHI==1
        v_phi_step=0.5*(v_phi+v_phi_prev);
        alphas_vphi_grad_psi_star=(v_phi_step.*alphas_grad_psi_star);
		delta_ps=1.5*delta_ps-0.5*delta_ps_prev;
		alphas_grad_Phi(INNER_PART)=-Ephi(INNER_PART)+(1/h)*delta_ps(INNER_PART)./alphas_Rpos(INNER_PART);
		alphas_pphi0(INNER_PART)=alphas_pphi0(INNER_PART)-...
			(h*ZHe).*(alphas_Rpos(INNER_PART).*(alphas_grad_Phi(INNER_PART))+alphas_vphi_grad_psi_star(INNER_PART));
	end
    %alphas_psi_star(INNER_PART)=interp2(1:size_r,scale_theta,psi_star_omega_map_half_step',alphas_psi(INNER_PART), alphas_omega(INNER_PART),'*linear');
    %delta_ps=sign(delta_ps).*min(abs(delta_ps),abs(alphas_psi_star-alphas_psi_star_half_prev));
    %alphas_grad_Phi(INNER_PART)=-Ephi(INNER_PART)+(1/h)*delta_ps(INNER_PART);
    %alphas_pphi0(INNER_PART)=alphas_pphi0(INNER_PART)-...
    %    (h*ZHe).*(alphas_Rpos(INNER_PART).*alphas_grad_Phi(INNER_PART)+alphas_vphi_grad_psi_star(INNER_PART));
%     estimate_half_time_step_positions_G;
%     INNER_PART=find(alphas_psi<=(size_r-1));
%     OUTER_PART=find(alphas_psi>(size_r-1));
%     alphas_omega=wrap2pi(alphas_theta-alphas_pos_phi_wrapped);
%     alphas_psi_star(INNER_PART)=interp2(1:size_r,scale_theta,psi_star_omega_map',alphas_psi(INNER_PART), alphas_omega(INNER_PART),'*linear');
%     if ~isempty(OUTER_PART)
%         alphas_psi_value_corr(OUTER_PART)=alphas_psi_value(OUTER_PART);
%     end
%     alphas_psi_value_corr(INNER_PART)=alphas_psiH_value(INNER_PART)+alphas_psi_star(INNER_PART);
%    
%     alphas_pphi0(INNER_PART)=0.5*(alphas_pphi0(INNER_PART)+((mHe/eV)*alphas_Rpos_int(INNER_PART)).*v_phi_step(INNER_PART)-(ZHe)*alphas_psi_value_corr(INNER_PART));
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
    alphas_pos_phi_wrapped=wrap2pi(alphas_pos_phi);
    
    % limited map range
    alphas_pos_x=min(alphas_pos_x,scale_X(end-3));
    alphas_pos_z=min(alphas_pos_z,scale_Z(end-3));
    alphas_pos_x=max(alphas_pos_x,scale_X(4));
    alphas_pos_z=max(alphas_pos_z,scale_Z(4));


    %     toc
    % estimating psi and theta values
    interpolate_theta_psi_fromXZ;
    alphas_omega=wrap2pi(alphas_theta-alphas_pos_phi_wrapped);
    % sorting out which particles are inside mixing region and which are out
	% clear INNER_PART OUTER_PART;
    INNER_PART=find(alphas_psi<=(size_r-1));
    OUTER_PART=find(alphas_psi>(size_r-1));       
    [IL3D_1 IL3D_2 IL3D_3 IL3D_4 IL3D_5 IL3D_6 IL3D_7 IL3D_8 slopex slopey slopez] = ...
        build_3Dinterp_indexarrays(scale_phi, scale_theta, scale_psi, DPHI,DTHETA,1,alphas_pos_phi_wrapped(INNER_PART), alphas_theta(INNER_PART), alphas_psi(INNER_PART));
    
    alphas_Epot_prev(INNER_PART)=interp2_omega_map(scale_psi(1:size_r),scale_theta,E_potential_omega_map,1,DTHETA,alphas_psi, alphas_omega,INNER_PART);
    if PSI_STAR_SIGN>0
        alphas_psi_star_prev(INNER_PART)=max(interp2_omega_map(scale_psi(1:size_r),scale_theta,psi_star_omega_map,1,DTHETA,alphas_psi, alphas_omega,INNER_PART),0);
    else
        alphas_psi_star_prev(INNER_PART)=min(interp2_omega_map(scale_psi(1:size_r),scale_theta,psi_star_omega_map,1,DTHETA,alphas_psi, alphas_omega,INNER_PART),0);
    end
    % alphas_Epot_prev(INNER_PART)=lininterp3( E_potential_PR_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
    % alphas_psi_star_prev(INNER_PART)=max(interp2(1:size_r,scale_theta,psi_star_omega_map',alphas_psi(INNER_PART), alphas_omega(INNER_PART),'*linear'),0);    
    if ~isempty(OUTER_PART)
        alphas_psi_star_prev(OUTER_PART)=zeros(length(OUTER_PART),1);
        alphas_Epot_prev(OUTER_PART)=zeros(length(OUTER_PART),1);
    end

    update_E_b_fields;
    % interpolate_theta_psi_fromXZ;

    %alphas_psi=interp2_XZ(interp_x,interp_z,psi_norm_XZsmall_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
    %INNER_PART=find(alphas_psi<=(size_r-1));
    %OUTER_PART=find(alphas_psi>(size_r-1));
%     [IL3D_1 IL3D_2 IL3D_3 IL3D_4 IL3D_5 IL3D_6 IL3D_7 IL3D_8 slopex slopey slopez] = ...
%             build_3Dinterp_indexarrays(scale_phi, scale_theta, scale_psi, DPHI,DTHETA,1,alphas_pos_phi_wrapped(INNER_PART), alphas_theta(INNER_PART), alphas_psi(INNER_PART));

%     alphas_Epot_step_prev=alphas_Epot;
%     alphas_psi_star_part_prev=alphas_psi_star;
    
    alphas_Epot(INNER_PART)=interp2_omega_map(scale_psi(1:size_r),scale_theta,E_potential_omega_map,1,DTHETA,alphas_psi, alphas_omega,INNER_PART);
    if PSI_STAR_SIGN>0
        alphas_psi_star(INNER_PART)=max(interp2_omega_map(scale_psi(1:size_r),scale_theta,psi_star_omega_map,1,DTHETA,alphas_psi, alphas_omega,INNER_PART),0);
    else
        alphas_psi_star(INNER_PART)=min(interp2_omega_map(scale_psi(1:size_r),scale_theta,psi_star_omega_map,1,DTHETA,alphas_psi, alphas_omega,INNER_PART),0);
    end
    % alphas_Epot(INNER_PART)=lininterp3( E_potential_PR_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
    % alphas_psi_star(INNER_PART)=max(interp2(1:size_r,scale_theta,psi_star_omega_map',alphas_psi(INNER_PART), alphas_omega(INNER_PART),'*linear'),0);
%     delta_E=v_X.*EX+v_Z.*EZ+v_phi.*Ephi;
    if ~isempty(OUTER_PART)
        alphas_Epot(OUTER_PART)=zeros(length(OUTER_PART),1);
        alphas_psi_star(OUTER_PART)=zeros(length(OUTER_PART),1);
    end

%     alphas_Epot_step_prev=0.5*(alphas_Epot+alphas_Epot_step_prev);

    % amplitude of the B field and toroidal gradients at local positions
    if ~isempty(OUTER_PART)
        if CALCULATE_TRUE_PPHI==1
            alphas_grad_Phi(OUTER_PART)=zeros(length(OUTER_PART),1);
            alphas_grad_psi_star(OUTER_PART)=zeros(length(OUTER_PART),1);
	    end
        delta_ps(OUTER_PART)=zeros(length(OUTER_PART),1);
        delta_E(OUTER_PART)=zeros(length(OUTER_PART),1);
    end
%     alphas_grad_Phi(INNER_PART)=lininterp3( grad_Phi_tor_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
    if CALCULATE_TRUE_PPHI==1
        alphas_grad_psi_star(INNER_PART)=lininterp3( grad_psi_star_map_phi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
    end
%     Ephi_prev=Ephi;
    update_Gfields_collapse;
%     Ephi_half=0.5*(Ephi+Ephi_prev);
    
    % psi star evolution
    delta_ps_prev=delta_ps;
    delta_ps(INNER_PART)=(alphas_psi_star(INNER_PART)-alphas_psi_star_prev(INNER_PART));
    %delta_ps_prev=1.5*delta_ps-0.5*delta_ps_prev;
%     alphas_grad_Phi(INNER_PART)=-Ephi(INNER_PART)+(1/h)*delta_ps(INNER_PART);
%     alphas_grad_Phi(INNER_PART)=-Ephi(INNER_PART);

    % Total energy evolution
    delta_E(INNER_PART)=ZHe*(alphas_Epot(INNER_PART)-alphas_Epot_prev(INNER_PART)+(v_phi(INNER_PART)./alphas_Rpos_int(INNER_PART)).*delta_ps(INNER_PART));
%     alphas_Etot_prev=alphas_Etot;
    alphas_Etot(INNER_PART)=alphas_Etot(INNER_PART)+delta_E(INNER_PART);
    

    % Kinetic energy evolution
%     alphas_Etot_prev=0.5*(alphas_Etot_prev+alphas_Etot);
%     alphas_Ekin_prev=alphas_Ekin;
%     alphas_Ekin=max(alphas_Etot_prev-ZHe*(alphas_Epot_step_prev),0);
    alphas_Ekin_prev=alphas_Ekin;
    alphas_Ekin(INNER_PART)=max(alphas_Etot(INNER_PART)-ZHe*(alphas_Epot(INNER_PART)),MIN_ALPHAS_EKIN);
	
	% work of the centrifugal force
    alphas_Ekin=max(alphas_Ekin+((mHe/eV)*h)*v_X.*alphas_Rpos_int.*alphas_momentum.^2,MIN_ALPHAS_EKIN);


    %evaluating the time step value
    v_X_step=(1.5*v_X-0.5*v_X_prev);
    v_Z_step=(1.5*v_Z-0.5*v_Z_prev);
    v_phi_step=(1.5*v_phi-0.5*v_phi_prev);

    % an estimate of the half time step Ekin
%     alphas_Ekin=alphas_Ekin+ZHe*h*(v_X.*EX+v_Z.*EZ+v_phi.*Ephi);

    
    v_X_prev=v_X;
    v_Z_prev=v_Z;
    v_phi_prev_prev=v_phi_prev;
    v_phi_prev=v_phi;

