
    if time<(tau_cr/FAST_SAWTOOTH)
        frame_rank_precise=1000*(time*FAST_SAWTOOTH/tau_cr)+1;
    else
        frame_rank_precise=1001
    end
    frame_rank_prev=10*floor((frame_rank_precise-1)/10)+1;
    frame_rank_next=10*ceil((frame_rank_precise-1)/10)+1;
    if (frame_rank_next==frame_rank_prev)
        frame_rank_next=frame_rank_next+10;
    end
    
    if frame_rank_next~=frame_rank_next_prev
        
        %frame_rank=frame_rank_prev;
        Efield_X_map_phi_prev=Efield_X_map_phi_next;
        Efield_Z_map_phi_prev=Efield_Z_map_phi_next;
        %     vExB_phi_map_phi_prev=vExB_phi_map_phi_next;
        bX_map_phi_prev=bX_map_phi_next;
        bZ_map_phi_prev=bZ_map_phi_next;
        %     bphi_map_phi_prev=bphi_map_phi_next;
        B_PR_map_phi_prev=B_PR_map_phi_next;
        
        psi_star_map_prev=psi_star_map_next;
        grad_Phi_map_phi_prev=grad_Phi_map_phi_next;
        Epot_map_phi_prev=Epot_map_phi_next;

        grad_psi_star_map_phi_prev=grad_psi_star_map_phi_next;
        
        frame_rank=frame_rank_next;
        if frame_rank<=1001
            disp('loading new frame rank#');
            disp(frame_rank);
            load_Bmaps_frame_rank_gpu;
            load_Emaps_frame_rank_gpu;
            update_psi_star_frame_rank;

			E_potential_omega_map(:,:)=Epot_evol(round((frame_rank-1)/10+1),:,:);
            
            Efield_X_map_phi_next=FAST_SAWTOOTH*Efield_X_map_phi;
            Efield_Z_map_phi_next=FAST_SAWTOOTH*Efield_Z_map_phi;
            %         vExB_phi_map_phi_next=FAST_SAWTOOTH*vExB_phi_map_phi;
            % grad_Phi_map_phi_next=FAST_SAWTOOTH*grad_Phi_tor_map_phi;
            Epot_map_phi_next=FAST_SAWTOOTH*E_potential_omega_map;
            grad_psi_star_map_phi_next=grad_psi_star_map_phi;
        else
            Efield_X_map_phi_next=zeros(NB_PHI,NB_THETA,size_r);
            Efield_Z_map_phi_next=zeros(NB_PHI,NB_THETA,size_r);
            %         vExB_phi_map_phi_next=zeros(NB_PHI,NB_THETA,size_r);
            % grad_Phi_map_phi_next=zeros(NB_PHI,NB_THETA,size_r);
            Epot_map_phi_next=zeros(size_r,NB_THETA);
            grad_psi_star_map_phi_next=zeros(NB_PHI,NB_THETA,size_r);
       end
        bX_map_phi_next=bX_map_phi;
        bZ_map_phi_next=bZ_map_phi;
        %     bphi_map_phi_next=bphi_map_phi;
        B_PR_map_phi_next=Btot_map_phi;
        psi_star_map_next=psi_star_omega_map_rank;
        
        disp('time = ');
        disp(time);
        frame_rank_next_prev=frame_rank_next;
        
        Efield_X_map_phi_gpu=gpuArray(Efield_X_map_phi_prev);
        Efield_Z_map_phi_gpu=gpuArray(Efield_Z_map_phi_prev);
        bX_map_phi_gpu=gpuArray(bX_map_phi_prev);
        bZ_map_phi_gpu=gpuArray(bZ_map_phi_prev);
        Btot_map_phi_gpu=gpuArray(B_PR_map_phi_prev);
        E_potential_omega_map_gpu=gpuArray(Epot_map_phi_prev);
        psi_star_omega_map_gpu=gpuArray(psi_star_map_prev);
        % grad_Phi_tor_map_phi_gpu=gpuArray(grad_Phi_map_phi_prev);
        grad_psi_star_map_phi_gpu=gpuArray(grad_psi_star_map_phi_prev);

        D_Efield_X_map_phi=gpuArray(delta_t_coef*(Efield_X_map_phi_next-Efield_X_map_phi_prev));
        D_Efield_Z_map_phi=gpuArray(delta_t_coef*(Efield_Z_map_phi_next-Efield_Z_map_phi_prev));
        D_bX_map_phi=gpuArray(delta_t_coef*(bX_map_phi_next-bX_map_phi_prev));
        D_bZ_map_phi=gpuArray(delta_t_coef*(bZ_map_phi_next-bZ_map_phi_prev));
        D_Btot_map_phi=gpuArray(delta_t_coef*(B_PR_map_phi_next-B_PR_map_phi_prev));
        D_E_potential_omega_map=gpuArray(delta_t_coef*(Epot_map_phi_next-Epot_map_phi_prev));
        % D_E_potential_PR_map_phi=gpuArray(delta_t_coef*(Epot_map_phi_next-Epot_map_phi_prev));
        % D_grad_Phi_tor_map_phi=gpuArray(delta_t_coef*(grad_Phi_map_phi_next-grad_Phi_map_phi_prev));
        D_grad_psi_star_tor_map_phi=gpuArray(delta_t_coef*(grad_psi_star_map_phi_next-grad_psi_star_map_phi_prev));

        % we also use this delta to calculate half time steps
        D_psi_star_omega_map=gpuArray(delta_t_coef*(psi_star_map_next-psi_star_map_prev));

    end
    
    
    Efield_X_map_phi_gpu=Efield_X_map_phi_gpu+D_Efield_X_map_phi;
    Efield_Z_map_phi_gpu=Efield_Z_map_phi_gpu+D_Efield_Z_map_phi;
    bX_map_phi_gpu=bX_map_phi_gpu+D_bX_map_phi;
    bZ_map_phi_gpu=bZ_map_phi_gpu+D_bZ_map_phi;
    Btot_map_phi_gpu=Btot_map_phi_gpu+D_Btot_map_phi;
    E_potential_omega_map_gpu=E_potential_omega_map_gpu+D_E_potential_omega_map;
    % grad_Phi_tor_map_phi_gpu=grad_Phi_tor_map_phi_gpu+D_grad_Phi_tor_map_phi;
    grad_psi_star_map_phi_gpu=grad_psi_star_map_phi_gpu+D_grad_psi_star_tor_map_phi;
    psi_star_omega_map_gpu=psi_star_omega_map_gpu+D_psi_star_omega_map;
    
    % psi_star_omega_map=gather(psi_star_omega_map_gpu);

    frame_rank_precise_prev=frame_rank_precise;

