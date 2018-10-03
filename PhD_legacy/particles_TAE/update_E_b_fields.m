% if (time<0)
%     Efield_X_map_phi_next=zeros(NB_PHI,NB_THETA,size_r);
%     Efield_Z_map_phi_next=zeros(NB_PHI,NB_THETA,size_r);
%     E_potential_PR_map_phi=zeros(NB_PHI,NB_THETA,size_r);
%     iEpot_map_phi=zeros(NB_PHI,NB_THETA,size_r);
%     Epot_map_phi_next=E_potential_PR_map_phi;
% else


time_oscill=mod(time,TAE_PERIOD);
if time<SIMULATION_TIME
    frame_rank_precise=(NB_FRAME)*(time_oscill/TAE_PERIOD);
else
    frame_rank_precise=0
end
frame_rank_prev=floor(frame_rank_precise);
frame_rank_next=ceil(frame_rank_precise);

if (frame_rank_prev==0) && (f_counter>1)
    frame_rank_prev=NB_FRAME;
    frame_rank_next=1;
end

% frame_rank_prev=mod(frame_rank_prev,NB_FRAME+1);
% frame_rank_next=mod(frame_rank_next,NB_FRAME+1);
if (frame_rank_next==frame_rank_prev)
    frame_rank_next=frame_rank_next+1;
end

if frame_rank_next~=frame_rank_next_prev
    disp('**********************************************************')
    f_counter=f_counter+1

    if f_counter>1
        %frame_rank=frame_rank_prev;
        Efield_X_map_phi_prev=Efield_X_map_phi_next;
        Efield_Z_map_phi_prev=Efield_Z_map_phi_next;
        %     vExB_phi_map_phi_prev=vExB_phi_map_phi_next;
        BsX_map_phi_prev=BsX_map_phi_next;
        BsZ_map_phi_prev=BsZ_map_phi_next;
        %     bphi_map_phi_prev=bphi_map_phi_next;
        B_PR_map_phi_prev=B_PR_map_phi_next;
        
%         psi_star_map_prev=psi_star_map_next;
%         iEpot_map_phi_prev=iEpot_map_phi_next;
%         Epot_map_phi_prev=Epot_map_phi_next;
        

    else
        
        % set initial perturbed fields to 0
        BsX_map_phi_prev=0*BsX_map_phi;
        BsZ_map_phi_prev=0*BsZ_map_phi;
        B_PR_map_phi_prev=Btot_map_phi_ini;

        Efield_X_map_phi_prev=0*Efield_X_map_phi;
        Efield_Z_map_phi_prev=0*Efield_Z_map_phi;
%         psi_star_map_prev=0*Epot_map_phi;
%         Epot_map_phi_prev=0*Epot_map_phi;
%         iEpot_map_phi_prev=0*iEpot_map_phi;
        

    end
%     if f_counter==2
%         disp('initialize global Ekin calculation')
%         Epart_tot_vD_rel_ini=0;
%         Epart_tot_vD_rel_ini_prev=0;
%         Ekin_vD=alphas_Ekin;
%         Ekin_tot_part=NB_PART_RESCALE*eV*sum(alphas_Ekin)
%     end
    
    frame_rank=frame_rank_next;
    %         if frame_rank<=(NB_FRAME*NB_OSCILLATIONS+1)
    disp('loading new frame rank#');
    disp(frame_rank);
    load_Bmaps_frame_rank;
    load_Emaps_frame_rank;
    
    Efield_X_map_phi_next=Efield_X_map_phi;
    Efield_Z_map_phi_next=Efield_Z_map_phi;
%     Epot_map_phi_next=Epot_map_phi;
%     iEpot_map_phi_next=iEpot_map_phi;

    BsX_map_phi_next=BsX_map_phi;
    BsZ_map_phi_next=BsZ_map_phi;
    Btot_map_phi=sqrt((BsX_map_phi_next+BpolX_map_phi_ini).^2+(BsZ_map_phi_next+BpolZ_map_phi_ini).^2+(Bsphi_map_phi+Btor_map_phi_ini).^2);
%     Btot_map_phi=sqrt((BsX_map_phi_next+BpolX_map_phi_ini).^2+(BsZ_map_phi_next+BpolZ_map_phi_ini).^2+(Btor_map_phi_ini).^2);

    B_PR_map_phi_next=Btot_map_phi;
%     for p=1:NB_PHI
%         psi_star_map_next(p,:,:)=-(kTAE/omega_TAE)*squeeze(Epot_map_phi(p,:,:)).*Rpos_PR_map(:,pTAE_inf:pTAE_sup);
%     end
    
    %     disp('frame_rank_precise = ');
    %     disp(frame_rank_precise);
    disp('time = ');
    disp(time);
    frame_rank_next_prev=frame_rank_next;
    
    Efield_X_map_phi=Efield_X_map_phi_prev;
    Efield_Z_map_phi=Efield_Z_map_phi_prev;
    BsX_map_phi=BsX_map_phi_prev;
    BsZ_map_phi=BsZ_map_phi_prev;

%     E_potential_map_phi=Epot_map_phi_prev;
%     psi_star_phi_map=psi_star_map_prev;
    % grad_Phi_tor_map_phi=grad_Phi_map_phi_prev;
%     iEpot_map_phi=iEpot_map_phi_prev;
    Btot_map_phi=B_PR_map_phi_prev;
    
    D_Efield_X_map_phi=delta_t_coef*(Efield_X_map_phi_next-Efield_X_map_phi_prev);
    D_Efield_Z_map_phi=delta_t_coef*(Efield_Z_map_phi_next-Efield_Z_map_phi_prev);
    D_bX_map_phi=delta_t_coef*(BsX_map_phi_next-BsX_map_phi_prev);
    D_bZ_map_phi=delta_t_coef*(BsZ_map_phi_next-BsZ_map_phi_prev);
    D_Btot_map_phi=delta_t_coef*(B_PR_map_phi_next-B_PR_map_phi_prev);
%     D_E_potential_map_phi=delta_t_coef*(Epot_map_phi_next-Epot_map_phi_prev);
    % D_grad_Phi_tor_map_phi=delta_t_coef*(grad_Phi_map_phi_next-grad_Phi_map_phi_prev);
%     D_iEpot_map_phi=delta_t_coef*(iEpot_map_phi_next-iEpot_map_phi_prev);
    
    % we also use this delta to calculate half time steps
%     D_psi_star_phi_map=delta_t_coef*(psi_star_map_next-psi_star_map_prev);

    
else
    
    Efield_X_map_phi=Efield_X_map_phi+D_Efield_X_map_phi;
    Efield_Z_map_phi=Efield_Z_map_phi+D_Efield_Z_map_phi;
    BsX_map_phi=BsX_map_phi+D_bX_map_phi;
    BsZ_map_phi=BsZ_map_phi+D_bZ_map_phi;
    Btot_map_phi=Btot_map_phi+D_Btot_map_phi;
%     E_potential_map_phi=E_potential_map_phi+D_E_potential_map_phi;
    % grad_Phi_tor_map_phi=grad_Phi_tor_map_phi+D_grad_Phi_tor_map_phi;
%     iEpot_map_phi=iEpot_map_phi+D_iEpot_map_phi;
    
%     psi_star_phi_map=psi_star_phi_map+D_psi_star_phi_map;
end



