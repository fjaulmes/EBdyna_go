
    
    disp('frame_rank = ');
    disp(frame_rank);
    
    reference_frame=min((frame_rank-1)*10+1,1001)
    
    f=reference_frame;
%     f0=f;
%     
%     rx_value=rx_evol_interp(f)
%     xpos_sep=interp1(radial_r_value_flux,1:Nradial,rx_value);
%     ksi0_value=ksi0_evol_interp(f);
%     
%     psi_limit13=max(interp1(1:Nradial,psi_star_initial_profile,rx_precise),0)
%     
%     x_pos_sep_final=interp1(psi_star_final_profile,1:size_r,psi_limit13);
    
    
    
    psi_star_omega_map_half(:,:)=psi_star_2D_evol_lin(frame_rank,:,:);
    % Using symmetry to reconstruct a poloidal turn
    psi_star_omega_map_rank=zeros(size_r,NB_THETA);
    psi_star_omega_map_rank(:,1:round(0.5*NB_THETA))=psi_star_omega_map_half(:,:);
    psi_star_omega_map_rank(:,round(0.5*NB_THETA):NB_THETA)=psi_star_omega_map_half(:,end:-1:1);
    
    psi_star_omega_map=zeros(size_r,NB_THETA);
    psi_star_omega_map(:,:)=psi_star_omega_map_rank(:,1:NB_THETA);
    psi_star_PR_map=psi_star_omega_map';
    
%     if (f<10)
%         frame_name='00';
%     elseif (f<100)
%         frame_name='0';
%     else
%         frame_name='';
%     end
%     frame_name=strcat(frame_name,num2str(f));
%     filename='../calculate_Epot/reconnection_maps/t0';
%     filename=strcat(filename,frame_name,'.mat');
%     
%     load(filename);
    
    psi_star_dot_PR_map=zeros(NP,size_r);
    E_potential_PR_map=zeros(NP,size_r);
    BstarX_PR_map=zeros(NP,size_r);
    BstarZ_PR_map=zeros(NP,size_r);
    
%     BstarX_PR_map_phi(NB_PHI,:,:)=BstarX_PR_map_phi(1,:,:);
%     BstarZ_PR_map_phi(NB_PHI,:,:)=BstarZ_PR_map_phi(1,:,:);
    
    
    
    % positions phi_rank=1 and phi_rank=NB_PHI have to be treated differently
    phi_index=1;
    clear psi_star_PR_map_rank_next psi_star_PR_map_rank_prev;
    minus_phi_rank=2;
    psi_star_PR_map_rank_next(:,:)=rotate_map_phi(psi_star_omega_map_rank,minus_phi_rank);
    psi_star_PR_map_rank_next=psi_star_PR_map_rank_next';
    
    minus_phi_rank=NB_THETA-1;
    psi_star_PR_map_rank_prev(:,:)=rotate_map_phi(psi_star_omega_map_rank,minus_phi_rank);
    psi_star_PR_map_rank_prev=psi_star_PR_map_rank_prev';
    
%     BstarX_PR_map=zeros(NP,size_r);
%     BstarZ_PR_map=zeros(NP,size_r);
%     
%     BstarX_PR_map(:,:)=BstarX_PR_map_phi(phi_rank,:,:);
%     BstarZ_PR_map(:,:)=BstarZ_PR_map_phi(phi_rank,:,:);
    phi_index=1;
    if PHI_OMEGA_RATIO==1
        phi_rank=phi_index
    else
        phi_rank=PHI_OMEGA_RATIO*phi_index-1
    end   
    calculate_Btot_phi_rank;
    
    if CALCULATE_VD_DATA_FILE==1
        calculate_vD_drift_speeds_phi_rank;
        vD_X_map_phi(phi_index,:,:)=vD_X_PR_map;
        vD_Z_map_phi(phi_index,:,:)=vD_Z_PR_map;
        vD_phi_map_phi(phi_index,:,:)=vD_phi_PR_map;
        vD_X_map_phi(NB_PHI,:,:)=vD_X_PR_map;
        vD_Z_map_phi(NB_PHI,:,:)=vD_Z_PR_map;
        vD_phi_map_phi(NB_PHI,:,:)=vD_phi_PR_map;
    end
    
    calculate_grad_psi_phi_rank;
    
    Btot_map_phi(phi_index,:,:)=Btot_PR_map;
    bX_map_phi(phi_index,:,:)=BpolX_PR_map./Btot_PR_map;
    bZ_map_phi(phi_index,:,:)=BpolZ_PR_map./Btot_PR_map;
    bphi_map_phi(phi_index,:,:)=Btor_PR_map(:,1:size_r)./Btot_PR_map;
    grad_psi_star_map_phi(phi_index,:,:)=grad_psi_star_tor_PR_map;
    
    Btot_map_phi(NB_PHI,:,:)=Btot_PR_map;
    bX_map_phi(NB_PHI,:,:)=BpolX_PR_map./Btot_PR_map;
    bZ_map_phi(NB_PHI,:,:)=BpolZ_PR_map./Btot_PR_map;
    bphi_map_phi(NB_PHI,:,:)=Btor_PR_map(:,1:size_r)./Btot_PR_map;
    grad_psi_star_map_phi(NB_PHI,:,:)=grad_psi_star_tor_PR_map;
    
    
    
    % Looping around in the toroidal direction
    for (phi_index=2:NB_PHI-1)
        %    for (phi_rank=12:12)
        
        disp('**************************************************************');
        if PHI_OMEGA_RATIO==1
            phi_rank=phi_index
        else
            phi_rank=PHI_OMEGA_RATIO*phi_index-1
        end
        
        psi_star_omega_map(:,:)=rotate_map_phi(psi_star_omega_map_rank,phi_rank);
        psi_star_PR_map=psi_star_omega_map';
        
        clear psi_star_PR_map_rank_next psi_star_PR_map_rank_prev;
        
%         minus_phi_rank=max(NB_THETA-(phi_rank-1)+1-1,1);
        psi_star_PR_map_rank_next(:,:)=rotate_map_phi(psi_star_omega_map_rank,phi_rank+1);
        psi_star_PR_map_rank_next=psi_star_PR_map_rank_next';
        
%         minus_phi_rank=max(NB_THETA-(phi_rank-1)+1+1,1);
        psi_star_PR_map_rank_prev(:,:)=rotate_map_phi(psi_star_omega_map_rank,phi_rank-1);
        psi_star_PR_map_rank_prev=psi_star_PR_map_rank_prev';
        
        
%         BstarX_PR_map=zeros(NP,size_r);
%         BstarZ_PR_map=zeros(NP,size_r);
%         
%         BstarX_PR_map(:,:)=BstarX_PR_map_phi(phi_rank,:,:);
%         BstarZ_PR_map(:,:)=BstarZ_PR_map_phi(phi_rank,:,:);
        
        calculate_Btot_phi_rank;
        if CALCULATE_VD_DATA_FILE==1
            calculate_vD_drift_speeds_phi_rank;
            vD_X_map_phi(phi_index,:,:)=vD_X_PR_map;
            vD_Z_map_phi(phi_index,:,:)=vD_Z_PR_map;
            vD_phi_map_phi(phi_index,:,:)=vD_phi_PR_map;
        end
        calculate_grad_psi_phi_rank;
        
        disp('****************calculating B and psi maps DONE!**********************');
        
        Btot_map_phi(phi_index,:,:)=Btot_PR_map;
        bX_map_phi(phi_index,:,:)=BpolX_PR_map./Btot_PR_map;
        bZ_map_phi(phi_index,:,:)=BpolZ_PR_map./Btot_PR_map;
        bphi_map_phi(phi_index,:,:)=Btor_PR_map(:,1:size_r)./Btot_PR_map;
        grad_psi_star_map_phi(phi_index,:,:)=grad_psi_star_tor_PR_map;
         
    end
    

    
    delta_energy_frame_rank=sum(energy_star_phi_evol(round(frame_rank),1:end-1),2)
    energy_frame_rank=sum(energy_phi_evol(round(frame_rank),1:end-1),2)
    B2_tot_evol(round(frame_rank))=energy_frame_rank;
    deltaB2_tot_evol(round(frame_rank))=delta_energy_frame_rank;

