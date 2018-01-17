
    
    disp('frame_rank = ');
    disp(frame_rank);
    
    reference_frame=min((frame_rank-1)*10+1,1001)
    
    f=reference_frame;
%     f0=f;
%     rx_value=rx_evol_interp(f)
%     ksi0_value=ksi0_evol_interp(f);
%     xpos_sep=interp1(radial_r_value_flux,1:Nradial,rx_value);
%     
%     psi_limit13=max(interp1(1:Nradial,Psih,rx_precise),0)
%     
%     x_pos_sep_final=interp1(Psih_final,1:Nradial,psi_limit13);
    
    
    Epot_omega_map(:,:)=Epot_evol(frame_rank,:,:);
    
    E_potential_PR_map_phi=zeros(NB_THETA,NB_THETA,size_r);
    for phi_rank=1:NB_THETA

        E_potential_PR_map_phi(phi_rank,:,:)=rotate_map_phi(Epot_omega_map,phi_rank)';
    end
    
    psi_star_dot_omega_map(:,:)=psi_star_dot_evol(frame_rank,:,:);
    
    psi_star_dot_PR_map_phi=zeros(NB_THETA,NB_THETA,size_r);
    for phi_rank=1:NB_THETA

        psi_star_dot_PR_map_phi(phi_rank,:,:)=rotate_map_phi(psi_star_dot_omega_map,phi_rank)';
    end
    
    
    psi_star_dot_PR_map=zeros(NP,size_r);
    E_potential_PR_map=zeros(NP,size_r);
%     BstarX_PR_map=zeros(NP,size_r);
%     BstarZ_PR_map=zeros(NP,size_r);
    
    phi_index=1;
	phi_rank=1;
    psi_star_dot_PR_map(:,:)=psi_star_dot_PR_map_phi(phi_rank,:,:);
    E_potential_PR_map(:,:)=E_potential_PR_map_phi(phi_rank,:,:);
%     BstarX_PR_map(:,:)=BstarX_PR_map_phi(phi_rank,:,:);
%     BstarZ_PR_map(:,:)=BstarZ_PR_map_phi(phi_rank,:,:);
      
    % positions phi_rank=1 and phi_rank=NB_PHI have to be treated
    % separately and have same values
    E_potential_PR_map_phi(NB_THETA,:,:)=E_potential_PR_map_phi(1,:,:);
%     BstarX_PR_map_phi(NB_PHI,:,:)=BstarX_PR_map_phi(1,:,:);
%     BstarZ_PR_map_phi(NB_PHI,:,:)=BstarZ_PR_map_phi(1,:,:);
    
    
    E_potential_PR_map_prev(:,:)=E_potential_PR_map_phi(end-1,:,:);
    E_potential_PR_map_next(:,:)=E_potential_PR_map_phi(2,:,:);
    calculate_Efield_phi_rank;
    
    if CALCULATE_EXB_DATA_FILE==1
        calculate_ExB_drift_speeds_phi_rank;
        vExB_X_map_phi(1,:,:)=vExB_X_PR_map;
        vExB_Z_map_phi(1,:,:)=vExB_Z_PR_map;
        vExB_phi_map_phi(1,:,:)=vExB_phi_PR_map;
        vExB_X_map_phi(NB_PHI,:,:)=vExB_X_PR_map;
        vExB_Z_map_phi(NB_PHI,:,:)=vExB_Z_PR_map;
        vExB_phi_map_phi(NB_PHI,:,:)=vExB_phi_PR_map;
    end
    
    Efield_X_map_phi(1,:,:)=EX_PR_map;
    Efield_Z_map_phi(1,:,:)=EZ_PR_map;
    Efield_phi_map_phi(1,:,:)=Ephi_PR_map;
    Efield_X_map_phi(NB_PHI,:,:)=EX_PR_map;
    Efield_Z_map_phi(NB_PHI,:,:)=EZ_PR_map;
    Efield_phi_map_phi(NB_PHI,:,:)=Ephi_PR_map;
    grad_Phi_tor_map_phi(1,:,:)=grad_Phi_tor_PR_map;
    grad_Phi_tor_map_phi(NB_PHI,:,:)=grad_Phi_tor_PR_map;
    
    
    
    
    % Looping around in the toroidal direction
    for (phi_index=2:NB_PHI-1)
        %    for (phi_rank=12:12)
        
        disp('**************************************************************');
        
        E_potential_PR_map=zeros(NP,size_r);
        E_potential_PR_map_prev=zeros(NP,size_r);
        E_potential_PR_map_next=zeros(NP,size_r);
        grad_phi_tor_num=zeros(NP,size_r);
        EX_PR_map=zeros(NP,size_r);
        EZ_PR_map=zeros(NP,size_r);
        Ephi_PR_map=zeros(NP,size_r);
        
        if PHI_OMEGA_RATIO==1
            phi_rank=phi_index
        else
            phi_rank=PHI_OMEGA_RATIO*phi_index-1
        end
        
        psi_star_dot_PR_map(:,:)=psi_star_dot_PR_map_phi(phi_rank,:,:);
        E_potential_PR_map(:,:)=E_potential_PR_map_phi(phi_rank,:,:);
%         BstarX_PR_map(:,:)=BstarX_PR_map_phi(phi_rank,:,:);
%         BstarZ_PR_map(:,:)=BstarZ_PR_map_phi(phi_rank,:,:);
        E_potential_PR_map_prev(:,:)=E_potential_PR_map_phi(phi_rank-1,:,:);
        E_potential_PR_map_next(:,:)=E_potential_PR_map_phi(phi_rank+1,:,:);
        
        calculate_Efield_phi_rank;
        
        if CALCULATE_EXB_DATA_FILE==1
            calculate_ExB_drift_speeds_phi_rank;
            vExB_X_map_phi(phi_index,:,:)=vExB_X_PR_map;
            vExB_Z_map_phi(phi_index,:,:)=vExB_Z_PR_map;
            vExB_phi_map_phi(phi_index,:,:)=vExB_phi_PR_map;
        end
        
        
        
        disp('****************calculate_Efield DONE!**********************');
        
        Efield_X_map_phi(phi_index,:,:)=EX_PR_map;
        Efield_Z_map_phi(phi_index,:,:)=EZ_PR_map;
        Efield_phi_map_phi(phi_index,:,:)=Ephi_PR_map;
        grad_Phi_tor_map_phi(phi_index,:,:)=grad_Phi_tor_PR_map;
        
    end
    
