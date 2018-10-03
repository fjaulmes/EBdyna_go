
close all

REINIT_ALL_TOKAMAK_DATA=1;

if REINIT_ALL_TOKAMAK_DATA==1
    clear all
    initialize_folder_names;
    initialize_collapse_map_calculation_context
    rescaling_to_XZsmall_maps
    DT_INTERPOLATION_METHOD='quadratic'     % by default
    CALCULATE_VD_DATA_FILE=0;
end
CALCULATE_VD_DATA_FILE=1;

% load('../B_maps/B0311.mat');
% load('../E_maps/E0311.mat');
%


frame_rank=32
SAVENAME='vExB_avg_evol.mat'

vExB_radial_evol=zeros(100,NP,size_r);
vExB_helical_evol=zeros(100,NP,size_r);

% load(SAVENAME)

frame_rank=31


% for (frame_rank=2:8:86)
    phi_index=1
    
    
    % for frame_rank=1:101
    disp('----------------------------------------');
    
    reference_frame=min((frame_rank-1)*10+1,size(rx_evol_interp,2))
    
    f=reference_frame;
    
    if (f<10)
        frame_name='00';
    elseif (f<100)
        frame_name='0';
    else
        frame_name='';
    end
    
    frame_name=strcat(frame_name,num2str(f));
    filename_B='../B_maps/B0';
    filename_B=strcat(filename_B,frame_name,'.mat');
    load(filename_B);
    filename_E='../E_maps/E0';
    filename_E=strcat(filename_E,frame_name,'.mat');
    load(filename_E);
    
    % phi_index=12
    vExB_radial_PR_map_phi=zeros(NB_PHI,NP,size_r);
    vExB_helical_PR_map_phi=zeros(NB_PHI,NP,size_r);
    
    for (phi_index=1:NB_PHI)
        
        
        psi_star_omega_map_half(:,:)=psi_star_2D_evol_lin(frame_rank,:,:);
        % Using symmetry to reconstruct a poloidal turn
        psi_star_omega_map_rank=zeros(size_r,NB_THETA);
        psi_star_omega_map_rank(:,1:round(0.5*NB_THETA))=psi_star_omega_map_half(:,:);
        psi_star_omega_map_rank(:,round(0.5*NB_THETA):NB_THETA)=psi_star_omega_map_half(:,round(0.5*NB_THETA):-1:1);
        
        psi_star_omega_map=zeros(size_r,NB_THETA);
        psi_star_PR_map=psi_star_omega_map';
        
        if PHI_OMEGA_RATIO==1
            phi_rank=phi_index
        else
            phi_rank=PHI_OMEGA_RATIO*phi_index-1
        end
        
        % if phi_index==1
        %     phi_rank=1
        % end
        
        psi_star_omega_map(:,:)=rotate_map_phi(psi_star_omega_map_rank,phi_rank);
        psi_star_PR_map=psi_star_omega_map';
        
        delta_phi_rank=phi_rank+1;
        psi_star_omega_map_rank_next(:,:)=rotate_map_phi(psi_star_omega_map_rank,delta_phi_rank);
        psi_star_PR_map_rank_next=psi_star_omega_map_rank_next';
        
        delta_phi_rank=phi_rank-1;
        psi_star_omega_map_rank_prev(:,:)=rotate_map_phi(psi_star_omega_map_rank,delta_phi_rank);
        psi_star_PR_map_rank_prev=psi_star_omega_map_rank_prev';
        
        calculate_Btot_phi_rank;
        
        
        Epot_omega_map(:,:)=Epot_evol(frame_rank,:,:);
        
        E_potential_PR_map_phi=zeros(NB_THETA,NB_THETA,size_r);
        for index=1:NB_THETA
            %         phi_rank=round(2*phi_index-1);
            E_potential_PR_map_phi(index,:,:)=rotate_map_phi(Epot_omega_map,index)';
        end
        
        psi_star_dot_omega_map(:,:)=psi_star_dot_evol(frame_rank,:,:);
        
        psi_star_dot_PR_map_phi=zeros(NB_THETA,NB_THETA,size_r);
        for index=1:NB_THETA
            %         phi_rank=round(2*phi_index-1);
            psi_star_dot_PR_map_phi(index,:,:)=rotate_map_phi(psi_star_dot_omega_map,index)';
        end
        
        psi_star_dot_PR_map(:,:)=psi_star_dot_PR_map_phi(phi_rank,:,:);
        % load psi_star_dot_PR_map_ref3.mat
        % psi_star_dot_PR_map(:,:)=psi_star_dot_PR_map_ref3(:,1:size_r);
        
        %     E_data=reshape(psi_star_dot_PR_map(:,1:size_r),NP*size_r,1);
        %     % Efield_3_XZ_map=griddata(finesse_data_X(1:size_r*NP),finesse_data_Z(1:size_r*NP),E_data,XX_small,ZZ_small,'cubic');
        %     Efield_3_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,E_data(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
        %     Efield_3_XZ_map=Efield_3_XZ_map';
        %     Efield_3_XZ_map(isnan(Efield_3_XZ_map))=0;
        
        %     Ephi_XZ_map=Efield_3_XZ_map./Rpos_XZsmall_map;
        % %
        % Bfield_3_XZmap=Bphi_XZsmall_map./Rpos_XZsmall_map;
        %
        E_potential_PR_map(:,:)=E_potential_PR_map_phi(phi_rank,:,:);
        if phi_rank==1
            E_potential_PR_map_prev(:,:)=E_potential_PR_map_phi(end-1,:,:);
            E_potential_PR_map_next(:,:)=E_potential_PR_map_phi(phi_rank+1,:,:);
        elseif phi_rank==257
            E_potential_PR_map_prev(:,:)=E_potential_PR_map_phi(phi_rank-1,:,:);
            E_potential_PR_map_next(:,:)=E_potential_PR_map_phi(2,:,:);
        else
            E_potential_PR_map_prev(:,:)=E_potential_PR_map_phi(phi_rank-1,:,:);
            E_potential_PR_map_next(:,:)=E_potential_PR_map_phi(phi_rank+1,:,:);
        end
        
        %     calculate_Efield_phi_rank;
        Epot_data=reshape(E_potential_PR_map(:,1:size_r),NP*size_r,1);
        
        Epot_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,Epot_data(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
        Epot_XZ_map=Epot_XZ_map';
        Epot_XZ_map(isnan(Epot_XZ_map))=0;
        
        calculate_Efield_XZ;
        
        E_potential_PR_map_DPHI=E_potential_PR_map_next-E_potential_PR_map_prev;
        grad_Phi_tor_PR_map=(E_potential_PR_map_DPHI./Rpos_PR_map(:,1:size_r))/(2*DOMEGA);
        Ephi_PR_map=psi_star_dot_PR_map./Rpos_PR_map(:,1:size_r)-grad_Phi_tor_PR_map;
        
        EX_map(:,:)=Efield_X_map_phi(phi_index,:,:);
        EZ_map(:,:)=Efield_Z_map_phi(phi_index,:,:);
        
        
        %     vExB_X_map=EZ_XZ_map.*Bphi_XZsmall_map-Ephi_XZ_map.*BpolZ_XZ_map;
        %     vExB_Z_map=Ephi_XZ_map.*BpolX_XZ_map-Bphi_XZsmall_map.*EX_XZ_map;
        %     vExB_phi_map=EX_XZ_map.*BpolZ_XZ_map-BpolX_XZ_map.*EZ_XZ_map;
        %     vExB_X_map=vExB_X_map./(Btot_XZ_map.^2);
        %     vExB_Z_map=vExB_Z_map./(Btot_XZ_map.^2);
        %     vExB_phi_map=vExB_phi_map./(Btot_XZ_map.^2);
        
        
        % dr_data=reshape(dr_X_PR_map(:,1:size_r),NP*size_r,1);
        % drX_map=dtinterp(finesse_mesh,finesse_mesh_dtri,dr_data(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
        % drX_map=drX_map';
        % dr_data=reshape(dr_Z_PR_map(:,1:size_r),NP*size_r,1);
        % drZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,dr_data(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
        % drZ_map=drZ_map';
        %
        % vExB_radial=vExB_X_map.*drX_map+vExB_Z_map.*drZ_map;
        
        %     Btot_PR_map=sqrt(BpolX_PR_map.^2+BpolZ_PR_map.^2+Btor_PR_map(:,1:size_r).^2);
        
        vExB_X_PR_map=EZ_PR_map.*Btor_PR_map(:,1:size_r)-Ephi_PR_map.*BpolZ_PR_map;
        vExB_Z_PR_map=Ephi_PR_map.*BpolX_PR_map-Btor_PR_map(:,1:size_r).*EX_PR_map;
        vExB_phi_PR_map=EX_PR_map.*BpolZ_PR_map-BpolX_PR_map.*EZ_PR_map;
        vExB_X_PR_map=vExB_X_PR_map./(Btot_PR_map(:,1:size_r).^2);
        vExB_Z_PR_map=vExB_Z_PR_map./(Btot_PR_map(:,1:size_r).^2);
        vExB_phi_PR_map=vExB_phi_PR_map./(Btot_PR_map(:,1:size_r).^2);
        
        vExB_radial_PR_map=vExB_X_PR_map.*dr_X_PR_map(:,1:size_r)+vExB_Z_PR_map.*dr_Z_PR_map(:,1:size_r);
        vExB_poloidal_PR_map=vExB_X_PR_map.*grad_theta_PR_map_X(:,1:size_r)+vExB_Z_PR_map.*grad_theta_PR_map_Z(:,1:size_r);
        vExB_toroidal_PR_map=vExB_phi_PR_map;
        
        vExB_radial_PR_map_phi(phi_index,:,:)=vExB_radial_PR_map;
        vExB_helical_PR_map_phi(phi_index,:,:)=vExB_poloidal_PR_map./sqrt(X_PR_map(:,1:size_r).^2+Z_PR_map(:,1:size_r).^2)-vExB_toroidal_PR_map./Rpos_PR_map(:,1:size_r);
        
    end
    
    vExB_radial_omegaR_map_phi=vExB_radial_PR_map_phi;
    vExB_helical_omegaR_map_phi=vExB_helical_PR_map_phi;
    
    for (phi_index=2:NB_PHI-1)
        if PHI_OMEGA_RATIO==1
            phi_rank=phi_index;
        else
            phi_rank=PHI_OMEGA_RATIO*phi_index-1;
        end
        vExB_radial_omegaR_map_phi(phi_index,:,:)=[squeeze(vExB_radial_PR_map_phi(phi_index,phi_rank:end,:)) ; squeeze(vExB_radial_PR_map_phi(phi_index,2:phi_rank,:))];
        vExB_helical_omegaR_map_phi(phi_index,:,:)=[squeeze(vExB_helical_PR_map_phi(phi_index,phi_rank:end,:)) ; squeeze(vExB_helical_PR_map_phi(phi_index,2:phi_rank,:))];
    end
    
    vExB_radial_PR_map_avg=squeeze(mean(vExB_radial_omegaR_map_phi,1));
    vExB_helical_PR_map_avg=squeeze(mean(vExB_helical_omegaR_map_phi,1));
    
    vExB_radial_evol(frame_rank,:,:)=vExB_radial_PR_map_avg;
    vExB_helical_evol(frame_rank,:,:)=vExB_helical_PR_map_avg;
    
% end

save(SAVENAME,'vExB_radial_evol','vExB_helical_evol')



figure(1)
subplot(2,1,1)
set(gca,'FontSize',26);
imagesc(radial_r_value_flux(1:size_r),((0:NP-1)/(NP-1))*2*pi , vExB_radial_PR_map_avg);
xlim([0 0.44])
xlabel('r')
ylabel('\omega')
% colorbar


subplot(2,1,2)
set(gca,'FontSize',26);
hold on
grid on
plot(radial_r_value_flux(1:size_r),mean(vExB_radial_PR_map_avg,1))
xlim([0 0.44])
xlabel('r')
ylabel('<v_E_r> (m/s)')
% colorbar


figure(2)
set(gca,'FontSize',26);
imagesc(radial_r_value_flux(1:size_r),((0:NP-1)/(NP-1))*2*pi , abs(vExB_helical_PR_map_avg));
xlim([0 0.44])
xlabel('r')
ylabel('\omega')
colorbar


figure(3)
subplot(2,1,1)
set(gca,'FontSize',26);
imagesc(((0:NP-1)/(NP-1))*2*pi, radial_r_value_flux(1:size_r), (vExB_helical_PR_map_avg)');
ylim([0 0.44])
xlim([0 2*pi])
xlabel('\omega')
ylabel('r')
colorbar


subplot(2,1,2)
set(gca,'FontSize',26);
hold on
grid on
dr_hel1=round(0.3*size_r);
plot(((0:NP-1)/(NP-1))*2*pi,mean(vExB_helical_PR_map_avg(:,end-dr_hel1:end),2),'b','linewidth',2)
dr_hel2=round(0.5*size_r);
plot(((0:NP-1)/(NP-1))*2*pi,mean(vExB_helical_PR_map_avg(:,end-dr_hel2:end),2),'k--','linewidth',3)
dr_hel3=round(0.7*size_r);
plot(((0:NP-1)/(NP-1))*2*pi,mean(vExB_helical_PR_map_avg(:,end-dr_hel3:end),2),'r-.','linewidth',4)
width1=radial_r_value_flux(size_r)-radial_r_value_flux(size_r-dr_hel1)
width2=radial_r_value_flux(size_r)-radial_r_value_flux(size_r-dr_hel2)
width3=radial_r_value_flux(size_r)-radial_r_value_flux(size_r-dr_hel3)
legend(strcat('\delta_r=',num2str(width1,3)),strcat('\delta_r=',num2str(width2,3)),strcat('\delta_r=',num2str(width3,3)));
xlim([0 2*pi])
xlabel('\omega')
ylabel('\omega_{v_E} (rad/s)')
colorbar
