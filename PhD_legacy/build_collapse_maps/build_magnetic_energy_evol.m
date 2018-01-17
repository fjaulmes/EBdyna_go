% reset_data_analysis_environment
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

PHI_OMEGA_RATIO=64;
NB_PHI=round((NP-1)/PHI_OMEGA_RATIO);
DPHI=(2*pi/(NB_PHI));
DOMEGA=DPHI/PHI_OMEGA_RATIO;

CALCULATE_VD_DATA_FILE=1;
CALCULATE_PR_DATA_FILE=0;
CALCULATE_BSTAR_DATA=0;

TIME_RESOLUTION=2


 load('../data_tokamak/volume_flux_geometry.mat');


SAVENAME='AUG30382_2p9_mag_energy_evol.mat'
psi_max=size_r-1
psi_max_value=interp1(1:Nradial,psi_scale,psi_max+2)
psi_data=reshape(psi_PR_map(:,1:NR),NB_THETA*NR,1);
% psi_XZ_map=griddata(finesse_data_X_extended,finesse_data_Z_extended,psi_data,XX_small,ZZ_small,'cubic');
psi_XZ_map_ini=dtinterp(finesse_mesh_extended,finesse_mesh_extended_dtri,psi_PR_map(IDTRI_EXT),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
psi_XZ_map_ini(isnan(psi_XZ_map_ini))=0;
psi_XZ_map_ini=psi_XZ_map_ini';
psi_XZ_map_mask=psi_XZ_map_ini*0;

for (x=3:sizeX-2)
    for (z=3:sizeZ-2)
        if psi_XZ_map_ini(x,z)>psi_max_value
            psi_XZ_map_mask(x,z)=1;
        end
    end
end

Bphi_XZsmall_map=Bphi_XZsmall_map.*psi_XZ_map_mask;

% B2_tot_evol=zeros(100,1);
% energy_phi_evol=zeros(100,NB_PHI);

% load(SAVENAME)

DPHI=(2*pi/(NB_PHI));

% for (frame_rank=2:8:86)
    phi_index=1
   
    

    
for frame_rank=TIME_RESOLUTION:TIME_RESOLUTION:100
    disp('----------------------------------------');
    %%
    reference_frame=min((frame_rank-1)*10+1,1001)
    
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
            phi_rank=PHI_OMEGA_RATIO*(phi_index-1)+1
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
        
        energy_mag=0;
        energy_star_mag=0;
        for (x=Xinf+1:Xsup-1)
            Z_inf=round(interp1(scale_Z,1:length(scale_Z),Z_psi_fit_down(psi_max,x)));
            Z_sup=round(interp1(scale_Z,1:length(scale_Z),Z_psi_fit_up(psi_max,x)));
            for (z=Z_inf:Z_sup)
                energy_mag=energy_mag+0.5*(Btot_XZ_map(x-Xinf,z)^2)*(scale_X(x-Xinf)+R0)*(DPHI*DX*DX/mu0);
                if CALCULATE_BSTAR_DATA==1
                    energy_star_mag=energy_star_mag+0.5*(Bstartot_XZ_map(x-Xinf,z)^2)*(scale_X(x-Xinf)+R0)*(DPHI*DX*DX/mu0);
                end
            end
        end
       
        energy_phi_evol(round(frame_rank/TIME_RESOLUTION),phi_index)=energy_mag;
        if CALCULATE_BSTAR_DATA==1
            energy_star_phi_evol(round(frame_rank/TIME_RESOLUTION),phi_index)=energy_star_mag;
        end
        
        
    end
    
%     figure(3)
%     imagesc(psi_star_PR_map',[-18 1]*1e-3);
%     colorbar
%     pause(0.1)
%     energy_mag
    
    energy_frame_rank=sum(energy_phi_evol(round(frame_rank/TIME_RESOLUTION),1:end-1),2)
    B2_tot_evol(round(frame_rank/TIME_RESOLUTION))=energy_frame_rank;
    if CALCULATE_BSTAR_DATA==1
        energy_star_frame_rank=sum(energy_star_phi_evol(round(frame_rank/TIME_RESOLUTION),1:end-1),2)
        Bstar2_tot_evol(round(frame_rank/TIME_RESOLUTION))=energy_star_frame_rank;
    end
end   
    

%save(SAVENAME,'Bstar2_tot_evol','energy_star_phi_evol','B2_tot_evol','energy_phi_evol','time_scale_lin')

if CALCULATE_BSTAR_DATA==1
    save(SAVENAME,'Bstar2_tot_evol','energy_star_phi_evol','B2_tot_evol','energy_phi_evol','time_scale_lin')
else
    save(SAVENAME,'B2_tot_evol','energy_phi_evol','time_scale_lin')
end

%%
figure(2)
grid on
hold on
set(gca,'FontSize',26);
plot(time_scale_lin(TIME_RESOLUTION:TIME_RESOLUTION:end),B2_tot_evol(1:1:end)-B2_tot_evol(end-1),'b','LineWidth',2)
xlim([0.0 0.9])
ylabel('J')
xlabel('t/\tau_{coll}')

%%
if CALCULATE_BSTAR_DATA==1
    
    figure(3)
    grid on
    hold on
    set(gca,'FontSize',26);
    plot(time_scale_lin(TIME_RESOLUTION:TIME_RESOLUTION:87),Bstar2_tot_evol,'b','LineWidth',2)
    xlim([0.0 0.9])
    ylabel('J')
    xlabel('t/\tau_{coll}')
end