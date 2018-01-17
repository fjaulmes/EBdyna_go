filename=strcat(DATA_FOLDER,'volume_flux_geometry.mat');
load(filename);

psi_bin_size=16;
psi_bin_size_half=0.5*psi_bin_size
max_psi_binned=round((size_r+10)/psi_bin_size)*psi_bin_size+1+psi_bin_size_half;
psi_bin_pos=(psi_bin_size_half+1:psi_bin_size:max_psi_binned);
N_psi_bins=length(psi_bin_pos)
psi_sup_binned=max_psi_binned+psi_bin_size_half-1

X_bin_size=10;
X_bin_size_half=0.5*X_bin_size

X_inf=X1_Nradial(max_psi_binned)+X_bin_size_half
X_max=X2_Nradial(max_psi_binned)-X_bin_size_half;
N_X_bins=round((X_max-X_inf)/X_bin_size)+1
X_max=X_inf+N_X_bins*X_bin_size;
if X_max<X2_Nradial(max_psi_binned)-X_bin_size_half
    disp('problem with the LFS boundary......')
end

% check if this does not make too small bins
X_bin_pos=(X_inf:X_bin_size:X_max);

% now calculate the approximate volumes of all thes bins
clear volume_bin_avg_up volume_bin_avg_down
for(p=1:N_psi_bins)
    for(x=1:N_X_bins)
        psi=psi_bin_pos(p);
        X_pos=X_bin_pos(x);
        % for the moment, take same volumes on both 
        % sides of horizontal axis....
        % need to be improved later ?
        volume_bin_avg_up(p,x)=0.5*sum(sum(volume_tor_diff(psi-psi_bin_size_half:psi+psi_bin_size_half-1, X_pos-X_bin_size_half:X_pos+X_bin_size_half-1)));
        volume_bin_avg_down(p,x)=0.5*sum(sum(volume_tor_diff(psi-psi_bin_size_half:psi+psi_bin_size_half-1, X_pos-X_bin_size_half:X_pos+X_bin_size_half-1)));
    end
end

alphas_Xrank=alphas_pos_x/DX+801;
Nb_part=0;
clear N_volume_up_binned N_volume_down_binned
% counting the number of particles in each bin
for(x=1:N_X_bins)
    X_pos=X_bin_pos(x);
    particles_in_bin_up=find((alphas_Xrank>=X_pos-X_bin_size_half).*(alphas_Xrank<X_pos+X_bin_size_half).*(alphas_pos_z>=Z_axis));
    particles_in_bin_down=find((alphas_Xrank>=X_pos-X_bin_size_half).*(alphas_Xrank<X_pos+X_bin_size_half).*(alphas_pos_z<Z_axis));
    for(p=1:N_psi_bins)
        psi=psi_bin_pos(p);
        Z_up_sup=Z_psi_fit_up(psi+psi_bin_size_half-1,X_pos-X_bin_size_half:X_pos+X_bin_size_half);
        Z_up_inf=Z_psi_fit_up(psi-psi_bin_size_half,X_pos-X_bin_size_half:X_pos+X_bin_size_half);
        Z_down_sup=Z_psi_fit_down(psi-psi_bin_size_half,X_pos-X_bin_size_half:X_pos+X_bin_size_half);
        Z_down_inf=Z_psi_fit_down(psi+psi_bin_size_half,X_pos-X_bin_size_half:X_pos+X_bin_size_half);
        %upper bin
        clear Z_up_sup_part Z_up_inf_part
        Z_up_sup_part=interp1(X_pos-X_bin_size_half:X_pos+X_bin_size_half,Z_up_sup,alphas_Xrank(particles_in_bin_up),'*linear');
        Z_up_inf_part=interp1(X_pos-X_bin_size_half:X_pos+X_bin_size_half,Z_up_inf,alphas_Xrank(particles_in_bin_up),'*linear');
        PARTS_BINNED=find((Z_up_sup_part-alphas_pos_z(particles_in_bin_up)>0).*(alphas_pos_z(particles_in_bin_up)-Z_up_inf_part>=0));
        PARTS_BINNED=particles_in_bin_up(PARTS_BINNED);
        Nb_part=Nb_part+length(PARTS_BINNED);
        N_volume_up_binned(p,x)=length(PARTS_BINNED);
        %lower bin
        clear Z_down_sup_part Z_down_inf_part
        Z_down_sup_part=interp1(X_pos-X_bin_size_half:X_pos+X_bin_size_half,Z_down_sup,alphas_Xrank(particles_in_bin_down),'*linear');
        Z_down_inf_part=interp1(X_pos-X_bin_size_half:X_pos+X_bin_size_half,Z_down_inf,alphas_Xrank(particles_in_bin_down),'*linear');
        PARTS_BINNED=find((Z_down_sup_part-alphas_pos_z(particles_in_bin_down)>0).*(alphas_pos_z(particles_in_bin_down)-Z_down_inf_part>=0));
        PARTS_BINNED=particles_in_bin_down(PARTS_BINNED);
        Nb_part=Nb_part+length(PARTS_BINNED);
        N_volume_down_binned(p,x)=length(PARTS_BINNED);
    end
end

density_binned_up=N_volume_up_binned./volume_bin_avg_up;
density_binned_down=N_volume_down_binned./volume_bin_avg_down;

density_binned_up(isnan(density_binned_up))=0;
density_binned_down(isnan(density_binned_down))=0;
density_binned_up(isinf(density_binned_up))=0;
density_binned_down(isinf(density_binned_down))=0;

clear density_x_profile
density_psi_x=(density_binned_up+density_binned_down);
for(p=1:N_psi_bins)
    density_psi_profile=density_psi_x(p,:);
    density_x_value=mean(density_psi_profile(density_psi_profile~=0));
    density_x_profile(p)=density_x_value;
end

figure(1)
set(gca,'FontSize',22);
hold on;
grid on;
plot(density_x_profile);axis xy

