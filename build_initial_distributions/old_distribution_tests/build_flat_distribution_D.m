clear all;
initialize_folder_names;
filename=strcat(DATA_FOLDER,'physics_constants.mat');
load(filename);
mr=mDT*mHe/(mDT+mHe);
filename=strcat(DATA_FOLDER,'pressure_profile.mat');
load(filename);
filename=strcat(DATA_FOLDER,'motions_map_dimensions.mat');
load(filename);

%pause
NB_PROCESS=16


for n=1:NB_PROCESS
    close all
    
    % Main important values for description of the distribution
    PROCESS_NUMBER=n
    PROCESS_VAL=(round(PROCESS_NUMBER)-0.5)*10
    w0=PROCESS_VAL*(10^3)*eV;
    v0=sqrt(2*w0/mHe);
    
    
    delta_psi=14;
    Npart_min=450;
    
    mHe=mD
    ZHe=1
    
    
    % vpll distribution
    N_vparallel_bins=31;
    
    N_vparallel_bins_half=round(0.5*(N_vparallel_bins-1));
    
    vpll_range=(v0)*(-N_vparallel_bins_half:N_vparallel_bins_half)*(1/N_vparallel_bins_half);
    vpll_bin_size=v0/N_vparallel_bins_half
    
    vpll_range_inf=(v0)*(-N_vparallel_bins_half:N_vparallel_bins_half)*(1/N_vparallel_bins_half)-0.5*vpll_bin_size;
    vpll_range_sup=(v0)*(-N_vparallel_bins_half:N_vparallel_bins_half)*(1/N_vparallel_bins_half)+0.5*vpll_bin_size;
    
    energy_bin_size=0.01*w0;
    
    w0_inf=(w0-0.5*energy_bin_size)*eV;
    w0_sup=(w0+0.5*energy_bin_size)*eV;
    % theta0=sqrt(2*w0/mHe);
    % theta0_inf=sqrt(2*w0_inf/mHe);
    % theta0_sup=sqrt(2*w0_sup/mHe);
    
    
    
    %Normalize the alpha density
    Nalpha_binned=1.2e4;
    Nalpha_max=1.2e4;
    Nalpha_binned_norm=Nalpha_binned;
    
    
    % this number need to be distributed
    % according to the density of particles
    
    filename=strcat(DATA_FOLDER,'volume_flux_geometry.mat');
    load(filename);
    filename=strcat(DATA_FOLDER,'XZsmall_fields_tokamak_pre_collapse.mat');
    load(filename);
    
    
    % rescaling flux surface curves on small map
    psi_rank=1;
    psi_rank_max=size_r+delta_psi
    
    X_min=mid_Xaxis_large-0.5*size(scale_X,2);
    X_max=mid_Xaxis_large+0.5*size(scale_X,2)-1;
    X_bin_size=24;
    X_bin_half_std=round(0.5*(X_bin_size));
    
    X1_Nradial=X1_Nradial-mid_Xaxis_large+mid_X;
    X2_Nradial=X2_Nradial-mid_Xaxis_large+mid_X;
    Z_psi_fit_up_small=Z_psi_fit_up(:,X_min:X_max);
    Z_psi_fit_down_small=Z_psi_fit_down(:,X_min:X_max);
    
    disp('Adjust here in case part of the distribution is missing......')
    
    % X_min=X1_Nradial(size_r+3*X_bin_size);
    % X_max=X2_Nradial(size_r+3*X_bin_size);
    
    
    X_offset=-round(0.5*mod(X_max-X_min,X_bin_size))+1+X_bin_half_std;
    if mod(X_max-X_min,X_bin_size)==0
        N_X_bins=floor((X_max-X_min)/X_bin_size);
    else
        N_X_bins=floor((X_max-X_min)/X_bin_size)+1;
    end
    N_X_bins=N_X_bins-2
    
    X_bin_pos=8+X_offset+X_bin_size*(0:N_X_bins-1)
    X_bin_pos_remap=X_bin_pos+mid_Xaxis_large-mid_X;
    
    N_radial_bins=1
    
    % volume of vertical sclices of toroidal dougnnuts shaped bins
    % (yummy)
    
    volume_bin=zeros(1,N_X_bins);
    Npart_binned=zeros(1,N_X_bins);
    
    
    for X_bin_rank=1:N_X_bins
        Xinf=min(X_bin_pos_remap(X_bin_rank)-X_bin_half_std,size(volume_tor_diff,2));
        Xsup=min(X_bin_pos_remap(X_bin_rank)+X_bin_half_std-1,size(volume_tor_diff,2));
        volume_bin(psi_rank,X_bin_rank)=sum(sum(volume_tor_diff(1:psi_rank_max,Xinf:Xsup)));
    end
    
    
    
    volume_bin_max=max(max(volume_bin));
    
    
    Npart_binned(1,:)=round(volume_bin(1,:)*Nalpha_binned_norm);
    
    
    % there should be no less than 2 particles for each energy in a volume
    Npart_binned=round(Npart_min*Npart_binned/min(Npart_binned(Npart_binned~=0)));
    
    DX_MIN_BIN=0;
    % compensating for the edge bins
    
    % if (X_rank_inf(psi_rank))<(X_rank_sup(psi_rank)-2)
    %     for X_bin_rank=1:N_X_bins
    %         if (X_rank_inf(psi_rank)==X_bin_rank)
    %             DX_MIN_BIN=-1;
    %             Npart_binned(psi_rank,X_bin_rank-DX_MIN_BIN)=Npart_binned(psi_rank,X_bin_rank-DX_MIN_BIN)+Npart_binned(psi_rank,X_bin_rank);
    %             Npart_binned(psi_rank,X_bin_rank)=0;
    %         end
    %         if (X_rank_sup(psi_rank)==X_bin_rank)
    %             DX_MIN_BIN=1;
    %             Npart_binned(psi_rank,X_bin_rank-DX_MIN_BIN)=Npart_binned(psi_rank,X_bin_rank-DX_MIN_BIN)+Npart_binned(psi_rank,X_bin_rank);
    %             Npart_binned(psi_rank,X_bin_rank)=0;
    %         end
    %     end
    % end
    
    disp('smallest number of particles');
    Npart_min=min(Npart_binned(Npart_binned~=0))
    
    % this gives us the ratio of particles to density
    density_part_ratio=(4e15)/sum(Npart_binned)
    
    Nalphas_simulated=sum(sum(Npart_binned));
    
    disp('number of particles generated');
    Nalphas_simulated
    
    alphas_pos_x=zeros(Nalphas_simulated,1);
    alphas_pos_z=zeros(Nalphas_simulated,1);
    alphas_pos_phi=zeros(Nalphas_simulated,1);
    alphas_Ekin=zeros(Nalphas_simulated,1);
    alphas_mm=zeros(Nalphas_simulated,1);
    alphas_vpll=zeros(Nalphas_simulated,1);
    
    Npart_rank_begin=1;
    Npart_rank_flux_surface_begin=1;
    
    
    
    %Z_bin_vector=zeros(N_radial_bins*N_X_bins*2,1);
    Z_index=1;
    
    psi_rank=1;
    
    psi0=size_r;
    X_min=X1_Nradial(psi0+delta_psi);
    X_max=X2_Nradial(psi0+delta_psi);
    
    Z_sup_max=Z_psi_fit_up_small(psi0+delta_psi,:);
    Z_sup_min=Z_psi_fit_up_small(1,:);
    
    Z_inf_min=Z_psi_fit_down_small(1,:);
    Z_inf_max=Z_psi_fit_down_small(psi0+delta_psi,:);
    
    Npart_flux_surface(psi_rank)=sum(Npart_binned(psi_rank,:));
    Npart_rank_flux_surface_end=Npart_rank_flux_surface_begin+sum(Npart_binned(psi_rank,:))-1;
    disp('Number of particles distributed in volume of flux surface element:');
    disp(Npart_flux_surface(psi_rank));
    
    
    % display particles
    figure(1);
    set(gca,'FontSize',22);
    
    hold on;
    axis xy square
    
    for X_bin_rank=1:N_X_bins
        
        
        X0=X_bin_pos(X_bin_rank);
        X_bin_half=X_bin_half_std;
        
        Z_sup_max_value=max(Z_sup_max(X0-X_bin_half-1:X0+X_bin_half+1));
        Z_sup_min_value=min(Z_sup_min(X0-X_bin_half-1:X0+X_bin_half+1));
        Z_sup_pos=0.5*(Z_sup_max_value+Z_sup_min_value);
        
        Z_inf_max_value=max((Z_inf_min(X0-X_bin_half-1:X0+X_bin_half+1)));
        Z_inf_min_value=min((Z_inf_max(X0-X_bin_half-1:X0+X_bin_half+1)));
        Z_inf_pos=-0.5*(Z_inf_max_value+Z_inf_min_value);
        
        x_min=interp1(1:size(scale_X,2),scale_X,X0-X_bin_half,'*linear');
        x_max=interp1(1:size(scale_X,2),scale_X,X0+X_bin_half,'*linear');
        
        z_min=1;
        z_max=-1;
        z=0;
        Npart_rank_end=min(Npart_rank_begin+Npart_binned(psi_rank,X_bin_rank),Nalphas_simulated);
        
        if (Npart_rank_begin<Npart_rank_end)
            Npart_binned(1,X_bin_rank)
            
            for(n=Npart_rank_begin:Npart_rank_end)
                % random coordinate on [x_min ; x_max[
                while (z<z_min)||(z>=z_max)
                    x=my_rand(1)*(x_max-x_min)+x_min;
                    z_max=interp1(scale_X(X0-X_bin_half-1:X0+X_bin_half+1),Z_sup_max(X0-X_bin_half-1:X0+X_bin_half+1),x,'*linear');
                    z_min=interp1(scale_X(X0-X_bin_half-1:X0+X_bin_half+1),Z_inf_max(X0-X_bin_half-1:X0+X_bin_half+1),x,'*linear');
                    z=(my_rand(1)*(Z_sup_max_value-Z_inf_min_value)+Z_inf_min_value);
                end
                alphas_pos_x(n)=x;
                alphas_pos_z(n)=z;
                
                z_min=1;
                z_max=-1;
                z=0;
            end
            
            figure(1);
            plot(alphas_pos_x(Npart_rank_begin:Npart_rank_end),alphas_pos_z(Npart_rank_begin:Npart_rank_end),'.');
            pause(0.05);
            
        end
        
        Npart_rank_begin=Npart_rank_begin+Npart_binned(psi_rank,X_bin_rank);
        Z_volume_bin(psi_rank,X_bin_rank,1)=Z_sup_pos;
        Z_volume_bin(psi_rank,X_bin_rank,2)=Z_inf_pos;
        Z_index=Z_index+2;
        
    end
    
    DNpart=Nalphas_simulated;
    
    Nrank=1;
    Epart_vector=ones(DNpart,1)*w0/eV;
    
    % distribution in v_parallel range
    
    %vpll_vector=ones(1,1);
    
    n_vpll_count=0;
    Nrank=1;
    %vpll_vector=ones(DNpart,1);
    vpll_vector=ones(1,1);
    vpll_max=sqrt(2*(eV/mHe)*Epart_vector);
    
    
    vpll_count=1;
    Nrank_Ebegin=Nrank;
    
    
    max_vpll_Erank=vpll_max;
    vpll_vector=rand(Nalphas_simulated,1).*(2*max_vpll_Erank)-max_vpll_Erank;
    
    
    % happily shuffle everything
    permutation_vector=randperm(DNpart);
    Epart_vector=Epart_vector(permutation_vector);
    vpll_vector=vpll_vector(permutation_vector);
    
    alphas_Ekin(Npart_rank_flux_surface_begin:Npart_rank_flux_surface_end)=Epart_vector;
    alphas_vpll(Npart_rank_flux_surface_begin:Npart_rank_flux_surface_end)=vpll_vector;
    
    Npart_rank_flux_surface_begin=Npart_rank_flux_surface_begin+sum(Npart_binned(psi_rank,:));
    
    figure(1);
    plot(scale_X,Z_sup_max,'r');
    plot(scale_X,Z_sup_min,'r');
    plot(scale_X,Z_inf_max,'r');
    plot(scale_X,Z_inf_min,'r');
    pause(0.1);
    
    
    
    % filling up the uniformly distributes toroidal position values
    for(n=1:Nalphas_simulated)
        alphas_pos_phi(n)=my_rand(1)*2*pi;
    end
    
    alphas_mm=(alphas_Ekin-0.5*(mHe*alphas_vpll.^2)/eV)./interp2(scale_X,scale_Z,Btot_XZ_map',alphas_pos_x,alphas_pos_z,'*linear');
    
    
    
    FILENAME=strcat('initial_flatD_',num2str(PROCESS_VAL),'keV_distribution.mat')
    save (FILENAME,'alphas_pos_x','alphas_pos_z','alphas_pos_phi','alphas_Ekin','alphas_mm','alphas_vpll','Nalphas_simulated','X_bin_pos','X_bin_size','volume_bin','mHe','ZHe');
    
    
    Nalphas_simulated
    
end
