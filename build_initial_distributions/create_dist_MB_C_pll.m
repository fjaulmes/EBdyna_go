clear all;
build_dist_MB_C;


filename=strcat(DATA_FOLDER,'plasma_current.mat');
load(filename, 'JPHI_XZ_map')

% 
% % co field rotation (co-current)
% mock up of rotation profile based on pressure
% ROTATION_profile_radial_ini=CORE_ROTATION*P_initial_profile/P_initial_profile(1);
rotation_profile=ROTATION_profile_radial_ini;

filename=strcat(DATA_FOLDER,'tokamak_PR_map.mat');
load(filename);
filename=strcat(DATA_FOLDER,'flux_geometry.mat');
load(filename);
filename=strcat(DATA_FOLDER,'volume_flux_geometry.mat');
load(filename);
filename=strcat(DATA_FOLDER,'XZsmall_fields_tokamak_pre_collapse.mat');
load(filename);


%recalculate volume elements
volume_flux=zeros(1,513);
volume_flux(2:end)=cumsum(sum(dVol,1));

Npart_simulated=1.9e6;

%pause
close all


%Normalize the alpha density
Nalpha_max=max(Nvalues);
Nalpha_binned_norm=Nvalues/Nalpha_max;


% this number need to be distributed
% according to the density of particles 


NTHETA=size(X_PR_map,2);

theta_scale=2*pi*((1:NTHETA)-1)/(NTHETA-1);

% rescaling flux surface curves on small map
X_min=mid_Xaxis_large-0.5*size(scale_X,2);
X_max=mid_Xaxis_large+0.5*size(scale_X,2)-1;

X1_Nradial=X1_Nradial-mid_Xaxis_large+mid_X;
X2_Nradial=X2_Nradial-mid_Xaxis_large+mid_X;
Z_psi_fit_up_small=Z_psi_fit_up(:,X_min:X_max);
Z_psi_fit_down_small=Z_psi_fit_down(:,X_min:X_max);

theta_bin_size=8;
theta_bin_pos=(1:theta_bin_size:NTHETA);
theta_bin_pos_values=((theta_bin_pos-1)+0.5*theta_bin_size);
theta_bin_pos_values=theta_bin_pos_values(1:end-1)
theta_bin_values=2*pi*(theta_bin_pos_values/(NTHETA-1));
theta_bin_values=theta_bin_values(1:end-1);
N_theta_bins=length(theta_bin_values);
theta_bin_bounds=2*pi*((1:N_theta_bins+1)-1)/(N_theta_bins);


psi_bin_pos=radial_bins
volume_bin_radial_edges=interp1(1:Nradial,volume_flux,radial_bins_edges);
volume_bin_radial=volume_bin_radial_edges(2:end)-volume_bin_radial_edges(1:end-1);




    close all
%     clear  alphas_pos_x alphas_pos_z alphas_pos_phi alphas_Ekin alphas_mm alphas_vpll

    volume_bin=zeros(N_radial_bins,N_theta_bins);
    Npart_binned=zeros(N_radial_bins,N_theta_bins);
    Npart_theta_binned=zeros(N_radial_bins,N_theta_bins);
    Fpart_theta_binned=zeros(N_radial_bins,N_theta_bins+1);
    
    for psi_rank=psi_pos_inf:psi_pos_sup
        for theta_rank=1:N_theta_bins
            psi_pos_avg=psi_bin_pos(psi_rank);
            theta_pos_inf=theta_bin_pos(theta_rank);
            theta_pos_sup=theta_bin_pos(theta_rank+1)-1;
            volume_bin(psi_rank,theta_rank)=sum(sum(dVol(theta_pos_inf:theta_pos_sup,psi_pos_avg-radial_bin_half_size:psi_pos_avg+radial_bin_half_size-1),1),2);
        end
    end
    
    for psi_rank=psi_pos_inf:psi_pos_sup
        
        psi_pos_avg=psi_bin_pos(psi_rank);
        Vtot_psi_rank=volume_bin_radial(psi_rank);
        Npart_theta_binned(psi_rank,:)=volume_bin(psi_rank,:)/Vtot_psi_rank;
        
    end
    
    volume_bin_max=max(max(volume_bin));
    


    
    for psi_rank=1:N_radial_bins
        Npart_binned(psi_rank,:)=volume_bin(psi_rank,:)*Nalpha_binned_norm(psi_rank);
    end
    
    Npart_binned_tot=sum(sum(Npart_binned));
    Npart_binned=Npart_simulated*Npart_binned/Npart_binned_tot;
    
    % should be integers
    Npart_binned=round(Npart_binned);
    for psi_rank=1:N_radial_bins
        Npart_theta_binned(psi_rank,:)=round(sum(Npart_binned(psi_rank,:),2)*Npart_theta_binned(psi_rank,:));
    end
    
    for psi_rank=1:N_radial_bins
        Fpart_theta_binned(psi_rank,1)=Npart_theta_binned(psi_rank,end)+Npart_theta_binned(psi_rank,1);
        for theta_rank=1:N_theta_bins-1
            Fpart_theta_binned(psi_rank,theta_rank+1)=Npart_theta_binned(psi_rank,theta_rank)+Npart_theta_binned(psi_rank,theta_rank+1);
        end
        Fpart_theta_binned(psi_rank,end)=Npart_theta_binned(psi_rank,end)+Npart_theta_binned(psi_rank,1);
    end
    Fpart_theta_binned=(0.5*Fpart_theta_binned);
    Fpart_theta_binned(end,:)=Fpart_theta_binned(end-1,:);
    
    
    disp('smallest number of particles');
    Npart_min=min(min(Npart_theta_binned(Npart_theta_binned~=0)))
    
    % this gives us the ratio of particles to density
    Ne_psi_bins=interp1(1:Nradial,Ne_profile/Ne_profile(1),psi_bin_pos);
    Ne_psi_bins=Ne_psi_bins/Ne_psi_bins(1);
    
    

    Npart_flux_surface=sum(Npart_theta_binned,2);

    density_flux_surface=Npart_flux_surface./sum(volume_bin,2);
    density_flux_surface(isnan(density_flux_surface))=0;
    density_flux_surface=density_flux_surface./density_flux_surface(psi_pos_inf);
    
    Npart_flux_surface(end+1)=0.95*Npart_flux_surface(end);
    density_flux_surface(end+1)=0.95*density_flux_surface(end);
    
    
    Nalphas_simulated=sum(sum(Npart_theta_binned));
    
    disp('number of particles generated');
    Nalphas_simulated;
    
    alphas_pos_x=zeros(Nalphas_simulated,1);
    alphas_pos_z=zeros(Nalphas_simulated,1);
    alphas_pos_phi=zeros(Nalphas_simulated,1);
    alphas_Ekin=zeros(Nalphas_simulated,1);
    alphas_mm=zeros(Nalphas_simulated,1);
    alphas_vpll=zeros(Nalphas_simulated,1);
    
    Npart_rank_flux_surface_begin=1;
    
    
    % display options
    figure(1);
    set(gca,'FontSize',22);
    
    hold on;
    axis xy square
    
    for psi_rank=psi_pos_inf:psi_pos_sup
        
        psi_bin_inf=psi_bin_pos(psi_rank)-radial_bin_half_size;
        psi_bin_sup=psi_bin_pos(psi_rank)+radial_bin_half_size;
        
        %     Npart_flux_surface(psi_rank)=sum(Npart_theta_binned(psi_rank,:));
        Npart_rank_flux_surface_end=Npart_rank_flux_surface_begin+sum(Npart_theta_binned(psi_rank,:))-1;
        disp('Number of particles distributed in volume of flux surface element:');
        disp(Npart_flux_surface(psi_rank));
        
        
        DNpart=Npart_rank_flux_surface_end-Npart_rank_flux_surface_begin+1;
        psi_vector=ones(DNpart,1);
        theta_vector=zeros(DNpart,1);
        
        for(n=1:DNpart)
            psi_vector(n)=my_rand_linear_dist(Npart_flux_surface(psi_rank),Npart_flux_surface(psi_rank+1),psi_bin_inf,psi_bin_sup);
        end
        Nrank=1;
        for theta_rank=1:N_theta_bins
            DNpart_theta=Npart_theta_binned(psi_rank,theta_rank);
            for(n=Nrank:Nrank+DNpart_theta-1)
                theta_vector(n)=my_rand_linear_dist(Fpart_theta_binned(psi_rank,theta_rank),Fpart_theta_binned(psi_rank,theta_rank+1),theta_bin_bounds(theta_rank),theta_bin_bounds(theta_rank+1));
            end
            Nrank=Nrank+DNpart_theta;
        end
        if (length(find(isnan(theta_vector)))~=0)||(length(find(isnan(psi_vector)))~=0)
            find(isnan(theta_vector))
            find(isnan(psi_vector))
        end
        %     theta_vector=rand(Npart_flux_surface(psi_rank),1)*2*pi;
        
        figure(3);
        hold on;
        theta_vector_plot=histc(theta_vector,theta_bin_bounds);
        plot(theta_vector_plot(1:end-1))
        figure(4);
        hold on;
        psi_vector_plot=histc(psi_vector,psi_bin_inf:psi_bin_sup);
        plot(psi_bin_inf:psi_bin_sup-1,psi_vector_plot(1:end-1))
        
        Nrank=1;
        Epart_vector=ones(DNpart,1);
        Nalpha_energy_rank=ones(N_energy_bins,1);
        
        for E_rank=1:N_energy_bins-1
            Nalpha_energy_rank(E_rank)=round(f_D_energy_percentage(psi_rank,E_rank)*sum(Npart_theta_binned(psi_rank,:),2));
            for(n=Nrank:Nrank+Nalpha_energy_rank(E_rank)-1)
                % Evalues and energy_bin_size are in eV
                % uniform distribution on the bin
                % this an approximation but we take small enough bins so
                % that this should not impact on the final result
                Epart_vector(n)=my_rand_linear_dist(f_alpha_E_bin_edge_inf(psi_rank,E_rank),f_alpha_E_bin_edge_sup(psi_rank,E_rank),energy_D_range(E_rank)-0.5*energy_bin_size,energy_D_range(E_rank)+0.5*energy_bin_size);
                %Epart_vector(n)=Evalues(E_rank)+(rand(1)*(energy_bin_size)-0.5*energy_bin_size);
            end
            Nrank=Nrank+Nalpha_energy_rank(E_rank);
        end
        
        % remaining particles go into the highest bin
        E_rank=N_energy_bins;
        Nalpha_energy_rank(E_rank)=DNpart-Nrank+1;
        for(n=Nrank:Nrank+Nalpha_energy_rank(E_rank)-1)
            %Epart_vector(n)=Evalues(E_rank)+(rand(1)*(energy_bin_size)-0.5*energy_bin_size);
            Epart_vector(n)=my_rand_linear_dist(f_alpha_E_bin_edge_inf(psi_rank,E_rank),f_alpha_E_bin_edge_sup(psi_rank,E_rank),energy_D_range(E_rank)-0.5*energy_bin_size,energy_D_range(E_rank)+0.5*energy_bin_size);
        end
        
        if (n~=DNpart)
            disp('MISSING PARTICLES in E!!!!');
            n
        end
        % distribution in v_parallel range
        
        %vpll_vector=ones(1,1);
        vpll_sup_values=sqrt(2*(energy_D_range+0.5*energy_bin_size)*eV/mHe);
        
        n_vpll_count=0;
        Nrank=1;
        %vpll_vector=ones(DNpart,1);
        vpll_vector=ones(1,1);
        vpll_max=sqrt(2*(eV/mHe)*Epart_vector);
        
        for E_rank=1:N_energy_bins
            Nrank_Ebegin=Nrank;
            if Nalpha_energy_rank(E_rank)>0
                if vpll_vector(1)==1
                    max_vpll_Erank=vpll_max(1:Nalpha_energy_rank(E_rank));
                    vpll_vector=rand(Nalpha_energy_rank(E_rank),1).*(2*max_vpll_Erank)-max_vpll_Erank;
                    Nrank=Nrank+Nalpha_energy_rank(E_rank);
                else
                    max_vpll_Erank=vpll_max(Nrank:Nrank+Nalpha_energy_rank(E_rank)-1);
                    vpll_vector=[vpll_vector ; rand(Nalpha_energy_rank(E_rank),1).*(2*max_vpll_Erank)-max_vpll_Erank];
                    Nrank=Nrank+Nalpha_energy_rank(E_rank);
                end
            end
            
        end
        
        
        
        % happily shuffle everything
        permutation_vector=randperm(DNpart);
        psi_vector=psi_vector(permutation_vector);
        theta_vector=theta_vector(permutation_vector);
        Epart_vector=Epart_vector(permutation_vector);
        vpll_vector=vpll_vector(permutation_vector);
        
        alphas_pos_x(Npart_rank_flux_surface_begin:Npart_rank_flux_surface_end)=interp2(theta_scale,1:Nradial,X_PR_map',theta_vector,psi_vector);
        alphas_pos_z(Npart_rank_flux_surface_begin:Npart_rank_flux_surface_end)=interp2(theta_scale,1:Nradial,Z_PR_map',theta_vector,psi_vector);
        alphas_Ekin(Npart_rank_flux_surface_begin:Npart_rank_flux_surface_end)=Epart_vector;
        alphas_vpll(Npart_rank_flux_surface_begin:Npart_rank_flux_surface_end)=vpll_vector;
        
        
        figure(1);
        plot(alphas_pos_x(Npart_rank_flux_surface_begin:Npart_rank_flux_surface_end),alphas_pos_z(Npart_rank_flux_surface_begin:Npart_rank_flux_surface_end),'b.')
        DPSI_POS=0.5*(psi_bin_pos(2)-psi_bin_pos(1));
        contour(scale_X,scale_Z,psi_XZsmall_map',[psi_scale(psi_bin_pos(psi_rank)+DPSI_POS) psi_scale(psi_bin_pos(psi_rank)+DPSI_POS)],'r','linewidth',2);
        
%         plot(scale_X,Z_psi_fit_up_small(psi_bin_inf,:),'r');
%         plot(scale_X,Z_psi_fit_up_small(psi_bin_sup,:),'r');
%         plot(scale_X,Z_psi_fit_down_small(psi_bin_inf,:),'r');
%         plot(scale_X,Z_psi_fit_down_small(psi_bin_sup,:),'r');
%         
        
        Npart_rank_flux_surface_begin=Npart_rank_flux_surface_begin+DNpart;
        
    end
    
    alphas_pos_psi=interp2(scale_X,scale_Z,psi_norm_XZsmall_map',alphas_pos_x,alphas_pos_z);
    alphas_momentum=interp1(1:Nradial,rotation_profile,alphas_pos_psi);
    alphas_current=interp2(X_scale,Z_scale,JPHI_XZ_map',alphas_pos_x,alphas_pos_z);
    
    
    % filling up the uniformly distributes toroidal position values
    for(n=1:Nalphas_simulated)
        alphas_pos_phi(n)=my_rand(1)*2*pi;
    end
    
    %load('XZsmall_fields_tokamak_pre_collapse.mat', 'Btot_XZ_map');
    alphas_mm=(alphas_Ekin-0.5*(mHe*alphas_vpll.^2)/eV)./interp2(scale_X,scale_Z,Btot_XZ_map',alphas_pos_x,alphas_pos_z,'*linear');
    
    psi_avg=interp2(scale_X,scale_Z,  psi_XZsmall_map',alphas_pos_x,alphas_pos_z);
    
     
    NB_BINS=40
    [density_pp,psi_edges]=histcounts(psi_avg,NB_BINS);
    volume_edge_bins=interp1(psi_scale,volume_flux,psi_edges);
    volume_bins=volume_edge_bins(2:end)-volume_edge_bins(1:end-1);
    
    psi_middle=interp1(1:NB_BINS+1,psi_edges,round(NB_BINS/2));
    density_ref=interp1(psi_scale,C_density_profile,psi_middle);

    psi_bins=psi_edges(1:end-1)+0.5*(psi_edges(2)-psi_edges(1));
    nb_part_ref=interp1(psi_bins,density_pp,psi_middle);
    volume_bin_ref=interp1(psi_bins,volume_bins,psi_middle);
    
    density_part_ratio=density_ref/(nb_part_ref/volume_bin_ref);
    density_recalc=density_part_ratio*density_pp./volume_bins;
    
    
    
    FILENAME=strcat('initial_alphas_MB_CZ4_distribution',num2str(round(Nalphas_simulated/1e3)),'k.mat')
    save (FILENAME,'mHe','ZHe','rotation_profile','density_part_ratio','alphas_pos_x','alphas_pos_z','alphas_pos_phi','alphas_Ekin','alphas_mm','alphas_vpll','alphas_momentum','alphas_current','Nalphas_simulated');
    
    
    Nalphas_simulated
    
    figure;hold on
    grid on
    plot(psi_scale,C_density_profile,'r')
    plot(psi_bins,density_recalc,'b')


    Nalphas_simulated


