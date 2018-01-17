clear all;

initialize_folder_names;
build_radial_alpha_particles_binning_lowEkin;

%pause
close all

clear pos_alpha_x pos_alpha_z

%Normalize the alpha density
Nalpha_max=max(Nalpha_binned);
Nalpha_binned_norm=Nalpha_binned/Nalpha_max;


% this number need to be distributed
% according to the density of particles 

filename=strcat(DATA_FOLDER,'volume_flux_geometry.mat');
load(filename);
filename=strcat(DATA_FOLDER,'XZsmall_fields_tokamak_pre_collapse.mat');
load(filename);

% rescaling flux surface curves on small map
X_min=mid_Xaxis_large-0.5*size(scale_X,2);
X_max=mid_Xaxis_large+0.5*size(scale_X,2)-1;

X1_Nradial=X1_Nradial-mid_Xaxis_large+mid_X;
X2_Nradial=X2_Nradial-mid_Xaxis_large+mid_X;
Z_psi_fit_up_small=Z_psi_fit_up(:,X_min:X_max);
Z_psi_fit_down_small=Z_psi_fit_down(:,X_min:X_max);


X_bin_size=24;
X_bin_half=0.5*(X_bin_size);
X_min=X1_Nradial(simulation_size_r);
X_max=X2_Nradial(simulation_size_r);

X_offset=-round(0.5*mod(X_max-X_min,X_bin_size))+1+X_bin_half;
if mod(X_max-X_min,X_bin_size)==0
    N_X_bins=ceil((X_max-X_min)/X_bin_size);
else
    N_X_bins=ceil((X_max-X_min)/X_bin_size)+1;
end

X_bin_pos=X_min+X_offset+X_bin_size*(0:N_X_bins-1)
X_bin_pos_remap=X_bin_pos+mid_Xaxis_large-mid_X;

% loop through all flux surface bins
delta_psi=radial_bin_half_size;
psi_bin_pos=1+radial_bin_half_size+radial_bin_size*(0:N_radial_bins-1)


% volume of vertical sclices of toroidal dougnnuts shaped bins
% (yummy)

for psi_rank=1:N_radial_bins
    for X_bin_rank=1:N_X_bins
        volume_bin(psi_rank,X_bin_rank)=sum(sum(volume_tor_diff(psi_bin_pos(psi_rank)-delta_psi:psi_bin_pos(psi_rank)+delta_psi,X_bin_pos_remap(X_bin_rank)-X_bin_half:X_bin_pos_remap(X_bin_rank)+X_bin_half),1),2);
    end
end

volume_bin_max=max(max(volume_bin));

for psi_rank=1:N_radial_bins
    Npart_binned(psi_rank,:)=volume_bin(psi_rank,:)*Nalpha_binned_norm(psi_rank);
end

% there should be no less than 2 particles for each energy in a volume
Npart_min=round(0.023*N_energy_bins)+1;
Npart_binned=round(Npart_min*Npart_binned/min(Npart_binned(Npart_binned~=0)));


% this gives us the ratio of particles to density
density_part_ratio=Nalpha_binned(N_radial_bins-1)/sum(Npart_binned(N_radial_bins-1,:),2)

Nalphas_simulated=sum(sum(Npart_binned))
alphas_pos_x=zeros(Nalphas_simulated,1);
alphas_pos_z=zeros(Nalphas_simulated,1);
alphas_pos_phi=zeros(Nalphas_simulated,1);
alphas_Ekin=zeros(Nalphas_simulated,1);
alphas_mm=zeros(Nalphas_simulated,1);
alphas_vpll=zeros(Nalphas_simulated,1);

Npart_rank_begin=1;
Npart_rank_flux_surface_begin=1;


% display options
figure(1);
hold on;
axis xy square

%Z_bin_vector=zeros(N_radial_bins*N_X_bins*2,1);
Z_index=1;

for psi_rank=1:N_radial_bins
    
    psi0=psi_bin_pos(psi_rank);
    X_min=X1_Nradial(psi0+delta_psi);
    X_max=X2_Nradial(psi0+delta_psi);
    
    Z_sup_max=Z_psi_fit_up_small(psi0+delta_psi,:);
    Z_sup_min=Z_psi_fit_up_small(psi0-delta_psi,:);
    
    Z_inf_min=Z_psi_fit_down_small(psi0-delta_psi,:);
    Z_inf_max=Z_psi_fit_down_small(psi0+delta_psi,:);
    
    Npart_flux_surface(psi_rank)=sum(Npart_binned(psi_rank,:));
    Npart_rank_flux_surface_end=Npart_rank_flux_surface_begin+sum(Npart_binned(psi_rank,:))-1;
    disp('Number of particles distributed in volume of flux surface element:');
    disp(Npart_flux_surface(psi_rank));
            
    for X_bin_rank=1:N_X_bins
        
        X0=X_bin_pos(X_bin_rank);
        
        if (Npart_binned(psi_rank,X_bin_rank)~=0)
            
            Z_sup_max_value=max(Z_sup_max(X0-X_bin_half-1:X0+X_bin_half+1));
            Z_sup_min_value=min(Z_sup_min(X0-X_bin_half-1:X0+X_bin_half+1));
            Z_sup_pos=0.5*(Z_sup_max_value+Z_sup_min_value);
            
            Z_inf_max_value=max(abs(Z_inf_max(X0-X_bin_half-1:X0+X_bin_half+1)));
            Z_inf_min_value=min(abs(Z_inf_min(X0-X_bin_half-1:X0+X_bin_half+1)));
            Z_inf_pos=-0.5*(Z_inf_max_value+Z_inf_min_value);
            
            x_min=interp1(1:size(scale_X,2),scale_X,X0-X_bin_half,'*linear');
            x_max=interp1(1:size(scale_X,2),scale_X,X0+X_bin_half,'*linear');
            
            z_min=1;
            z_max=-1;
            z=0;
            Npart_rank_end=Npart_rank_begin+Npart_binned(psi_rank,X_bin_rank)-1;
            

            
            for(n=Npart_rank_begin:Npart_rank_end)
                % random coordinate on [x_min ; x_max[
                while (z<z_min)||(z>=z_max)
                    x=my_rand(1)*(x_max-x_min)+x_min;
                    if round(rand(1))==1
                        % upper part
                        sign_z=1;
                        z_min=interp1(scale_X(X0-X_bin_half-1:X0+X_bin_half+1),Z_sup_min(X0-X_bin_half-1:X0+X_bin_half+1),x,'*linear');
                        z_max=interp1(scale_X(X0-X_bin_half-1:X0+X_bin_half+1),Z_sup_max(X0-X_bin_half-1:X0+X_bin_half+1),x,'*linear');
                        % upper part
                        z=sign_z*(my_rand(1)*(Z_sup_max_value-Z_sup_min_value)+Z_sup_min_value);
                    else
                        % lower part
                        sign_z=-1;
                        z_min=interp1(scale_X(X0-X_bin_half-1:X0+X_bin_half+1),Z_inf_max(X0-X_bin_half-1:X0+X_bin_half+1),x,'*linear');
                        z_max=interp1(scale_X(X0-X_bin_half-1:X0+X_bin_half+1),Z_inf_min(X0-X_bin_half-1:X0+X_bin_half+1),x,'*linear');
                         % lower part
                        z=sign_z*(my_rand(1)*(Z_inf_max_value-Z_inf_min_value)+Z_inf_min_value);
                    end
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
            
            Npart_rank_begin=Npart_rank_begin+Npart_binned(psi_rank,X_bin_rank);
            Z_volume_bin(psi_rank,X_bin_rank,1)=Z_sup_pos;
            Z_volume_bin(psi_rank,X_bin_rank,2)=Z_inf_pos;
            Z_bin_vector(Z_index)=Z_sup_pos;
            Z_bin_vector(Z_index+1)=Z_inf_pos;  
            Z_index=Z_index+2;
        end
        
    end
    Z_bin_vector=sort(Z_bin_vector);
    
    DNpart=Npart_rank_flux_surface_end-Npart_rank_flux_surface_begin+1;
     
    Nrank=1;
    Epart_vector=ones(DNpart,1);
    Nalpha_energy_rank=ones(N_energy_bins,1);
    
    for E_rank=1:N_energy_bins-1
        Nalpha_energy_rank(E_rank)=round(f_alpha_energy_percentage(psi_rank,E_rank)*sum(Npart_binned(psi_rank,:),2));
        for(n=Nrank:Nrank+Nalpha_energy_rank(E_rank)-1)
            % Evalues and energy_bin_size are in eV
            % uniform distribution on the bin
            % this an approximation but we take small enough bins so
            % that this should not impact on the final result
%             Epart_vector(n)=Evalues(E_rank)+(rand(1)*(energy_bin_size)-0.5*energy_bin_size);
            Epart_vector(n)=my_rand_linear_dist(f_alpha_E_bin_edge_inf(psi_rank,E_rank),f_alpha_E_bin_edge_sup(psi_rank,E_rank),Evalues(E_rank)-0.5*energy_bin_size,Evalues(E_rank)+0.5*energy_bin_size);
        end
        Nrank=Nrank+Nalpha_energy_rank(E_rank);
    end
    
    % remaining particles go into the highest bin
    E_rank=N_energy_bins;
    Nalpha_energy_rank(E_rank)=DNpart-Nrank+1;
    for(n=Nrank:Nrank+Nalpha_energy_rank(E_rank)-1)
        Epart_vector(n)=Evalues(E_rank)+(rand(1)*(energy_bin_size)-0.5*energy_bin_size);
    end
    
    
    % distribution in v_parallel range
    
    vpll_max=sqrt(2*(eV/mHe)*Epart_vector);
%     vpll_max=sqrt(2*(Epart_vector+0.5*energy_bin_size)*eV/mHe);
    %vpll_vector=ones(1,1);
    vpll_sup_values=sqrt(2*(Evalues+0.5*energy_bin_size)*eV/mHe);
    
    n_vpll_count=0;
    Nrank=1;
    vpll_vector=ones(1,1);
    
    E_rank=1;
    max_vpll_Erank=vpll_max(1:Nalpha_energy_rank(E_rank));
    vpll_vector=rand(Nalpha_energy_rank(E_rank),1).*(2*max_vpll_Erank)-max_vpll_Erank;
    Nrank=Nrank+Nalpha_energy_rank(E_rank);
    
    for E_rank=2:N_energy_bins
%         vpll_count=1;
%         Nrank_Ebegin=Nrank;

%         if size(vpll_vector,1)<2
%             max_vpll_Erank=vpll_max(1:Nalpha_energy_rank(E_rank));
%             vpll_vector=rand(Nalpha_energy_rank(E_rank),1).*(2*max_vpll_Erank)-max_vpll_Erank;
%             Nrank=Nrank+Nalpha_energy_rank(E_rank);
%         else
            max_vpll_Erank=vpll_max(Nrank:Nrank+Nalpha_energy_rank(E_rank)-1);
            vpll_vector=[vpll_vector ; rand(Nalpha_energy_rank(E_rank),1).*(2*max_vpll_Erank)-max_vpll_Erank];
            Nrank=Nrank+Nalpha_energy_rank(E_rank);
%         end        
%         w0=Evalues(E_rank)*eV;
%         theta0=sqrt(2*w0/mHe);
%         f_alpha_vpll=0.5*((4/3)*max(theta0^2-vpll_range.^2,0))/(theta0^3);
%         f_alpha_vpll=f_alpha_vpll/max(f_alpha_vpll);
%         f_alpha_vpll=f_alpha_vpll/mean(f_alpha_vpll)/size(f_alpha_vpll,2);

        % filling negative parallel speeds
    
%         for(vpll_rank=1:N_vparallel_bins_half)
%             Nalpha_vpll_rank=round(f_alpha_vpll_percentage(E_rank,vpll_rank)*Nalpha_energy_rank(E_rank));
%             %vpll_value=vpll_range(vpll_rank)-0.5*vpll_bin_size;
%             if (Nalpha_vpll_rank>=1) % && (abs(vpll_value)<vpll_sup_values(E_rank))
%                 for (n=Nrank:Nrank+Nalpha_vpll_rank-1)
%                     % uniform distribution on the bin
%                     % taking into account the energy value of the bin
%                     if sqrt(Epart_vector(n)/Evalues(E_rank))*abs(vpll_range(vpll_rank))-0.5*vpll_bin_size<vpll_max(n)
% %                         vpll_sat=sqrt(2*(eV/mHe)*alphas_Ekin(Nrank));
% %                         vpll_value=vpll_sat;
%                         vpll_value=vpll_max(n);
%                         while(abs(vpll_value)>=vpll_max(n))
%                             vpll_value=sqrt(Epart_vector(n)/Evalues(E_rank))*vpll_range(vpll_rank)-(my_rand(1)*(vpll_bin_size)-0.5*vpll_bin_size);
%                         end
%                         vpll_vector(n)=vpll_value;
%                     else
%                         %add to the middle bin
%                         vpll_vector(n)=vpll_range(N_vparallel_bins_half+1)+(my_rand(1)*(vpll_bin_size)-0.5*vpll_bin_size);
%                         disp('missing particles put in middle bin!')
%                     end
%                     vpll_count=vpll_count+1;
%                 end
%                 Nrank=Nrank+Nalpha_vpll_rank;
%             end
%             
%         end
%         % filling positive parallel speeds
%         for(vpll_rank=N_vparallel_bins_half+2:N_vparallel_bins)
%             Nalpha_vpll_rank=round(f_alpha_vpll_percentage(E_rank,vpll_rank)*Nalpha_energy_rank(E_rank));
%             %vpll_value=vpll_range(vpll_rank);
%             if (Nrank<=DNpart) && (Nalpha_vpll_rank>=1)  && (Nrank+Nalpha_vpll_rank-1<=length(Epart_vector))
%                 for (n=Nrank:Nrank+Nalpha_vpll_rank-1)
%                     % uniform distribution on the bin
%                     % this an approximation but we take small enough bins so
%                     % that this should not impact on the final result
%                     if sqrt(Epart_vector(n)/Evalues(E_rank))*vpll_range(vpll_rank)-0.5*vpll_bin_size<vpll_max(n)
% %                         vpll_sat=(2*(eV/mHe)*alphas_Ekin(Nrank));
% %                         vpll_value=vpll_sat;
%                         vpll_value=vpll_max(n);
%                         while(vpll_value>=vpll_max(n))
%                             vpll_value=sqrt(Epart_vector(n)/Evalues(E_rank))*vpll_range(vpll_rank)+(my_rand(1)*(vpll_bin_size)-0.5*vpll_bin_size);
%                         end
%                         vpll_vector(n)=vpll_value;
%                     else
%                         %add to the middle bin
%                         vpll_vector(n)=vpll_range(N_vparallel_bins_half+1)+(my_rand(1)*(vpll_bin_size)-0.5*vpll_bin_size);
%                         disp('missing particles put in middle bin!')
%                     end
%                     vpll_count=vpll_count+1;
%                 end
%                 Nrank=Nrank+Nalpha_vpll_rank;
%             end
%             
%         end
%         % remaining particles go in middle bin
%         vpll_rank=N_vparallel_bins_half+1;
%         Nalpha_vpll_rank=Nalpha_energy_rank(E_rank)-vpll_count+1;
%         for (n=Nrank:Nrank+Nalpha_vpll_rank-1)
%             % uniform distribution on the bin
%             vpll_value=vpll_max(n);
%             while(vpll_value>=vpll_max(n))
%                 vpll_value=vpll_range(vpll_rank)+(my_rand(1)*(vpll_bin_size)-0.5*vpll_bin_size);
%                 vpll_count=vpll_count+1;
%             end
%             vpll_vector(n)=vpll_value;
%         end
%         Nrank=Nrank+Nalpha_vpll_rank;
%         if size(vpll_vector,1)<2
%             vpll_vector=vpll_vector;
%         else
%             vpll_vector=[vpll_vector ; vpll_vector];
%         end
    end
    
%     disp('pause')
%     pause
    
    % remaining particles go in middle bin
    vpll_rank=N_vparallel_bins_half+1;
    for (n=Nrank:DNpart)
        % uniform distribution on the bin
        vpll_vector(n)=vpll_range(vpll_rank)+my_rand(1)*(vpll_bin_size)-0.5*vpll_bin_size;
        disp('missing particles put in middle bin!');
    end
    
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
end


% filling up the uniformly distributes toroidal position values
for(n=1:Nalphas_simulated)
    alphas_pos_phi(n)=my_rand(1)*2*pi;
end


alphas_mm=(alphas_Ekin-0.5*(mHe*alphas_vpll.^2)/eV)./interp2(scale_X,scale_Z,Btot_XZ_map',alphas_pos_x,alphas_pos_z,'*linear');

NEGATIVE_MM=find(alphas_mm<0);

disp('Beware of negative mm!!!!!!')
disp('length(NEGATIVE_MM)');
disp(length(NEGATIVE_MM));


% figure(1);
% axis xy square
% hold on
% 
% plot(scale_X(X0-X_bin_half:X0+X_bin_half),Z_sup_max(X0-X_bin_half:X0+X_bin_half));
% plot(scale_X(X0-X_bin_half:X0+X_bin_half),Z_sup_min(X0-X_bin_half:X0+X_bin_half));
% plot(scale_X(X0-X_bin_half:X0+X_bin_half),Z_inf_max(X0-X_bin_half:X0+X_bin_half));
% plot(scale_X(X0-X_bin_half:X0+X_bin_half),Z_inf_min(X0-X_bin_half:X0+X_bin_half));
% 
% figure(2);
% plot(pos_alpha_x,pos_alpha_z,'.');


save initial_alphas_lowEkin_distribution_200.mat alphas_pos_x alphas_pos_z alphas_pos_phi alphas_Ekin alphas_mm alphas_vpll Nalphas_simulated
% save speed_bins.mat Evalues N_energy_bins vpll_range energy_bin_size N_vparallel_bins vpll_bin_size
% save geomtry_bins.mat density_part_ratio Z_bin_vector Z_volume_bin volume_bin Npart_binned psi_bin_pos N_radial_bins N_X_bins X_bin_size radial_bin_size X_bin_pos Z_psi_fit_up_small Z_psi_fit_down_small
%save('map_dimensions.mat','mid_Xaxis_large','mid_Xzero','DPHI','DX','NP','NB_PHI','NB_PHI_DATA_HALF','R0','a','XX_small','ZZ_small','RR','Z_PR_map','scale_X','scale_Z','simulation_size_r','X_axis','Z_axis','mid_X','mid_Z','finesse_mesh','finesse_mesh_dtri','IDTRI','XZ_mesh','XZ_mesh_dtri','Raxis','simulation_size_r');

Nalphas_simulated
    close all
hist(alphas_vpll,20)