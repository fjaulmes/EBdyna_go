initialize_folder_names;

close all;

filename=strcat(DATA_FOLDER,'physics_constants.mat');
load(filename);
filename=strcat(DATA_FOLDER,'pressure_profile.mat');
load(filename);
filename=strcat(DATA_FOLDER,'motions_map_dimensions.mat');
if exist(filename)
    load(filename);
else
    run(strcat(EQ_FOLDER,'build_pre_collapse_XZ_maps.m'))
    filename=strcat(DATA_FOLDER,'motions_map_dimensions.mat');
    load(filename);
end



run('build_D_flat_density_profile.m')
simulation_size_r=512
simulation_start_r=420


radial_bin_size=2;
radial_bin_half_size=0.5*(radial_bin_size);

radial_bins=[simulation_start_r+1:radial_bin_size:simulation_size_r]';
radial_bins_range=[simulation_start_r+1-radial_bin_half_size:radial_bin_size:simulation_size_r+radial_bin_half_size]';
N_radial_bins=size(radial_bins,1)





psi_pos_inf=floor(interp1(radial_bins,1:N_radial_bins,radial_bins(1)))
psi_pos_sup=ceil(interp1(radial_bins,1:N_radial_bins,radial_bins(end)))

Nalpha_binned=zeros(N_radial_bins,1);

ENERGY_MIN=2.0 % minimum energy of 5 eV

energy_bin_size=60; %eV
energy_max=TempD+ENERGY_MIN;          %eV
energy_D(1)=0.5*energy_bin_size+ENERGY_MIN;
N_energy_bins=round((energy_max-ENERGY_MIN)/energy_bin_size)

if N_energy_bins>1
    for (x=2:N_energy_bins)
        energy_D(x)= energy_D(x-1)+energy_bin_size;
    end
end
speed_D=sqrt(2*(eV/mHe)*energy_D);

energy_D_range(1)=ENERGY_MIN;  
for (x=2:N_energy_bins+1)
    energy_D_range(x)= energy_D_range(x-1)+energy_bin_size;
end
speed_D_range=sqrt(2*(eV/mHe)*energy_D_range);


for psi_pos=psi_pos_inf:psi_pos_sup
    
    Nvalues(psi_pos)=mean(D_density_profile(radial_bins(psi_pos)-radial_bin_half_size:radial_bins(psi_pos)+radial_bin_half_size));
    
    
    Nalpha_binned(psi_pos)=Nvalues(psi_pos);    
    f_D_E_dist(psi_pos,1)=1;
    f_D_E_dist(psi_pos,2)=1;
    f_D_E_dist(psi_pos,3)=1;
    f_D_E_dist(psi_pos,4)=1;
    
end

f_D_energy_int=f_D_E_dist(:,1:end-1)*0;
f_D_energy_percentage=f_D_E_dist(:,1:end-1)*0;

for psi_pos=psi_pos_inf:psi_pos_sup
    for en=1:size(f_D_E_dist,2)-1
        Delta_en=energy_D_range(en+1)-energy_D_range(en);
        f_D_energy_int(psi_pos,en)=Delta_en*(f_D_E_dist(psi_pos,en+1)+0.5*(f_D_E_dist(psi_pos,en)-f_D_E_dist(psi_pos,en+1)));
    end
    
    % all bins should add up to 1
    max_density_norm=max(f_D_E_dist(psi_pos,:));
    f_D_energy_percentage(psi_pos,:)=f_D_energy_int(psi_pos,:)/sum(f_D_energy_int(psi_pos,:),2);

    f_alpha_E_bin_edge_inf(psi_pos,1)=1;
    f_alpha_E_bin_edge_sup(psi_pos,1)=1;
    f_alpha_E_bin_edge_inf(psi_pos,2)=1;
    f_alpha_E_bin_edge_sup(psi_pos,2)=1;
    f_alpha_E_bin_edge_inf(psi_pos,3)=1;
    f_alpha_E_bin_edge_sup(psi_pos,3)=1;

%%    
%     figure(1);
%     plot(energy_D,f_D_energy_percentage(psi_pos,:));
%     grid on;
%     hold on;
%     xlabel('E');
%     ylabel('dN_{D} / Ni_0 (%)');
    
end


% completing the initial distribution with the distribution of the parallel
% speeds for a given energy level
% N_vparallel_bins=61;
% N_vparallel_bins_half=round(0.5*(N_vparallel_bins-1));
% w0=max(energy_D_range*eV);
% v0=sqrt(2*w0/mHe);
% 
% vpll_range=(v0)*(-N_vparallel_bins_half:N_vparallel_bins_half)*(1/N_vparallel_bins_half);
% vpll_bin_size=v0/N_vparallel_bins_half
% 
% 
% vpll_range_inf=(v0)*(-N_vparallel_bins_half:N_vparallel_bins_half)*(1/N_vparallel_bins_half)-0.5*vpll_bin_size;
% vpll_range_sup=(v0)*(-N_vparallel_bins_half:N_vparallel_bins_half)*(1/N_vparallel_bins_half)+0.5*vpll_bin_size;
% 


figure(2);
axis xy;
imagesc(radial_bins,energy_D_range,f_D_energy_percentage');
xlabel('r');
ylabel('E');
title('distribution of particles (normalized)')
colorbar;




