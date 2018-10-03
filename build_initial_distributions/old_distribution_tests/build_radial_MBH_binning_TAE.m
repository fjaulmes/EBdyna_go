
close all;
filename=strcat(DATA_FOLDER,'flux_geometry.mat');
load(filename);
filename=strcat(DATA_FOLDER,'physics_constants.mat');
load(filename);
%     mr=mHe*mHe/(mHe+mHe);
filename=strcat(DATA_FOLDER,'pressure_profile.mat');
load(filename);
filename=strcat(DATA_FOLDER,'motions_map_dimensions.mat');
load(filename);
filename=strcat(DATA_FOLDER,'tokamak_PR_map.mat');
load(filename);
%     filename=strcat(DATA_FOLDER,'TAE_data.mat');
%     load(filename);

%lowest possible value so that distribution matches n=10 case

rho_tor=(tor_flux_profile/max(tor_flux_profile));
nH0=1.44131e17
fast_density_profile=nH0*0.521298*exp(-(0.198739/0.298228)*tanh((sqrt(rho_tor)-0.49123)/0.198739));
mHe=mH
ZHe=1

density_on_axis=max(fast_density_profile)
% rho_pol=sqrt(1-(psi_scale/min(psi_scale)));
% fast_density_profile=(0.75e17)*exp(-(0.2/0.3) * tanh((rho_pol-0.5)/0.2));

MB_TEMP=500  % keV

radial_bin_size=8;
radial_bin_half_size=0.5*(radial_bin_size);
% simulation_size_r=round(1.38*size_r);
% simulation_size_r=round(simulation_size_r-mod(simulation_size_r,radial_bin_size))

radial_bins=[radial_bin_half_size+1:radial_bin_size:Nradial-4]';
N_radial_bins=size(radial_bins,1)


%typical JET like optimistic parameters
%no idea what could be actual effective charge here
%depends on the result of this script, right???


energy_bin_size=10.0*1e3; %keV
energy_max=4000*1e3; %keV
energy_H(1)=0.5*energy_bin_size;
N_energy_bins=round(energy_max/energy_bin_size)

for (x=2:N_energy_bins)
    energy_H(x)= energy_H(x-1)+energy_bin_size;
end
speed_H=sqrt(2*(eV/mH)*energy_H);

energy_H_range(1)=0;
for (x=2:N_energy_bins+1)
    energy_H_range(x)= energy_H_range(x-1)+energy_bin_size;
end
speed_H_range=sqrt(2*(eV/mH)*energy_H_range);

energy_range_global=(1:energy_max+2*energy_bin_size);
Ti_value=MB_TEMP*1e3
f_H_E_dist_global=(2*sqrt(energy_range_global*eV/pi).*(1/(Ti_value*eV))^(3/2)).*exp(-energy_range_global/Ti_value);
for n=1:N_energy_bins+1
    f_H_E_dist(n)=sum(f_H_E_dist_global((n-1)*energy_bin_size+1:n*energy_bin_size));
end

psi_pos_inf=1
psi_pos_sup=length(radial_bins)-1
Nalpha_binned=zeros(N_radial_bins,1);

for psi_pos=psi_pos_inf:length(radial_bins)-1
    
    disp('-------------------------------------------------')
    
    
    Ni_value=mean(fast_density_profile(radial_bins(psi_pos)-radial_bin_half_size:radial_bins(psi_pos)+radial_bin_half_size))
    
    disp('total density of H');
    disp(Ni_value);
    Nalpha_binned(psi_pos)=Ni_value;
    Ti_value=MB_TEMP*1e3;
    
%     for n=1:N_energy_bins+1
%         f_H_speed_dist(psi_pos,n)=((mHe/(2*pi*Ti_value*eV))^(3/2))*(4*pi*speed_H_range(n)^2)*exp(-(mHe*speed_H_range(n)^2)/(2*Ti_value*eV));
%         f_H_E_dist(psi_pos,n)=(2*sqrt(energy_H_range(n)*eV/pi).*(1/mHe*Ti_value*eV)^(3/2)).*exp(-energy_H_range(n)/Ti_value);
%     end
    for n=1:N_energy_bins
        f_H_E_dist_bin(psi_pos,n)=(2*sqrt(energy_H(n)*eV/pi).*(1/(Ti_value*eV))^(3/2)).*exp(-energy_H(n)/Ti_value);
    end
    
end

f_H_energy_bin_percentage=f_H_E_dist_bin*0;
f_H_energy_percentage=f_H_E_dist*0;
f_alpha_E_bin_edge_inf=f_H_E_dist*0;
f_alpha_E_bin_edge_sup=f_H_E_dist*0;

% for psi_pos=psi_pos_inf:psi_pos_sup
    
    
    % all bins should add up to 1
    max_density_norm=max(f_H_E_dist(1,:));
    f_H_energy_percentage(1,:)=f_H_E_dist(1,:)/max_density_norm;
    f_H_energy_percentage(1,:)=f_H_energy_percentage(1,:)/sum(f_H_energy_percentage(1,:),2);
    %     f_DT_energy_percentage(psi_pos,:)=f_DT_energy_percentage(psi_pos,:)/max_density_norm;
    
    f_alpha_E_bin_edge_inf(1,1:end-1)=interp1(energy_range_global,f_H_E_dist_global,energy_H-0.5*energy_bin_size);
    f_alpha_E_bin_edge_sup(1,1:end-1)=interp1(energy_range_global,f_H_E_dist_global,energy_H+0.5*energy_bin_size);
    f_alpha_E_bin_edge_inf(1,1)=0.2*f_alpha_E_bin_edge_inf(1,2);
    f_alpha_E_bin_edge_sup(1,end)=0.1*f_alpha_E_bin_edge_inf(1,end);
    
    f_H_energy_bin_percentage(1,:)=f_H_E_dist_bin(1,:)/max_density_norm;
    f_H_energy_bin_percentage(1,:)=f_H_energy_bin_percentage(1,:)/sum(f_H_energy_bin_percentage(1,:),2);
    
    
    figure(1);
    plot(energy_H,f_H_energy_bin_percentage(1,:));
    grid on;
    hold on;
    xlabel('E');
    ylabel('dN_\alpha / Ne_0 (%)');
    
% end


figure(4);
plot(radial_bins,Nalpha_binned);
xlabel('r');
ylabel('N_\alpha');


% completing the initial distribution with the distribution of the parallel
% speeds for a given energy level
N_vparallel_bins=31;
N_vparallel_bins_half=round(0.5*(N_vparallel_bins-1));
w0=max(energy_H_range)*eV
v0=sqrt(2*w0/mHe);

vpll_range=(v0)*(-N_vparallel_bins_half:N_vparallel_bins_half)*(1/N_vparallel_bins_half);
vpll_bin_size=v0/N_vparallel_bins_half

vpll_range_inf=(v0)*(-N_vparallel_bins_half:N_vparallel_bins_half)*(1/N_vparallel_bins_half)-0.5*vpll_bin_size;
vpll_range_sup=(v0)*(-N_vparallel_bins_half:N_vparallel_bins_half)*(1/N_vparallel_bins_half)+0.5*vpll_bin_size;


% 
% figure(2);
% axis xy;
% imagesc(radial_bins,energy_H,f_H_E_dist_bin');
% xlabel('r');
% ylabel('E');
% title('distribution of particles (normalized)')
% colorbar;




