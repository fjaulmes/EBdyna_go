
close all;

filename=strcat(DATA_FOLDER,'physics_constants.mat');
load(filename);
filename=strcat(DATA_FOLDER,'pressure_profile.mat');
load(filename);
filename=strcat(DATA_FOLDER,'motions_map_dimensions.mat');
load(filename);

   
    
    
% TE_profile_radial_ini=interp1((0:Nradial-1)/(Nradial-1),TE_profile_interp_ini,rho_tor_scale);
% TI_profile_radial_ini=interp1((0:Nradial-1)/(Nradial-1),TI_profile_interp_ini,rho_tor_scale);
% NE_profile_radial_ini=interp1((0:Nradial-1)/(Nradial-1),NE_profile_interp_ini,rho_tor_scale);
% NI_profile_radial_ini=interp1((0:Nradial-1)/(Nradial-1),NI_profile_interp_ini,rho_tor_scale);
% PTOT_profile_radial_ini=interp1((0:Nradial-1)/(Nradial-1),PTOT_profile_interp_ini,rho_tor_scale);
% OMEGA_profile_radial_ini=interp1((0:Nradial-1)/(Nradial-1),OMEGA_profile_interp_ini,rho_tor_scale);


Ti_profile=Te_profile;
Ni_profile=Ne_profile;
Ne0=Ne_profile(1)
Te0=Te_profile(1)

PI_radial_profile=Ni_profile.*Ti_profile;


radial_bin_size=8;
radial_bin_half_size=0.5*(radial_bin_size);
% simulation_size_r=round(1.38*size_r);
% simulation_size_r=round(simulation_size_r-mod(simulation_size_r,radial_bin_size))

radial_bins=[radial_bin_half_size+1:radial_bin_size:simulation_size_r+3*radial_bin_half_size+1]';
N_radial_bins=size(radial_bins,1)

    %lowest possible value so that distribution matches n=10 case
    pTAE_inf=radial_bins(1)
    pTAE_sup=radial_bins(end)

%typical JET like optimistic parameters
frac_fusion=0.001;
Zeff=frac_fusion*2+(1-frac_fusion)*1;
%no idea what could be actual effective charge here
%depends on the result of this script, right???


Pvalues=zeros(N_radial_bins,1);

psi_pos_inf=floor(interp1(radial_bins,1:N_radial_bins,pTAE_inf))
psi_pos_sup=ceil(interp1(radial_bins,1:N_radial_bins,pTAE_sup))-1

Nalpha_binned=zeros(N_radial_bins,1)

energy_bin_size=1*1e3; %keV
energy_max=40*1e3; %keV
energy_D(1)=0.5*energy_bin_size;
N_energy_bins=energy_max/energy_bin_size

for (x=2:N_energy_bins)
    energy_D(x)= energy_D(x-1)+energy_bin_size;
end
speed_D=sqrt(2*(eV/mD)*energy_D);

energy_D_range(1)=0;
for (x=2:N_energy_bins+1)
    energy_D_range(x)= energy_D_range(x-1)+energy_bin_size;
end
speed_D_range=sqrt(2*(eV/mD)*energy_D_range);


for psi_pos=psi_pos_inf:psi_pos_sup
    clear vb speed_range vbc f_alpha w
    
    disp('-------------------------------------------------')
    Pvalues(psi_pos)=mean(P_initial_profile(radial_bins(psi_pos)-radial_bin_half_size:radial_bins(psi_pos)+radial_bin_half_size));
    Nvalues(psi_pos)=mean(Ne_profile(radial_bins(psi_pos)-radial_bin_half_size:radial_bins(psi_pos)+radial_bin_half_size));
    
%     Ne0=4*10^19;
    Ne0_psi=Ne_profile(radial_bins(psi_pos));
    Ni0_psi=(1/Zeff)*Ne0_psi;
    Te_psi=0.5*Pvalues(psi_pos)/Ne0_psi/eV  %eV
    Ti_psi=0.5*Pvalues(psi_pos)/((1/Zeff)*Ne0_psi)/eV; %eV
	vthe=sqrt(2*eV*Te_psi/me);
    vthi=sqrt(2*eV*Te_psi/mD);
    
    
    disp('total density of D');
    disp(Ni0_psi);
    Nalpha_binned(psi_pos)=Ni0_psi;    
    
    for n=1:N_energy_bins+1
        f_D_speed_dist(psi_pos,n)=((mD/(2*pi*Ti_psi*eV))^(3/2))*(4*pi*speed_D_range(n)^2)*exp(-(mD*speed_D_range(n)^2)/(2*Ti_psi*eV));
        f_D_E_dist(psi_pos,n)=eV*energy_bin_size*(2*sqrt(energy_D_range(n)*eV/pi).*(1/(Ti_psi*eV))^(3/2)).*exp(-energy_D_range(n)/Ti_psi);
    end
    
end



for psi_pos=psi_pos_inf:psi_pos_sup

    % all bins should add up to 1
    max_density_norm=max(f_D_E_dist(psi_pos,:));
    f_D_energy_percentage(psi_pos,:)=f_D_E_dist(psi_pos,:)/max_density_norm;    
    f_D_energy_percentage(psi_pos,:)=f_D_energy_percentage(psi_pos,:)/sum(f_D_energy_percentage(psi_pos,:),2);
%     f_D_energy_percentage(psi_pos,:)=f_D_energy_percentage(psi_pos,:)/max_density_norm;

    f_alpha_E_bin_edge_inf(psi_pos,:)=interp1(energy_D_range,f_D_energy_percentage(psi_pos,:),energy_D-0.5*energy_bin_size);
    f_alpha_E_bin_edge_sup(psi_pos,:)=interp1(energy_D_range,f_D_energy_percentage(psi_pos,:),energy_D+0.5*energy_bin_size);
    f_alpha_E_bin_edge_inf(psi_pos,1)=f_alpha_E_bin_edge_inf(psi_pos,2);
    f_alpha_E_bin_edge_sup(psi_pos,end)=0.1*f_alpha_E_bin_edge_inf(psi_pos,end);

    
    figure(1);
    plot(energy_D_range,f_D_energy_percentage(psi_pos,:));
    grid on;
    hold on;
    xlabel('E');
    ylabel('dN_{D} / Ni_0 (%)');
    
end


figure(4);
plot(radial_bins,Nalpha_binned);
xlabel('r');
ylabel('N_\alpha');


% completing the initial distribution with the distribution of the parallel
% speeds for a given energu level
N_vparallel_bins=61;
N_vparallel_bins_half=round(0.5*(N_vparallel_bins-1));
w0=max(energy_D_range*eV);
v0=sqrt(2*w0/mD);

vpll_range=(v0)*(-N_vparallel_bins_half:N_vparallel_bins_half)*(1/N_vparallel_bins_half);
vpll_bin_size=v0/N_vparallel_bins_half

% vpll_range_inf=vpll_range(1:end-1);
% vpll_range_sup=vpll_range(2:end);

vpll_range_inf=(v0)*(-N_vparallel_bins_half:N_vparallel_bins_half)*(1/N_vparallel_bins_half)-0.5*vpll_bin_size;
vpll_range_sup=(v0)*(-N_vparallel_bins_half:N_vparallel_bins_half)*(1/N_vparallel_bins_half)+0.5*vpll_bin_size;


for n=1:N_energy_bins
    w0=energy_D(n)*eV;
    w0_inf=(energy_D(n)-0.5*energy_bin_size)*eV;
    w0_sup=(energy_D(n)+0.5*energy_bin_size)*eV;
    theta0=sqrt(2*w0/mD);
    theta0_inf=sqrt(2*w0_inf/mD);
    theta0_sup=sqrt(2*w0_sup/mD);
    %f_alpha_vpll=0.5*((4/3)*max(theta0^2-vpll_range.^2,0))/(theta0^3);
    build_falpha_vpll_D;
    f_alpha_vpll=f_alpha_vpll/max(f_alpha_vpll);
    f_alpha_vpll=f_alpha_vpll/mean(f_alpha_vpll)/size(f_alpha_vpll,2);
    f_alpha_vpll_percentage(n,:)=f_alpha_vpll;
    figure(3)
    plot(vpll_range,f_alpha_vpll);
    xlabel('v_{||} (m/s)');
    ylabel('fraction of total number of particles');
    grid on
    hold on
    
    f_alpha_vpll=0.5*((4/3)*max(theta0^2-vpll_range_inf.^2,0))/(theta0^3);
    f_alpha_vpll=f_alpha_vpll/max(f_alpha_vpll);
    f_alpha_vpll=f_alpha_vpll/mean(f_alpha_vpll)/size(f_alpha_vpll,2);
    f_alpha_vpll_edge_inf(n,:)=f_alpha_vpll;
    f_alpha_vpll=0.5*((4/3)*max(theta0^2-vpll_range_sup.^2,0))/(theta0^3);
    f_alpha_vpll=f_alpha_vpll/max(f_alpha_vpll);
    f_alpha_vpll=f_alpha_vpll/mean(f_alpha_vpll)/size(f_alpha_vpll,2);
    f_alpha_vpll_edge_sup(n,:)=f_alpha_vpll;
     
end
MIN_PROBA_VPLL=0.5*min(min(f_alpha_vpll_edge_inf(f_alpha_vpll_edge_inf~=0)));
% correcting edges



%f_alpha_percentage=f_alpha_percentage/max_alpha_density;

figure(2);
axis xy;
imagesc(radial_bins,energy_D_range,f_D_energy_percentage');
xlabel('r');
ylabel('E');
title('distribution of particles (normalized)')
colorbar;

% figure(2);
% grid on; hold on;
% plot(speed_range,max(nubi,nube),'k','LineWidth',2);
% plot(speed_range,nubi,'r--');
% plot(speed_range,nube,'b--');
% legend('\nu_{max}','\nu_i','\nu_e');
% xlabel('v_\alpha/v_c');
% ylabel('\nu');
% xlim([0.5 3]);


