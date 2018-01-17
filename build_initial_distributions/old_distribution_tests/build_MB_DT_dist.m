close all;
clear f_DT_speed_dist f_DT_energy_percentage
clear f_alpha_E_bin_edge_inf f_alpha_E_bin_edge_sup
clear f_alpha_vpll_edge_inf f_alpha_vpll_edge_sup f_alpha_vpll_percentage vpll_range



filename=strcat(DATA_FOLDER,'physics_constants.mat');
load(filename);
mr=mDT*mHe/(mDT+mHe);
%mDT=mD
mHe=mDT
ZHe=1

filename=strcat(DATA_FOLDER,'pressure_profile.mat');
load(filename);
filename=strcat(DATA_FOLDER,'motions_map_dimensions.mat');
load(filename);

BETA_ALPHAS=mean(P_initial_profile./(Te_profile.*Ne_profile))-2
TE_FRAC=1/(2+BETA_ALPHAS)

radial_bin_half_size=0.5*(radial_bin_size);

radial_bins=[radial_bin_half_size+1:radial_bin_size:simulation_size_r-30]';
N_radial_bins=size(radial_bins,1)

%typical JET like optimistic parameters
Zeff=1.1

Pvalues=zeros(N_radial_bins,1);

energy_bin_size=1.3*1e3; %keV
energy_max=130*1e3; %keV
energy_DT(1)=0.5*energy_bin_size;
N_energy_bins=energy_max/energy_bin_size

for (x=2:N_energy_bins)
    energy_DT(x)= energy_DT(x-1)+energy_bin_size;
end
speed_DT=sqrt(2*(eV/mDT)*energy_DT);

energy_DT_range(1)=0;
for (x=2:N_energy_bins+1)
    energy_DT_range(x)= energy_DT_range(x-1)+energy_bin_size;
end
speed_DT_range=sqrt(2*(eV/mDT)*energy_DT_range);

for psi_pos=1:N_radial_bins
    disp('-------------------------------------------------')
    
    Pvalues(psi_pos)=mean(P_initial_profile(radial_bins(psi_pos)-radial_bin_half_size:radial_bins(psi_pos)+radial_bin_half_size));
    Nvalues(psi_pos)=mean(Ne_profile(radial_bins(psi_pos)-radial_bin_half_size:radial_bins(psi_pos)+radial_bin_half_size));
    
    Ni_value=(1/Zeff)*Nvalues(psi_pos);
    Te0=Te_profile(1)/eV;
    
    Te_value=TE_FRAC*(Pvalues(psi_pos)./Nvalues(psi_pos))/eV;  %eV
%     Te_value=Te0;
    Ti_value=Te_value; %eV     % Ti=Te
    Tvalues(psi_pos)=Ti_value;
	vthe=sqrt(2*eV*Te_value/me);
    vthi=sqrt(2*eV*Ti_value/mDT);
    
    disp('total density of DT');
    disp(Ni_value);
    Nalpha_binned(psi_pos)=Ni_value;    
    
    for n=1:N_energy_bins+1
        f_DT_speed_dist(psi_pos,n)=((mDT/(2*pi*Ti_value*eV))^(3/2))*(4*pi*speed_DT_range(n)^2)*exp(-(mDT*speed_DT_range(n)^2)/(2*Ti_value*eV));
        f_DT_E_dist(psi_pos,n)=(2*sqrt(energy_DT_range(n)*eV/pi).*(1/mDT*Ti_value*eV)^(3/2)).*exp(-energy_DT_range(n)/Ti_value);
    end

end



for psi_pos=1:N_radial_bins

    
    % all bins should add up to 1
    max_density_norm=max(f_DT_E_dist(psi_pos,:));
    f_DT_energy_percentage(psi_pos,:)=f_DT_E_dist(psi_pos,:)/max_density_norm;    
    f_DT_energy_percentage(psi_pos,:)=f_DT_energy_percentage(psi_pos,:)/sum(f_DT_energy_percentage(psi_pos,:),2);
%     f_DT_energy_percentage(psi_pos,:)=f_DT_energy_percentage(psi_pos,:)/max_density_norm;

    f_alpha_E_bin_edge_inf(psi_pos,:)=interp1(energy_DT_range,f_DT_energy_percentage(psi_pos,:),energy_DT-0.5*energy_bin_size);
    f_alpha_E_bin_edge_sup(psi_pos,:)=interp1(energy_DT_range,f_DT_energy_percentage(psi_pos,:),energy_DT+0.5*energy_bin_size);
    f_alpha_E_bin_edge_inf(psi_pos,1)=f_alpha_E_bin_edge_inf(psi_pos,2);
    f_alpha_E_bin_edge_sup(psi_pos,end)=0.1*f_alpha_E_bin_edge_inf(psi_pos,end);
    
    %max_alpha_density_norm=mean(f_alpha_E_bin_edge_inf(psi_pos,:))/max(f_alpha);
%     f_alpha_E_bin_edge_inf(psi_pos,:)=f_alpha_E_bin_edge_inf(psi_pos,:)*max_alpha_density_norm/max(f_alpha);
    %max_alpha_density_norm=mean(f_alpha_E_bin_edge_sup(psi_pos,:))/max(f_alpha);
%     f_alpha_E_bin_edge_sup(psi_pos,:)=f_alpha_E_bin_edge_sup(psi_pos,:)*max_alpha_density_norm/max(f_alpha);
%     f_alpha_E_bin_edge_inf(psi_pos,:)=f_alpha_E_bin_edge_inf(psi_pos,:)/sum(f_alpha_E_bin_edge_inf(psi_pos,:),2);
%     f_alpha_E_bin_edge_sup(psi_pos,:)=f_alpha_E_bin_edge_sup(psi_pos,:)/sum(f_alpha_E_bin_edge_sup(psi_pos,:),2);
    
    figure(2);
    plot(energy_DT_range,f_DT_energy_percentage(psi_pos,:));
    grid on;
    hold on;
    xlabel('E');
    ylabel('dN_{DT} / Ni_0 (%)');

%     figure(2);
%     plot(energy_DT,f_alpha_E_bin_edge_sup(psi_pos,:)-f_alpha_E_bin_edge_inf(psi_pos,:),'b');
% %     plot(energy_DT,f_alpha_E_bin_edge_sup(psi_pos,:),'r');
%     grid on;
%     hold on;
%     xlabel('E');
%     ylabel('dN_{DT} / Ni_0 (%)');

end

for psi_pos=1:N_radial_bins
    cum_f_DT_energy_percentage(1,psi_pos)=0;
    for n=2:N_energy_bins+1
%         f_DT_energy_percentage(psi_pos,:)=f_DT_speed_dist(psi_pos,:)/max_density_norm;
        cum_f_DT_energy_percentage(n,psi_pos)=cum_f_DT_energy_percentage(n-1,psi_pos)+f_DT_energy_percentage(psi_pos,n);
    end
    cum_f_DT_energy_percentage(:,psi_pos)=cum_f_DT_energy_percentage(:,psi_pos)/max(cum_f_DT_energy_percentage(:,psi_pos));
end


for psi_pos=1:N_radial_bins
    Tval_recalc(psi_pos)=(2/3)*interp1(cum_f_DT_energy_percentage(1:end-5,psi_pos),energy_DT_range(1:end-5),0.5);
end
energy_rescale=Tval_recalc(1)/Tvalues(1)
% energy_DT_range=energy_DT_range/energy_rescale;
% energy_bin_size=energy_bin_size/energy_rescale;
% energy_DT=energy_DT/energy_rescale;
% Evalues=energy_rescale*Evalues;

% completing the initial distribution with the distribution of the parallel
% speeds for a given energu level
N_vparallel_bins=61;
N_vparallel_bins_half=round(0.5*(N_vparallel_bins-1));
w0=max(energy_DT_range*eV);
v0=sqrt(2*w0/mDT);

vpll_range=(v0)*(-N_vparallel_bins_half:N_vparallel_bins_half)*(1/N_vparallel_bins_half);
vpll_bin_size=v0/N_vparallel_bins_half

% vpll_range_inf=vpll_range(1:end-1);
% vpll_range_sup=vpll_range(2:end);

vpll_range_inf=(v0)*(-N_vparallel_bins_half:N_vparallel_bins_half)*(1/N_vparallel_bins_half)-0.5*vpll_bin_size;
vpll_range_sup=(v0)*(-N_vparallel_bins_half:N_vparallel_bins_half)*(1/N_vparallel_bins_half)+0.5*vpll_bin_size;


for n=1:N_energy_bins
    w0=energy_DT(n)*eV;
    w0_inf=(energy_DT(n)-0.5*energy_bin_size)*eV;
    w0_sup=(energy_DT(n)+0.5*energy_bin_size)*eV;
    theta0=sqrt(2*w0/mDT);
    theta0_inf=sqrt(2*w0_inf/mDT);
    theta0_sup=sqrt(2*w0_sup/mDT);
    %f_alpha_vpll=0.5*((4/3)*max(theta0^2-vpll_range.^2,0))/(theta0^3);
    build_falpha_vpll_DT;
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




