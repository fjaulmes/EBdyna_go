
close all;

    filename=strcat(DATA_FOLDER,'physics_constants.mat');
    load(filename);
    mr=mDT*mHe/(mDT+mHe);
    filename=strcat(DATA_FOLDER,'pressure_profile.mat');
    load(filename);
    filename=strcat(DATA_FOLDER,'motions_map_dimensions.mat');
    load(filename);

    %lowest possible value 
    pTAE_inf=1
    pTAE_sup=size_r+30

radial_bin_size=8;
radial_bin_half_size=0.5*(radial_bin_size);
% simulation_size_r=round(1.38*size_r);
% simulation_size_r=round(simulation_size_r-mod(simulation_size_r,radial_bin_size))

radial_bins=[radial_bin_half_size+1:radial_bin_size:simulation_size_r-radial_bin_half_size+1]';
N_radial_bins=size(radial_bins,1)


%typical JET like optimistic parameters
frac_fusion=0.001;
Zeff=frac_fusion*2+(1-frac_fusion)*1;
%no idea what could be actual effective charge here
%depends on the result of this script, right???


Pvalues=zeros(N_radial_bins,1);
%Pvalue=interp1(1:size_r,P_initial_profile(1:size_r),radial_bins)
energy_bin_size=61;
half_energy_bin_size=0.5*(energy_bin_size-1);
energy_upper_cutoff=3.6*1e3;
%N_energy_bins_prev=4*round(energy_upper_cutoff/energy_bin_size)

psi_pos_inf=1
% psi_pos_inf=floor(interp1(radial_bins,1:N_radial_bins,pTAE_inf))
psi_pos_sup=ceil(interp1(radial_bins,1:N_radial_bins,pTAE_sup))
Nalpha_binned=zeros(N_radial_bins,1)

for psi_pos=psi_pos_inf:psi_pos_sup+1
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
    vthi=sqrt(2*eV*Te_psi/mDT);
    
    
    lambda_D=sqrt(epsilon0*Te_psi*eV/(Ne0_psi*eV^2));
    Lambda_ei=9*(4*pi/3*Ne0_psi*lambda_D^3);
    
    log_lambda=log(Lambda_ei);
    
    tau_e=(3*(2*pi)^1.5)*(epsilon0^2*sqrt(me)*(eV*Te_psi)^1.5)/(Ne0_psi*log_lambda*eV^4);
    tau_eHe=(3*(2*pi)^1.5)*(epsilon0^2*sqrt(me)*(eV*Te_psi)^1.5)/(Ne0_psi*log_lambda*4*eV^4);
    
    tau_s=(Ni0_psi*mDT/(Ne0_psi*me))*tau_eHe
    
    % Number of radial points
    NumberX=4000;
    
    % Main important values for description of the distribution
    w0=3.5*(10^6)*eV;
    v0=sqrt(2*w0/mHe);
    
    % approximative parameters
    vc=sqrt(44*2*Te_psi*eV/mHe);

    
    DeltaX=1/(NumberX-1);
	%Ti0 in keV
	averaged_cross_section=1.1*(1e-24)*(Ti_psi/1000)^2;
    Salpha=0.25*Ni0_psi*Ni0_psi*averaged_cross_section;
    Salpha=Salpha*Ti_psi*1e-4
    
    for (x=1:2*NumberX-1)
        speed_range(x)=(x-1)*DeltaX;
        vb(x)=vc*speed_range(x);
    end
    
    for (x=0.25*NumberX:2*NumberX-1)
        nubi(x)=(1/(4*pi))*(4*Ni0_psi*log_lambda*eV^4)/((epsilon0^2)*mr*mHe)*(1/(vb(x)^3+1.3*vthi^3));
        nube(x)=(1/(4*pi))*(4*Ne0_psi*log_lambda*eV^4)/((epsilon0^2)*me*mHe)*(1/(vb(x)^3+1.3*vthe^3));
    end
    
    % need to saturate collision rate
    % because the low speed part of the curve is not
    % what we have modelled and anyway we only use these curves
    % to find the critical speed transition from electrons to ions collision
    % dominated slowing down

    for (x=1:1000)
        nubi(x)=nubi(1000);
        nube(x)=nube(1000);
    end
    
    [vc_eps vc_pos]=min(abs(nubi-nube));
    vc=speed_range(vc_pos)*vc;
    wcrit=(0.5*mHe*vc^2)/eV;
    
    
    % redo the collision rate curves with corrected speed range
    for (x=1:2*NumberX-1)
        vb(x)=v0*speed_range(x);
        vbc(x)=vc*speed_range(x);
    end
    for (x=1000:2*NumberX-1)
        nubi(x)=(1/(4*pi))*(4*Ni0_psi*log_lambda*eV^4)/((epsilon0^2)*mr*mHe)*(1/(vbc(x)^3+1.3*vthi^3));
        nube(x)=(1/(4*pi))*(4*Ne0_psi*log_lambda*eV^4)/((epsilon0^2)*me*mHe)*(1/(vbc(x)^3+1.3*vthe^3));
    end
    for (x=1:1000)
        nubi(x)=nubi(1000);
        nube(x)=nube(1000);
    end
    
    
    ratio_v0_vc=v0/vc;
    %disp('ratio v0/vc =');
    %disp(ratio_v0_vc);
    
    % beyond critical speed
    v0_NumberX=ceil(ratio_v0_vc*NumberX);
    % beyond birth speed
    total_NumberX=ceil(1.12*v0_NumberX);
    vbc=zeros(total_NumberX,1);
    for (x=1:total_NumberX)
        speed_range(x)=(x-1)*DeltaX;
        vbc(x)=vc*speed_range(x);
    end
    
    s_var=Salpha*(mHe/mDT)*(tau_s/(2*pi*vc^3));
    % distribution below critical speed
    for (x=1:v0_NumberX)
        phi_s(x)=(2*mHe/mDT)*(1+speed_range(x)^3);
        w(x)=0.5*mHe*vc^2*speed_range(x)^2;
        f_alpha(x)=(s_var)/(phi_s(x));
    end
    
    
   % distribution for speed greater than critical speed
   for (x=v0_NumberX:total_NumberX)
        phi_s(x)=(2*mHe/mDT)*(1+speed_range(x)^3);
        w(x)=0.5*mHe*vc^2*speed_range(x)^2;
        f_alpha(x)=exp(-1/(Te_psi*eV)*(w(x)-w0))*(s_var)/(phi_s(x));
    end
    
    integ_alpha=0;
    for (x=1:total_NumberX)
        integ_alpha=integ_alpha+f_alpha(x)*DeltaX;
    end
    disp('total density of alphas');
    disp(integ_alpha);
    Nalpha_binned(psi_pos)=integ_alpha;
    
    NumberX_norm=round(2.5*NumberX)-1;
    % rescaling the distribution to have equal energy bins
    f_alpha_norm=interp1((1:total_NumberX)*NumberX_norm/total_NumberX,f_alpha,1:NumberX_norm);
    v_norm=interp1((1:total_NumberX)*NumberX_norm/total_NumberX,vbc,1:NumberX_norm);
    v_norm(1)=0;
    E_norm=(0.5*mHe/eV)*v_norm.^2;
    
    % binning in total kinetic energy values
    energy_bins_global=[1:energy_bin_size:round(2.04*NumberX)-1]';
    N_energy_bins=size(energy_bins_global,1);
    Evalues_global=E_norm(energy_bins_global);
    
    % 50 keV lower cutoff
    Erank_lower_cutoff=round(interp1(E_norm,1:size(E_norm,2), 50*1e3));
    energy_bins_global(1)=Erank_lower_cutoff;
    for (n_energy=2:N_energy_bins)
        energy_bins_global(n_energy)=energy_bins_global(n_energy-1)+energy_bin_size;
    end
    %N_energy_bins=size(energy_bins_global,1);
    Evalues_global=E_norm(energy_bins_global);
	
	% 3.8 MeV cutoff
    %Erank_cutoff=min(round(interp1(Evalues_global,1:N_energy_bins, 3.8*1e6))-2,N_energy_bins_prev);
    
    energy_bins=energy_bins_global(1:N_energy_bins);
    %N_energy_bins=size(energy_bins,1);
    %energy_bins=[1:N_energy_bins]*energy_bin_size;
    %E_norm=interp1((1:total_NumberX)*NumberX_norm/total_NumberX,vbc,1:NumberX_norm);
    Evalues_norm(psi_pos,:)=E_norm(energy_bins);
    
    for n=1:N_energy_bins
        f_alpha_energy_percentage_norm(psi_pos,n)=mean(f_alpha_norm(energy_bins(n)-half_energy_bin_size:energy_bins(n)+half_energy_bin_size));
    end
    alpha_cold_density(psi_pos)=max(f_alpha);
    max_alpha_density_norm=mean(f_alpha_energy_percentage_norm(psi_pos,:))/max(f_alpha);
    f_alpha_energy_percentage_norm(psi_pos,:)=f_alpha_energy_percentage_norm(psi_pos,:)*max_alpha_density_norm/max(f_alpha);

    f_alpha_energy_percentage_norm(psi_pos,:)=f_alpha_energy_percentage_norm(psi_pos,:)/sum(f_alpha_energy_percentage_norm(psi_pos,:),2);
%     frac_fusion=integ_alpha/Ne0_psi
    %N_energy_bins_prev=N_energy_bins;
end
N_energy_bins=N_energy_bins-1;
energy_bin_size=mean(mean(E_norm(:,energy_bins(2:end))-E_norm(:,energy_bins(1:end-1))));
Evalues=mean(E_norm(:,energy_bins(end-1)+half_energy_bin_size),1)-(0:N_energy_bins-1)*energy_bin_size;
Evalues(1:end)=Evalues(end:-1:1);

for psi_pos=psi_pos_inf:psi_pos_sup

    f_alpha_energy_percentage(psi_pos,:)=interp1(Evalues_norm(psi_pos,:),f_alpha_energy_percentage_norm(psi_pos,:),Evalues);
    max_alpha_density_norm=mean(f_alpha_energy_percentage(psi_pos,:))/max(f_alpha);
    f_alpha_energy_percentage(psi_pos,:)=f_alpha_energy_percentage(psi_pos,:)*max_alpha_density_norm/max(f_alpha);
    f_alpha_energy_percentage(psi_pos,:)=f_alpha_energy_percentage(psi_pos,:)/sum(f_alpha_energy_percentage(psi_pos,:),2);

    f_alpha_E_bin_edge_inf(psi_pos,:)=interp1(Evalues_norm(psi_pos,:),f_alpha_energy_percentage_norm(psi_pos,:),Evalues-0.5*energy_bin_size);
    f_alpha_E_bin_edge_sup(psi_pos,:)=interp1(Evalues_norm(psi_pos,:),f_alpha_energy_percentage_norm(psi_pos,:),Evalues+0.5*energy_bin_size);
    f_alpha_E_bin_edge_inf(psi_pos,1)=f_alpha_E_bin_edge_inf(psi_pos,2);
    f_alpha_E_bin_edge_sup(psi_pos,end)=0.1*f_alpha_E_bin_edge_inf(psi_pos,end);
    
    %max_alpha_density_norm=mean(f_alpha_E_bin_edge_inf(psi_pos,:))/max(f_alpha);
    f_alpha_E_bin_edge_inf(psi_pos,:)=f_alpha_E_bin_edge_inf(psi_pos,:)*max_alpha_density_norm/max(f_alpha);
    %max_alpha_density_norm=mean(f_alpha_E_bin_edge_sup(psi_pos,:))/max(f_alpha);
    f_alpha_E_bin_edge_sup(psi_pos,:)=f_alpha_E_bin_edge_sup(psi_pos,:)*max_alpha_density_norm/max(f_alpha);
    f_alpha_E_bin_edge_inf(psi_pos,:)=f_alpha_E_bin_edge_inf(psi_pos,:)/sum(f_alpha_E_bin_edge_inf(psi_pos,:),2);
    f_alpha_E_bin_edge_sup(psi_pos,:)=f_alpha_E_bin_edge_sup(psi_pos,:)/sum(f_alpha_E_bin_edge_sup(psi_pos,:),2);
    
    figure(1);
    plot(Evalues,f_alpha_energy_percentage(psi_pos,:));
    grid on;
    hold on;
    xlabel('E');
    ylabel('dN_\alpha / Ne_0 (%)');
    
end
energy_bin_size=mean(Evalues(2:end)-Evalues(1:end-1));
half_energy_bin_size=0.5*(energy_bin_size);

figure(4);
plot(radial_bins,Nalpha_binned);
xlabel('r');
ylabel('N_\alpha');


% completing the initial distribution with the distribution of the parallel
% speeds for a given energu level
N_vparallel_bins=31;
N_vparallel_bins_half=round(0.5*(N_vparallel_bins-1));
w0=3.5*1e6*eV;
v0=sqrt(2*w0/mHe);

vpll_range=(v0)*(-N_vparallel_bins_half:N_vparallel_bins_half)*(1/N_vparallel_bins_half);
vpll_bin_size=v0/N_vparallel_bins_half

vpll_range_inf=(v0)*(-N_vparallel_bins_half:N_vparallel_bins_half)*(1/N_vparallel_bins_half)-0.5*vpll_bin_size;
vpll_range_sup=(v0)*(-N_vparallel_bins_half:N_vparallel_bins_half)*(1/N_vparallel_bins_half)+0.5*vpll_bin_size;


for n=1:N_energy_bins
    w0=Evalues(n)*eV;
    w0_inf=(Evalues(n)-0.5*energy_bin_size)*eV;
    w0_sup=(Evalues(n)+0.5*energy_bin_size)*eV;
    theta0=sqrt(2*w0/mHe);
    theta0_inf=sqrt(2*w0_inf/mHe);
    theta0_sup=sqrt(2*w0_sup/mHe);
    %f_alpha_vpll=0.5*((4/3)*max(theta0^2-vpll_range.^2,0))/(theta0^3);
    build_falpha_vpll;
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
imagesc(radial_bins,Evalues,f_alpha_energy_percentage_norm');
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


