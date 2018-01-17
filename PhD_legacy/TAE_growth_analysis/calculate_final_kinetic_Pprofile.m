close all

load('initialG_DT_MB_all_pre_collapse.mat')
load('initialG_DT_MB_all_precession_stats.mat')
alphas_psi_ini=alphas_psi;
alphas_r_ini=interp1(1:Nradial,radial_r_value_flux,alphas_psi_ini);
alphas_pphi_ini=alphas_pphi0;
pphi_min=min(alphas_pphi_ini);
pphi_max=max(alphas_pphi_ini);
psi_min=min(alphas_psi_ini);
psi_max=max(alphas_psi_ini);
Ekin_inf=0*1e3
Ekin_sup=600*1e3

radial_bin_size_half=6;
radial_bin_size=radial_bin_size_half*2;
psi_bin_pos=(1:radial_bin_size:255-radial_bin_size);
nb_psi_bins=length(psi_bin_pos)-1

for (psi=1:256)
    volume_radial(psi)=sum(volume_tor_diff(psi,:));
end
volume_radial_bin=zeros(nb_psi_bins,1);
for (psi=1:nb_psi_bins)
    volume_radial_bin(psi)=sum(volume_radial((psi-1)*radial_bin_size+1:psi*radial_bin_size));
end

[Npsi0 init_psi_bins]=histc(alphas_psi,psi_bin_pos);
density_radial_ini=Npsi0(1:end-1)./volume_radial_bin;
n_rescale=Ne0/max(density_radial_ini)
Eini_profile=density_radial_ini*0;

for (psi=1:nb_psi_bins-1)
    if ~isempty(init_psi_bins==psi)
        Eini_profile(psi)=mean(alphas_Ekin(init_psi_bins==psi));
    end
end

% for frame_movie=2000:2000:20000
    load(strcat('DT_MBall_collapse_Glisa_fc1h2_G110414','.mat'));
    [Npsi end_psi_bins]=histc(alphas_psi,psi_bin_pos);
    density_radial_end=Npsi(1:end-1)./volume_radial_bin;
    Eend_profile=density_radial_end*0;
    for (psi=1:nb_psi_bins-1)
        if ~isempty(end_psi_bins==psi)
            Eend_profile(psi)=mean(alphas_Ekin(end_psi_bins==psi));
        end
    end
    figure(1)
    clf(1);
    set(gca,'FontSize',26);
    hold on;
    grid on;
    plot(radial_r_value_flux(round(psi_bin_pos(1:end-1)+0.5*radial_bin_size)),3*(2/3)*n_rescale*eV*Eini_profile.*density_radial_ini,'r','LineWidth',3);
    plot(radial_r_value_flux(round(psi_bin_pos(1:end-1)+0.5*radial_bin_size)),3*(2/3)*n_rescale*eV*Eend_profile.*density_radial_end,'b--','LineWidth',3);
    ylabel('P (Pa)')
    xlim([0 1.8])
    xlabel('r (m)')

%     alphas_r=interp1(1:Nradial,radial_r_value_flux,alphas_psi);
%     figure(2)
%     clf(2);
%     set(gca,'FontSize',22)
%     title(num2str(frame_movie/(0.36*20000)))
%     grid on;
%     hold on
%     ENERGY_POP=(alphas_Ekin>=Ekin_inf).*(alphas_Ekin<=Ekin_sup).*(~alphas_ejected);
%     tr_parts=find(ENERGY_POP.*ALL_TRAPPED_POP);
%     cp_parts=find(ENERGY_POP.*CO_PASSING_POP);
%     cntp_parts=find(ENERGY_POP.*COUNTER_PASSING_POP);
%     plot(alphas_r_ini(tr_parts),alphas_r(tr_parts),'g.')
%     plot(alphas_r_ini(cp_parts),alphas_r(cp_parts),'r.')
% %     pause
%     plot(alphas_r_ini(cntp_parts),alphas_r(cntp_parts),'b.')
% %     plot(pphi_min:pphi_min,pphi_max:pphi_max,'k--','linewidth',2);
% %     plot(psi_min:psi_max,psi_min:psi_max,'k--','linewidth',2);
%     plot([0  2*r_value_q1_mean],[0  2*r_value_q1_mean],'k--','linewidth',2);
%     xlim([0  2*r_value_q1_mean])
%     ylim([0  2*r_value_q1_mean])
%     xlabel('r_{ini} (m)')
%     ylabel('r (m)')
%     pause
% end

Ne_final_kprofile=interp1(radial_r_value_flux(round(psi_bin_pos(1:end-1))),n_rescale*density_radial_end,radial_r_value_flux);
Te_final_kprofile=interp1(radial_r_value_flux(round(psi_bin_pos(1:end-1))),(2/3)*eV*Eend_profile,radial_r_value_flux);
P_final_kprofile=3*Ne_final_kprofile.*Te_final_kprofile;

save Pkinetic_final_profiles.mat Ne_final_kprofile Te_final_kprofile P_final_kprofile