close all

figure(1)
load('initialG_flatD_5keV_pre_collapse.mat')
load('initialG_flatD_5keV_precession_stats.mat')
alphas_psi_ini=alphas_psi;
alphas_pphi_ini=alphas_pphi0;
pphi_min=min(alphas_pphi_ini);
pphi_max=max(alphas_pphi_ini);
psi_min=min(alphas_psi_ini);
psi_max=max(alphas_psi_ini);
Ekin_inf=0*1e3
Ekin_sup=600*1e3

for frame_movie=2000:2000:20000
    load(strcat('flatD_5keV_collapse_Glisa_fc2h2_G290114record_',num2str(frame_movie),'.mat'));
    clf(1);
    set(gca,'FontSize',22)
    title(num2str(frame_movie/(0.36*20000)))
    grid on;
    hold on
    ENERGY_POP=(alphas_Ekin>=Ekin_inf).*(alphas_Ekin<=Ekin_sup).*(~alphas_ejected);
    tr_parts=find(ENERGY_POP.*ALL_TRAPPED_POP);
    cp_parts=find(ENERGY_POP.*CO_PASSING_POP);
    cntp_parts=find(ENERGY_POP.*COUNTER_PASSING_POP);
    plot(alphas_psi_ini(tr_parts),alphas_psi(tr_parts),'g.')
    plot(alphas_psi_ini(cp_parts),alphas_psi(cp_parts),'r.')
%     pause
    plot(alphas_psi_ini(cntp_parts),alphas_psi(cntp_parts),'b.')
%     plot(pphi_min:pphi_min,pphi_max:pphi_max,'k--','linewidth',2);
    plot(psi_min:psi_min,psi_max:psi_max,'k--','linewidth',2);
    pause
end