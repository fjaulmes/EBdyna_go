filename=strcat(DATA_FOLDER,'q_profile.mat');
load(filename);
close all
load_particles_init_stats;

PART_POP=find(~alphas_ejected.*CO_PASSING_POP);

psiTAE=interp1(q_initial_profile,psi_scale,qTAE)
pphiTAE=(mHe/eV)*R0*vA_TAE-ZHe*psiTAE

tmax=time_stamp

Dpphi=pphi_output(tmax,:)-pphi_output(1,:);
% Dpphibis=pphi_evol(tmax,:)-pphi_evol(3,:);
DEkin=Ekin_output(tmax,:)-Ekin_output(2,:);
% DEkingc=Ekin_gc_output(tmax,:)-Ekin_gc_output(3,:);
DEpot=Epot_output(tmax,:)-Epot_output(1,:);
DE=Etot_output(tmax,:)-Etot_output(1,:);
DEbis=DEkin+ZHe*DEpot;
% DEth=alphas_Etot_th-Etot_output(2,:)';

RSIZE=0.05

radial_bins_lim=0.2:RSIZE:0.8;
radial_bins_pos=0.2+0.5*RSIZE:RSIZE:0.8-0.5*RSIZE;
DEkin_r=radial_bins_pos*0;
Dpphi_r=radial_bins_pos*0;

for r=1:length(radial_bins_pos)
    RPOP=find((r_avg>=radial_bins_lim(r)).*(r_avg<radial_bins_lim(r+1)));
    DEkin_r(r)=mean(DEkin(RPOP));
    Dpphi_r(r)=mean(Dpphi(RPOP));
end

%%
figure(2)
hold on
plot(radial_bins_pos,Dpphi_r)
plot([r1 r1],[-1.5 1.5],'g--','linewidth',3)
plot([r2 r2],[-1.5 1.5],'g--','linewidth',3)
plot([rTAE rTAE],[-1.5 1.5],'g--','linewidth',3)
% ylim([-1.2 1.2]*1e-3)


%%
figure(3)
set(gca,'fontsize',20)
hold on
grid on
plot(mean(vparallel_output(1:tmax,PART_POP),1),DE(PART_POP),'b.');
plot([vA_TAE/9 vA_TAE/9],[-4.5 4.5]*1e4,'g--','linewidth',3)
plot([vA_TAE/5 vA_TAE/5],[-4.5 4.5]*1e4,'g--','linewidth',3)
plot([vA3_TAE vA3_TAE],[-4.5 4.5]*1e4,'g--','linewidth',3)
plot([vA_TAE vA_TAE],[-4.5 4.5]*1e4,'g--','linewidth',3)
ylim([-4.5 4.5]*1e4)
xlim([0 2.0]*1e7)
xlabel('$$<v_{\parallel}>$$ (m/s)','Interpreter','latex')
ylabel('$$\Delta \mathcal{E}_{kin}$$ (eV)','Interpreter','latex')

%%
figure(4)
set(gca,'fontsize',24)
hold on
grid on

plot((-0.2:0.1:0.2),omega_TAE/nTAE*(-0.2:0.1:0.2),'k--','linewidth',3)
plot(Dpphi(PART_POP),DE(PART_POP),'g.')
plot((-0.2:0.1:0.2),omega_TAE/nTAE*(-0.2:0.1:0.2),'k--','linewidth',3)
ylabel('$$\Delta \mathcal{E}_{kin}$$ (eV)','Interpreter','latex')
xlabel('$$\Delta p_{\varphi}$$ (eV)','Interpreter','latex')
xlim([-0.2 0.2])

hl=legend('$$\Delta \mathcal{E}_{kin} = (\omega/n) \Delta p_{\varphi}$$');
set(hl,'Interpreter','latex')
% plot(Dpphi,DEth,'b.')
% plot(Dpphi,DE,'b.')
% plot(Dpphi,DEbis,'r.')
%%