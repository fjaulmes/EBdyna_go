% load('initial_MB_H300_precession_stats_super_light.mat')
% load('initial_MB_H150_precession_stats_all.mat')
% load('initial_MB_H400_precession_stats_all.mat')

% load('initial_MBD_400_precession_stats_all.mat')
% load('initial_MBD_400_pre_collapse_all.mat')
load('initial_MBD_250_precession_stats_all.mat')
load('initial_MBD_250_pre_collapse_all.mat')

% load('initial_MBD_300_precession_stats_super_light.mat')
% load('initial_MBD_300_pre_collapse_super_light.mat')

Nalphas_simulated=length(alphas_weight)

ALL_TRAPPED=find(ALL_TRAPPED_POP(1:Nalphas_simulated));
CO_PASSING=find(CO_PASSING_POP(1:Nalphas_simulated));
COUNTER_PASSING=find(COUNTER_PASSING_POP(1:Nalphas_simulated));
STAGNATION=find(STAGNATION_POP(1:Nalphas_simulated));
disp('---------------------------------------------------')

trapped_fraction=sum(alphas_weight(ALL_TRAPPED))/sum(alphas_weight)
copassing_fraction=sum(alphas_weight(CO_PASSING))/sum(alphas_weight)
counterpassing_fraction=sum(alphas_weight(COUNTER_PASSING))/sum(alphas_weight)
stagnation_fraction=sum(alphas_weight(STAGNATION))/sum(alphas_weight)

disp('---------------------------------------------------')
alphas_Bfield_ini=interp2(scale_X,scale_Z,Btot_XZ_map',pos_X_gc,pos_Z_gc);
alphas_Eperp=max(alphas_Ekin-0.5*(mHe/eV)*alphas_vpll.^2,0);
alphas_vperp=sqrt(2*(eV/mHe)*alphas_Eperp);
alphas_omega_ci=ZHe*eV*alphas_Bfield_ini/mHe;
alphas_rhoL=alphas_vperp./alphas_omega_ci;

average_rhoL=mean(alphas_weight.*alphas_rhoL)
average_deltar=mean(alphas_weight.*delta_r_avg)


disp('---------------------------------------------------')
% load('initial_MBD_500_pre_collapse_super_light.mat');
% load('initial_MBD_500_precession_stats_super_light.mat');

%%
NB_PART_RESCALE=3.96e12

mHem=mD
ZHe=1
tmax=time_stamp
TINI=2
POP_NON_EJ=find(~(alphas_ejected));
POP_PLOT=POP_NON_EJ;
% POP_PLOT=find(~(alphas_ejected).*STAGNATION_POP(1:length(alphas_ejected)));
% POP_PLOT=find(~(alphas_ejected).*ALL_PASSING_POP(1:length(alphas_ejected)));
% POP_PLOT=find(~(alphas_ejected).*ALL_TRAPPED_POP(1:length(alphas_ejected)));
% POP_PLOT=find(~(alphas_ejected).*ALL_TRAPPED_POP(1:length(alphas_ejected)));

pphi_recalc_output=(mHe/eV)*(R0+Xpos_part_output(:,:)).*vphi_output(:,:)-ZHe*psi_value_output(:,:);
Dpphi_recalc=pphi_recalc_output(tmax,:)-pphi_recalc_output(TINI,:);
Dpphi=pphi_output(tmax,:)-pphi_output(TINI,:);
DE=Etot_output(tmax,:)-Etot_output(TINI,:);
DE_gc=Ekin_gc_output(tmax,:)-Ekin_gc_output(TINI,:);
DE_gc_recalc=Ekin_gc_output(tmax,:)+Epot_output(tmax,:)-Ekin_gc_output(TINI,:)-Epot_output(TINI,:);

% DE=DE_gc;
alphas_weight=alphas_weight(1:length(DE));
alphas_weight=alphas_weight';
% DE_th=Etot_th_output(tmax,:)-Etot_th_output(TINI,:);

DE_straight=(omega_TAE/nTAE)*Dpphi;



% DE=DE_gc;

P1=POP_PLOT(1);
P2=POP_PLOT(10);
P3=POP_PLOT(20);
P4=POP_PLOT(30);

P5=100;
P6=150;
P7=200;
P8=300;

figure(1)
subplot(2,1,1)
grid on
hold on 
plot(pphi_output(:,P1))
plot(pphi_recalc_output(:,P1),'r--')
xlim([0 tmax])

subplot(2,1,2)
grid on
hold on 
plot(pphi_output(:,P2))
plot(pphi_recalc_output(:,P2),'r--')
xlim([0 tmax])



figure(2)
subplot(2,1,1)
grid on
hold on 
plot(pphi_output(:,P3))
plot(pphi_recalc_output(:,P3),'r--')
xlim([0 tmax])

subplot(2,1,2)
grid on
hold on 
plot(pphi_output(:,P4))
plot(pphi_recalc_output(:,P4),'r--')
xlim([0 tmax])


figure(3)
subplot(2,1,1)
grid on
hold on 
plot(pphi_output(:,P5))
plot(pphi_recalc_output(:,P5),'r--')
xlim([0 tmax])

subplot(2,1,2)
grid on
hold on 
plot(pphi_output(:,P6))
plot(pphi_recalc_output(:,P6),'r--')
xlim([0 tmax])


figure(4)
subplot(2,1,1)
grid on
hold on 
plot(pphi_output(:,P7))
plot(pphi_recalc_output(:,P7),'r--')
xlim([0 tmax])

subplot(2,1,2)
grid on
hold on 
plot(pphi_output(:,P8))
plot(pphi_recalc_output(:,P8),'r--')
xlim([0 tmax])

%%
POP_NON_EJ=find(~(alphas_ejected));
POP_PLOT_PASSING=find(~(alphas_ejected).*ALL_PASSING_POP(1:length(alphas_ejected)));
POP_PLOT_COPASSING=find(~(alphas_ejected).*CO_PASSING_POP(1:length(alphas_ejected)));
POP_PLOT_COUNTERPASSING=find(~(alphas_ejected).*COUNTER_PASSING_POP(1:length(alphas_ejected)));
POP_PLOT_TRAPPED=find(~(alphas_ejected).*ALL_TRAPPED_POP(1:length(alphas_ejected)));
POP_PLOT_STAGNATION=find(~(alphas_ejected).*STAGNATION_POP(1:length(alphas_ejected)));
POP_PLOT_POTATOES=find(~(alphas_ejected).*POTATOES_POP(1:length(alphas_ejected)));

figure(8)
set(gca,'fontsize',20)
hold on
grid on
% plot(Dpphi(POP_PLOT),DE_gc(POP_PLOT),'g.')
plot(Dpphi_recalc(POP_PLOT_TRAPPED(1:8:end)),DE_gc_recalc(POP_PLOT_TRAPPED(1:8:end)),'g.')
plot(Dpphi_recalc(POP_PLOT_POTATOES(1:8:end)),DE_gc_recalc(POP_PLOT_POTATOES(1:8:end)),'k.')
plot(Dpphi_recalc(POP_PLOT_COPASSING(1:8:end)),DE_gc_recalc(POP_PLOT_COPASSING(1:8:end)),'b.')
plot(Dpphi_recalc(POP_PLOT_COUNTERPASSING(1:8:end)),DE_gc_recalc(POP_PLOT_COUNTERPASSING(1:8:end)),'.','color',[0 0 0.5])
plot(Dpphi_recalc(POP_PLOT_STAGNATION(1:8:end)),DE_gc_recalc(POP_PLOT_STAGNATION(1:8:end)),'r.')


% deltaE_potatoes=length(POP_PLOT_POTATOES)*mean(alphas_weight(POP_PLOT_POTATOES).*DE(POP_PLOT_POTATOES))
% deltaE_stagnation=length(POP_PLOT_STAGNATION)*mean(alphas_weight(POP_PLOT_STAGNATION).*DE(POP_PLOT_STAGNATION))
% deltaE_passing=length(POP_PLOT_PASSING)*mean(alphas_weight(POP_PLOT_PASSING).*DE(POP_PLOT_PASSING))
% deltaE_copassing=length(POP_PLOT_COPASSING)*mean(alphas_weight(POP_PLOT_COPASSING).*DE(POP_PLOT_COPASSING))
% deltaE_counterpassing=length(POP_PLOT_COUNTERPASSING)*mean(alphas_weight(POP_PLOT_COUNTERPASSING).*DE(POP_PLOT_COUNTERPASSING))
% deltaE_trapped=length(POP_PLOT_TRAPPED)*mean(alphas_weight(POP_PLOT_TRAPPED).*DE(POP_PLOT_TRAPPED))

deltaE_potatoes=eV*NB_PART_RESCALE*sum(alphas_weight(POP_PLOT_POTATOES).*DE(POP_PLOT_POTATOES))
deltaE_stagnation=eV*NB_PART_RESCALE*sum(alphas_weight(POP_PLOT_STAGNATION).*DE(POP_PLOT_STAGNATION))
deltaE_passing=eV*NB_PART_RESCALE*sum(alphas_weight(POP_PLOT_PASSING).*DE(POP_PLOT_PASSING))
deltaE_counterpassing=eV*NB_PART_RESCALE*sum(alphas_weight(POP_PLOT_COUNTERPASSING).*DE(POP_PLOT_COUNTERPASSING))
deltaE_copassing=eV*NB_PART_RESCALE*sum(alphas_weight(POP_PLOT_COPASSING).*DE(POP_PLOT_COPASSING))
deltaE_trapped=eV*NB_PART_RESCALE*sum(alphas_weight(POP_PLOT_TRAPPED).*DE(POP_PLOT_TRAPPED))

% deltaE_potatoes =
%   -7.1722e-04
% deltaE_stagnation =
%    -5.7518
% deltaE_passing =
%   -55.3146
% deltaE_counterpassing =
%   -10.6643
% deltaE_copassing =
%   -44.6503
% deltaE_trapped =
%     0.5510
    
% plot(Dpphi(POP_PLOT),DE(POP_PLOT),'r.')
plot([-0.22 0.22],(omega_TAE/nTAE)*[-0.22 0.22],'k--')
% plot(Dpphi_recalc(POP_PLOT),DE_gc(POP_PLOT),'k.')
xl=xlabel('$\Delta p_\varphi /e$')
set(xl,'Interpreter','latex')
yl=ylabel('$\Delta \mathcal{E} (eV)$')
set(yl,'Interpreter','latex')
xlim([-0.22 0.22])
ylim([-1.5 1.5]*1e4)

hl=legend(strcat('$\Delta \mathcal{E} _{\rm{trapped}} = 391 \ \rm{J}$'),...
    strcat('$\Delta \mathcal{E} _{\rm{potatoes}} = -0.005 \ \rm{J}$'),...
    strcat('$\Delta \mathcal{E} _{\rm{counter-passing}} = -448 \ \rm{J}$'),...
    strcat('$\Delta \mathcal{E} _{\rm{co-passing}} = -1522 \ \rm{J}$'),...
    strcat('$\Delta \mathcal{E} _{\rm{stagnation}} = -507 \ \rm{J}$'),...
    '$\Delta \mathcal{E} = (\omega_{\rm{TAE}}/n )\Delta p_\varphi $');
set(hl,'Interpreter','latex')
set(hl,'fontsize',18)

%%
VPLLBINSIZE=0.04
vpll_bins=(-1.6:VPLLBINSIZE:1.6)*1e7;
vpll_values=vpll_bins(1:end-1)+0.5*VPLLBINSIZE*1e7;

DE_vpll_values_trapped=vpll_values*0;
DE_vpll_values_potatoes=vpll_values*0;
DE_vpll_values_copassing=vpll_values*0;
DE_vpll_values_counterpassing=vpll_values*0;
DE_vpll_values_stagnation=vpll_values*0;

for vb=1:length(vpll_values)
    POP_VPLL=find(~(alphas_ejected).*COUNTER_PASSING_POP(1:length(alphas_ejected)).*(vpll_avg(1:length(alphas_ejected))>vpll_bins(vb)).*(vpll_avg(1:length(alphas_ejected))<=vpll_bins(vb+1)));
    DE_vpll_values_counterpassing(vb)=sum(DE(POP_VPLL).*alphas_weight(POP_VPLL));
    POP_VPLL=find(~(alphas_ejected).*CO_PASSING_POP(1:length(alphas_ejected)).*(vpll_avg(1:length(alphas_ejected))>vpll_bins(vb)).*(vpll_avg(1:length(alphas_ejected))<=vpll_bins(vb+1)));
    DE_vpll_values_copassing(vb)=sum(DE(POP_VPLL).*alphas_weight(POP_VPLL));

    POP_VPLL=find(~(alphas_ejected).*STAGNATION_POP(1:length(alphas_ejected)).*(vpll_avg(1:length(alphas_ejected))>vpll_bins(vb)).*(vpll_avg(1:length(alphas_ejected))<=vpll_bins(vb+1)));
    DE_vpll_values_stagnation(vb)=sum(DE(POP_VPLL).*alphas_weight(POP_VPLL));
    POP_VPLL=find(~(alphas_ejected).*ALL_TRAPPED_POP(1:length(alphas_ejected)).*(vpll_avg(1:length(alphas_ejected))>vpll_bins(vb)).*(vpll_avg(1:length(alphas_ejected))<=vpll_bins(vb+1)));
    DE_vpll_values_trapped(vb)=sum(DE(POP_VPLL).*alphas_weight(POP_VPLL));
    POP_VPLL=find(~(alphas_ejected).*POTATOES_POP(1:length(alphas_ejected)).*(vpll_avg(1:length(alphas_ejected))>vpll_bins(vb)).*(vpll_avg(1:length(alphas_ejected))<=vpll_bins(vb+1)));
    DE_vpll_values_potatoes(vb)=sum(DE(POP_VPLL).*alphas_weight(POP_VPLL));

end

DE_vpll_values_counterpassing=DE_vpll_values_counterpassing/max(DE_vpll_values_trapped);
DE_vpll_values_copassing=DE_vpll_values_copassing/max(DE_vpll_values_trapped);
DE_vpll_values_stagnation=DE_vpll_values_stagnation/max(DE_vpll_values_trapped);
DE_vpll_values_trapped=DE_vpll_values_trapped/max(DE_vpll_values_trapped);
DE_vpll_values_potatoes=DE_vpll_values_potatoes/max(DE_vpll_values_trapped);

%%
figure1=figure(7)
xlim([-1.5 1.5]*1e7)
ylim([-2 2]*1e4)

set(gca,'fontsize',20)
grid on
hold on

% plot(mean(vparallel_output(1:tmax,POP_PLOT),1),DE(POP_PLOT),'b.')
plot(vpll_avg(POP_PLOT),DE(POP_PLOT),'b.')
Epll_output=0.5*(mHe/eV)*vparallel_output.^2;
Eperp_output=max(Ekin_output-Epll_output,0);
vperp_output=sqrt(2*(eV/mHe)*Eperp_output);

plot([vA_TAE/5 vA_TAE/5],[-4.5 4.5]*1e4,'g--','linewidth',3)
plot([vA3_TAE vA3_TAE],[-4.5 4.5]*1e4,'g--','linewidth',3)
plot([vA_TAE vA_TAE],[-4.5 4.5]*1e4,'g--','linewidth',3)
plot(-[vA_TAE/5 vA_TAE/5],[-4.5 4.5]*1e4,'g--','linewidth',3)
plot(-[vA3_TAE vA3_TAE],[-4.5 4.5]*1e4,'g--','linewidth',3)
plot(-[vA_TAE vA_TAE],[-4.5 4.5]*1e4,'g--','linewidth',3)

xl=xlabel('$\left< v_{\parallel} \right> (\rm{m/s})$')
set(xl,'Interpreter','latex')
yl=ylabel('$\Delta \mathcal{E} (a.u.)$')
set(yl,'Interpreter','latex')


% Create textbox
annotation(figure1,'textbox',...
    [0.845006568144499 0.848101265822784 0.0704285714285715 0.0770042194092819],...
    'String',{'v_A'},...
    'FontSize',36,...
    'FontName','Agency FB',...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'Color',[0 1 0]);

% Create textbox
annotation(figure1,'textbox',...
    [0.524579638752052 0.846075949367085 0.0704285714285715 0.0770042194092819],...
    'String',{'v_A/5'},...
    'FontSize',36,...
    'FontName','Agency FB',...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'Color',[0 1 0]);

% Create textbox
annotation(figure1,'textbox',...
    [0.6457947454844 0.85236286919831 0.0704285714285715 0.0770042194092819],...
    'String',{'v_A/3'},...
    'FontSize',36,...
    'FontName','Agency FB',...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'Color',[0 1 0]);

% Create textbox
annotation(figure1,'textbox',...
    [0.314694581280788 0.173122362869194 0.0704285714285715 0.0770042194092819],...
    'String',{'-v_A/3'},...
    'FontSize',36,...
    'FontName','Agency FB',...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'Color',[0 1 0]);

% Create textbox
annotation(figure1,'textbox',...
    [0.44083579638752 0.159367088607591 0.0704285714285715 0.0770042194092819],...
    'String',{'-v_A/5'},...
    'FontSize',36,...
    'FontName','Agency FB',...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'Color',[0 1 0]);

% Create textbox
annotation(figure1,'textbox',...
    [0.150720853858785 0.327046413502108 0.0704285714285715 0.0770042194092819],...
    'String',{'-v_A'},...
    'FontSize',36,...
    'FontName','Agency FB',...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'Color',[0 1 0]);



%%
figure(6)
set(gca,'fontsize',20)

hold on
% grid on
plot(vpll_values,DE_vpll_values_trapped,'k','linewidth',2)
% plot(vpll_values,DE_vpll_values_potatoes,'k','linewidth',2)
plot(vpll_values,DE_vpll_values_counterpassing,'color',[0 0 0.7],'linewidth',2)
plot(vpll_values,DE_vpll_values_copassing,'b','linewidth',2)
plot(vpll_values,DE_vpll_values_stagnation,'r','linewidth',2)

plot([vA_TAE/5 vA_TAE/5],[-4.5 4.5]*1e6,'g--','linewidth',3)
plot([vA3_TAE vA3_TAE],[-4.5 4.5]*1e6,'g--','linewidth',3)
plot([vA_TAE vA_TAE],[-4.5 4.5]*1e6,'g--','linewidth',3)
plot(-[vA_TAE/5 vA_TAE/5],[-4.5 4.5]*1e6,'g--','linewidth',3)
plot(-[vA3_TAE vA3_TAE],[-4.5 4.5]*1e6,'g--','linewidth',3)
plot(-[vA_TAE vA_TAE],[-4.5 4.5]*1e6,'g--','linewidth',3)


plot(vpll_values,DE_vpll_values_trapped,'k','linewidth',2)
% plot(vpll_values,DE_vpll_values_potatoes,'k','linewidth',2)
plot(vpll_values,DE_vpll_values_counterpassing,'color',[0 0 0.7],'linewidth',2)
plot(vpll_values,DE_vpll_values_copassing,'b','linewidth',2)
plot(vpll_values,DE_vpll_values_stagnation,'r','linewidth',2)

legend('trapped','counter-passing','co-passing','stagnation')

xlim([-1.5 1.5]*1e7)
ylim([-18.2 3.2])

xl=xlabel('$\left< v_{\parallel} \right> (\rm{m/s})$')
set(xl,'Interpreter','latex')
yl=ylabel('$\Delta \mathcal{E} (a.u.)$')
set(yl,'Interpreter','latex')

% figure(6)
% plot(mean(Ekin_output(1:tmax,POP_PLOT),1),DE(POP_PLOT),'b.')

%%
alphas_Eperp=max(alphas_Ekin-0.5*(mHe/eV)*vpll_avg.^2,0);
alphas_vperp=sqrt(2*(eV/mHe)*alphas_Eperp);

VPERPBINSIZE=0.05
vperp_bins=(0:VPERPBINSIZE:2)*1e7;
vperp_values=vperp_bins(1:end-1)+0.5*VPERPBINSIZE*1e7;
DE_vperp_values_trapped=vperp_values*0;
DE_vperp_values_potatoes=vperp_values*0;
DE_vperp_values_passing=vperp_values*0;
DE_vperp_values_stagnation=vperp_values*0;

for vb=1:length(vperp_values)
    POP_VPLL=find(~(alphas_ejected).*ALL_PASSING_POP(1:length(alphas_ejected)).*(alphas_vperp(1:length(alphas_ejected))>vperp_bins(vb)).*(alphas_vperp(1:length(alphas_ejected))<=vperp_bins(vb+1)));
    DE_vperp_values_passing(vb)=sum(DE(POP_VPLL).*alphas_weight(POP_VPLL));
    POP_VPLL=find(~(alphas_ejected).*STAGNATION_POP(1:length(alphas_ejected)).*(alphas_vperp(1:length(alphas_ejected))>vperp_bins(vb)).*(alphas_vperp(1:length(alphas_ejected))<=vperp_bins(vb+1)));
    DE_vperp_values_stagnation(vb)=sum(DE(POP_VPLL).*alphas_weight(POP_VPLL));
    POP_VPLL=find(~(alphas_ejected).*ALL_TRAPPED_POP(1:length(alphas_ejected)).*(alphas_vperp(1:length(alphas_ejected))>vperp_bins(vb)).*(alphas_vperp(1:length(alphas_ejected))<=vperp_bins(vb+1)));
    DE_vperp_values_trapped(vb)=sum(DE(POP_VPLL).*alphas_weight(POP_VPLL));
    POP_VPLL=find(~(alphas_ejected).*POTATOES_POP(1:length(alphas_ejected)).*(alphas_vperp(1:length(alphas_ejected))>vperp_bins(vb)).*(alphas_vperp(1:length(alphas_ejected))<=vperp_bins(vb+1)));
    DE_vperp_values_potatoes(vb)=sum(DE(POP_VPLL).*alphas_weight(POP_VPLL));

end

DE_vperp_values_passing=DE_vperp_values_passing/max(DE_vperp_values_trapped);
DE_vperp_values_stagnation=DE_vperp_values_stagnation/max(DE_vperp_values_trapped);
DE_vperp_values_trapped=DE_vperp_values_trapped/max(DE_vperp_values_trapped);
DE_vperp_values_potatoes=DE_vperp_values_potatoes/max(DE_vperp_values_trapped);


%%
figure(5)
set(gca,'fontsize',20)

hold on

hold on
grid on
plot(vperp_values,DE_vperp_values_trapped,'k','linewidth',2)
% plot(vpll_values,DE_vperp_values_potatoes,'k','linewidth',2)
plot(vperp_values,DE_vperp_values_passing,'b','linewidth',2)
plot(vperp_values,DE_vperp_values_stagnation,'r','linewidth',2)

plot([vA_TAE/5 vA_TAE/5],[-4.5 4.5],'g--','linewidth',3)
plot([vA3_TAE vA3_TAE],[-4.5 4.5],'g--','linewidth',3)
plot([vA_TAE vA_TAE],[-4.5 4.5],'g--','linewidth',3)


plot(vperp_values,DE_vperp_values_trapped,'k','linewidth',2)
% plot(vpll_values,DE_vpll_values_potatoes,'k','linewidth',2)
plot(vperp_values,DE_vperp_values_passing,'b','linewidth',2)
plot(vperp_values,DE_vperp_values_stagnation,'r','linewidth',2)

legend('trapped','passing','stagnation')

ylim([-4.0 1.05])

xl=xlabel('$\left< v_{\bot} \right> (\rm{m/s})$')
set(xl,'Interpreter','latex')
yl=ylabel('$\Delta \mathcal{E} (a.u.)$')
set(yl,'Interpreter','latex')





%%
VPLLBINSIZE=0.04
vpll_bins=(-1.6:VPLLBINSIZE:1.6)*1e7;
vpll_values=vpll_bins(1:end-1)+0.5*VPLLBINSIZE*1e7;

MMBINSIZE=0.15
mm_bins=(0:MMBINSIZE:3.0)*1e5;
mm_values=mm_bins(1:end-1)+0.5*MMBINSIZE*1e5;
DE_vspace_values_trapped=zeros(length(vpll_values),length(mm_values));
DE_vspace_values_potatoes=zeros(length(vpll_values),length(mm_values));
DE_vspce_values_copassing=zeros(length(vpll_values),length(mm_values));
DE_vspace_values_counterpassing=zeros(length(vpll_values),length(mm_values));
DE_vspace_values_stagnation=zeros(length(vpll_values),length(mm_values));
DE_vspace_values_all=zeros(length(vpll_values),length(mm_values));

for vb=1:length(vpll_values)
    for mmb=1:length(mm_values)
        POP_VPLL=find(~(alphas_ejected)...
            .*(vpll_avg(1:length(alphas_ejected))>vpll_bins(vb)).*(vpll_avg(1:length(alphas_ejected))<=vpll_bins(vb+1))...
            .*(alphas_mm(1:length(alphas_ejected))>mm_bins(mmb)).*(alphas_mm(1:length(alphas_ejected))<=mm_bins(mmb+1)));
%         DE_vpll_values_all(vb,mmb)=sum(DE(POP_VPLL).*alphas_weight(POP_VPLL));
        DE_vspace_values_all(vb,mmb)=sum(DE(POP_VPLL).*alphas_weight(POP_VPLL));
%         POP_VPLL=find(~(alphas_ejected).*COUNTER_PASSING_POP(1:length(alphas_ejected))...
%             .*(vpll_avg(1:length(alphas_ejected))>vpll_bins(vb)).*(vpll_avg(1:length(alphas_ejected))<=vpll_bins(vb+1))...
%             .*(alphas_mm(1:length(alphas_ejected))>mm_bins(mmb)).*(alphas_mm(1:length(alphas_ejected))<=mm_bins(mmb+1)));
%         DE_vpll_values_counterpassing(vb,mmb)=sum(DE(POP_VPLL).*alphas_weight(POP_VPLL));
%         POP_VPLL=find(~(alphas_ejected).*CO_PASSING_POP(1:length(alphas_ejected))...
%             .*(vpll_avg(1:length(alphas_ejected))>vpll_bins(vb)).*(vpll_avg(1:length(alphas_ejected))<=vpll_bins(vb+1))...
%             .*(alphas_mm(1:length(alphas_ejected))>mm_bins(mmb)).*(alphas_mm(1:length(alphas_ejected))<=mm_bins(mmb+1)));
%         DE_vpll_values_copassing(vb,mmb)=sum(DE(POP_VPLL).*alphas_weight(POP_VPLL));
        
%         POP_VPLL=find(~(alphas_ejected).*STAGNATION_POP(1:length(alphas_ejected))...
%             .*(vpll_avg(1:length(alphas_ejected))>vpll_bins(vb)).*(vpll_avg(1:length(alphas_ejected))<=vpll_bins(vb+1))...
%             .*(alphas_mm(1:length(alphas_ejected))>mm_bins(mmb)).*(alphas_mm(1:length(alphas_ejected))<=mm_bins(mmb+1)));
%         DE_vpll_values_stagnation(vb,mmb)=sum(DE(POP_VPLL).*alphas_weight(POP_VPLL));
%         POP_VPLL=find(~(alphas_ejected).*ALL_TRAPPED_POP(1:length(alphas_ejected))...
%             .*(vpll_avg(1:length(alphas_ejected))>vpll_bins(vb)).*(vpll_avg(1:length(alphas_ejected))<=vpll_bins(vb+1))...
%             .*(alphas_mm(1:length(alphas_ejected))>mm_bins(mmb)).*(alphas_mm(1:length(alphas_ejected))<=mm_bins(mmb+1)));
%         DE_vpll_values_trapped(vb,mmb)=sum(DE(POP_VPLL).*alphas_weight(POP_VPLL));
%         POP_VPLL=find(~(alphas_ejected).*POTATOES_POP(1:length(alphas_ejected))...
%             .*(vpll_avg(1:length(alphas_ejected))>vpll_bins(vb)).*(vpll_avg(1:length(alphas_ejected))<=vpll_bins(vb+1))...
%             .*(alphas_mm(1:length(alphas_ejected))>mm_bins(mmb)).*(alphas_mm(1:length(alphas_ejected))<=mm_bins(mmb+1)));
%         DE_vpll_values_potatoes(vb,mmb)=sum(DE(POP_VPLL).*alphas_weight(POP_VPLL));
    end
end

%%
figure(11)
set(gca,'fontsize',20)
hold on
imagesc(vpll_values,mm_values,-DE_vspace_values_all');
axis xy
plot([vA_TAE/5 vA_TAE/5],[-4.5 4.5]*1e6,'k--','linewidth',3)
plot([vA3_TAE vA3_TAE],[-4.5 4.5]*1e6,'k--','linewidth',3)
plot([vA_TAE vA_TAE],[-4.5 4.5]*1e6,'k--','linewidth',3)
plot(-[vA_TAE/5 vA_TAE/5],[-4.5 4.5]*1e6,'k--','linewidth',3)
plot(-[vA3_TAE vA3_TAE],[-4.5 4.5]*1e6,'k--','linewidth',3)
plot(-[vA_TAE vA_TAE],[-4.5 4.5]*1e6,'k--','linewidth',3)


xlim([-1.05 1.05]*1e7)
ylim([0 3]*1e5)

ylabel('\mu (eV/T)')
xl=xlabel('$\left< v_{\parallel} \right> (\rm{m/s})$')
set(xl,'Interpreter','latex')


%%

DRBINSIZE=0.035
dr_bins=(0:DRBINSIZE:0.7);
dr_values=dr_bins(1:end-1)+0.5*DRBINSIZE;
DE_dr_values_trapped=dr_values*0;
DE_dr_values_potatoes=dr_values*0;
DE_dr_values_counterpassing=dr_values*0;
DE_dr_values_copassing=dr_values*0;
DE_dr_values_stagnation=dr_values*0;

for vb=1:length(dr_values)
    POP_VPLL=find(~(alphas_ejected).*COUNTER_PASSING_POP(1:length(alphas_ejected))...
        .*(delta_r_avg(1:length(alphas_ejected))>dr_bins(vb)).*(delta_r_avg(1:length(alphas_ejected))<=dr_bins(vb+1)));
    DE_dr_values_counterpassing(vb)=sum(DE(POP_VPLL).*alphas_weight(POP_VPLL));
    POP_VPLL=find(~(alphas_ejected).*CO_PASSING_POP(1:length(alphas_ejected))...
        .*(delta_r_avg(1:length(alphas_ejected))>dr_bins(vb)).*(delta_r_avg(1:length(alphas_ejected))<=dr_bins(vb+1)));
    DE_dr_values_copassing(vb)=sum(DE(POP_VPLL).*alphas_weight(POP_VPLL));
    POP_VPLL=find(~(alphas_ejected).*ALL_TRAPPED_POP(1:length(alphas_ejected))...
        .*(delta_r_avg(1:length(alphas_ejected))>dr_bins(vb)).*(delta_r_avg(1:length(alphas_ejected))<=dr_bins(vb+1)));
    DE_dr_values_trapped(vb)=sum(DE(POP_VPLL).*alphas_weight(POP_VPLL));
    POP_VPLL=find(~(alphas_ejected).*STAGNATION_POP(1:length(alphas_ejected))...
        .*(delta_r_avg(1:length(alphas_ejected))>dr_bins(vb)).*(delta_r_avg(1:length(alphas_ejected))<=dr_bins(vb+1)));
    DE_dr_values_stagnation(vb)=sum(DE(POP_VPLL).*alphas_weight(POP_VPLL));
%     POP_VPLL=find(~(alphas_ejected).*COUNTER_PASSING_POP(1:length(alphas_ejected))...
%         .*(delta_r_avg(1:length(alphas_ejected))>dr_bins(vb)).*(delta_r_avg(1:length(alphas_ejected))<=dr_bins(vb+1)));
%     DE_dr_values_counterpassing(vb)=sum(DE(POP_VPLL));
%     POP_VPLL=find(~(alphas_ejected).*CO_PASSING_POP(1:length(alphas_ejected))...
%         .*(delta_r_avg(1:length(alphas_ejected))>dr_bins(vb)).*(delta_r_avg(1:length(alphas_ejected))<=dr_bins(vb+1)));
%     DE_dr_values_copassing(vb)=sum(DE(POP_VPLL));
%     POP_VPLL=find(~(alphas_ejected).*ALL_TRAPPED_POP(1:length(alphas_ejected))...
%         .*(delta_r_avg(1:length(alphas_ejected))>dr_bins(vb)).*(delta_r_avg(1:length(alphas_ejected))<=dr_bins(vb+1)));
%     DE_dr_values_trapped(vb)=sum(DE(POP_VPLL));
%     POP_VPLL=find(~(alphas_ejected).*STAGNATION_POP(1:length(alphas_ejected))...
%         .*(delta_r_avg(1:length(alphas_ejected))>dr_bins(vb)).*(delta_r_avg(1:length(alphas_ejected))<=dr_bins(vb+1)));
%     DE_dr_values_stagnation(vb)=sum(DE(POP_VPLL));
end

DE_dr_values_stagnation=DE_dr_values_stagnation/max(DE_dr_values_trapped);
DE_dr_values_counterpassing=DE_dr_values_counterpassing/max(DE_dr_values_trapped);
DE_dr_values_copassing=DE_dr_values_copassing/max(DE_dr_values_trapped);
DE_dr_values_trapped=DE_dr_values_trapped/max(DE_dr_values_trapped);







%%
DRBINSIZE=0.035
dr_bins=(0:DRBINSIZE:0.7);
dr_values=dr_bins(1:end-1)+0.5*DRBINSIZE;

figure(15)
set(gca,'fontsize',20)

hold on

hold on
grid on
plot(dr_values,DE_dr_values_trapped,'k','linewidth',2)
plot(dr_values,DE_dr_values_counterpassing,'--','color',[0 0 0.7],'linewidth',2)
plot(dr_values,DE_dr_values_copassing,'b','linewidth',2)
% plot(dr_values,DE_dr_values_stagnation,'r','linewidth',2)

xl=xlabel('$\delta r (\rm{m})$')
set(xl,'Interpreter','latex')
yl=ylabel('$\Delta \mathcal{E} (a.u.)$')
set(yl,'Interpreter','latex')
% ylim([-4.0 1.05])

legend('trapped','counter-passing','co-passing')
% legend('trapped','counter-passing','co-passing','stagnation')

xlim([0 0.66])






%%

alphas_Bfield_ini=interp2(scale_X,scale_Z,Btot_XZ_map',pos_X_gc,pos_Z_gc);

alphas_Eperp=max(alphas_Ekin-0.5*(mHe/eV)*alphas_vpll.^2,0);
alphas_vperp=sqrt(2*(eV/mHe)*alphas_Eperp);
alphas_omega_ci=ZHe*eV*alphas_Bfield_ini/mHe;
alphas_rhoL=alphas_vperp./alphas_omega_ci;


DRBINSIZE=0.006
dr_bins=(0:DRBINSIZE:0.08);
dr_values=dr_bins(1:end-1)+0.5*DRBINSIZE;
DE_rhoL_values_trapped=dr_values*0;
DE_rhoL_values_potatoes=dr_values*0;
DE_rhoL_values_counterpassing=dr_values*0;
DE_rhoL_values_copassing=dr_values*0;
DE_rhoL_values_stagnation=dr_values*0;

for vb=1:length(dr_values)
    POP_VPLL=find(~(alphas_ejected).*COUNTER_PASSING_POP(1:length(alphas_ejected))...
        .*(alphas_rhoL(1:length(alphas_ejected))>dr_bins(vb)).*(alphas_rhoL(1:length(alphas_ejected))<=dr_bins(vb+1)));
    DE_rhoL_values_counterpassing(vb)=sum(DE(POP_VPLL).*alphas_weight(POP_VPLL));
    POP_VPLL=find(~(alphas_ejected).*CO_PASSING_POP(1:length(alphas_ejected))...
        .*(alphas_rhoL(1:length(alphas_ejected))>dr_bins(vb)).*(alphas_rhoL(1:length(alphas_ejected))<=dr_bins(vb+1)));
    DE_rhoL_values_copassing(vb)=sum(DE(POP_VPLL).*alphas_weight(POP_VPLL));
    POP_VPLL=find(~(alphas_ejected).*ALL_TRAPPED_POP(1:length(alphas_ejected))...
        .*(alphas_rhoL(1:length(alphas_ejected))>dr_bins(vb)).*(alphas_rhoL(1:length(alphas_ejected))<=dr_bins(vb+1)));
    DE_rhoL_values_trapped(vb)=sum(DE(POP_VPLL).*alphas_weight(POP_VPLL));
    POP_VPLL=find(~(alphas_ejected).*STAGNATION_POP(1:length(alphas_ejected))...
        .*(alphas_rhoL(1:length(alphas_ejected))>dr_bins(vb)).*(alphas_rhoL(1:length(alphas_ejected))<=dr_bins(vb+1)));
    DE_rhoL_values_stagnation(vb)=sum(DE(POP_VPLL).*alphas_weight(POP_VPLL));
%     POP_VPLL=find(~(alphas_ejected).*COUNTER_PASSING_POP(1:length(alphas_ejected))...
%         .*(alphas_rhoL(1:length(alphas_ejected))>dr_bins(vb)).*(alphas_rhoL(1:length(alphas_ejected))<=dr_bins(vb+1)));
%     DE_rhoL_values_counterpassing(vb)=sum(DE(POP_VPLL));
%     POP_VPLL=find(~(alphas_ejected).*CO_PASSING_POP(1:length(alphas_ejected))...
%         .*(alphas_rhoL(1:length(alphas_ejected))>dr_bins(vb)).*(alphas_rhoL(1:length(alphas_ejected))<=dr_bins(vb+1)));
%     DE_rhoL_values_copassing(vb)=sum(DE(POP_VPLL));
%     POP_VPLL=find(~(alphas_ejected).*ALL_TRAPPED_POP(1:length(alphas_ejected))...
%         .*(alphas_rhoL(1:length(alphas_ejected))>dr_bins(vb)).*(alphas_rhoL(1:length(alphas_ejected))<=dr_bins(vb+1)));
%     DE_rhoL_values_trapped(vb)=sum(DE(POP_VPLL));
%     POP_VPLL=find(~(alphas_ejected).*STAGNATION_POP(1:length(alphas_ejected))...
%         .*(alphas_rhoL(1:length(alphas_ejected))>dr_bins(vb)).*(alphas_rhoL(1:length(alphas_ejected))<=dr_bins(vb+1)));
%     DE_rhoL_values_stagnation(vb)=sum(DE(POP_VPLL));
end

DE_rhoL_values_stagnation=DE_rhoL_values_stagnation/max(DE_rhoL_values_trapped);
DE_rhoL_values_counterpassing=DE_rhoL_values_counterpassing/max(DE_rhoL_values_trapped);
DE_rhoL_values_copassing=DE_rhoL_values_copassing/max(DE_rhoL_values_trapped);
DE_rhoL_values_trapped=DE_rhoL_values_trapped/max(DE_rhoL_values_trapped);







%%
figure(16)
DRBINSIZE=0.006
dr_bins=(0:DRBINSIZE:0.08);
dr_values=dr_bins(1:end-1)+0.5*DRBINSIZE;

set(gca,'fontsize',20)

hold on

hold on
grid on
plot(dr_values,DE_rhoL_values_trapped,'k','linewidth',2)
plot(dr_values,DE_rhoL_values_counterpassing,'--','color',[0 0 0.7],'linewidth',2)
plot(dr_values,DE_rhoL_values_copassing,'b','linewidth',2)
% plot(dr_values,DE_dr_values_stagnation,'r','linewidth',2)

xl=xlabel('$\rho_L (\rm{m})$')
set(xl,'Interpreter','latex')
yl=ylabel('$\Delta \mathcal{E} (a.u.)$')
set(yl,'Interpreter','latex')
% ylim([-4.0 1.05])

legend('trapped','counter-passing','co-passing')
% legend('trapped','counter-passing','co-passing','stagnation')

xlim([0 0.066])