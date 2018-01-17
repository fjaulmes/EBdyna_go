
% Deuterium particles

ZHe=1
mHe=mD


load('../data_tokamak/q_profile.mat')
load('../data_tokamak/flux_geometry.mat', 'X_PR_map')
load('../data_tokamak/flux_geometry.mat', 'Z_PR_map')
theta_values=(0:256)*2*pi/256;

pphi_recalc_outputG=(mHe/eV)*(Xpos_outputG(:,1:end)+R0).*(vphi_outputG(:,1:end))-ZHe*(psi_value_outputG(:,1:end));

PPHICORE=min(min(pphi_recalc_outputG))
PPHIMIN=max(max(pphi_recalc_outputG))
pphi_recalc_outputG=(pphi_recalc_outputG-PPHIMIN);
PPHICORE=min(min(pphi_recalc_outputG))+0.029;
PPHICORE=min(min(pphi_recalc_outputG));
pphi_recalc_outputG=pphi_recalc_outputG/PPHICORE;
PPHICORE=max(max(pphi_recalc_outputG));
pphi_recalc_outputG=pphi_recalc_outputG-(PPHICORE-1);

PHI_SMALL_ANGLE=4.0e-2

POLOIDAL_PLOT=0;

figure(1)
set(gca,'fontsize',20)

hold on
xlabel('\theta')
yl=ylabel('$-p_\varphi$');
set(yl,'Interpreter','latex')

% psi_qm5n2=interp1(q_initial_profile,psi_scale,2.5);
% psi_qm2n1=interp1(q_initial_profile,psi_scale,2);
% psi_qm3n1=interp1(q_initial_profile,psi_scale,3);
% psi_qm3n2=interp1(q_initial_profile,psi_scale,1.5);
% psi_qm1n1=interp1(q_initial_profile,psi_scale,1);
% psi_qm7n2=interp1(q_initial_profile,psi_scale,3.5);
% psi_qm4n1=interp1(q_initial_profile,psi_scale,4);
% 
% psi_qm4n3=interp1(q_initial_profile,psi_scale,4/3);
% psi_qm6n5=interp1(q_initial_profile,psi_scale,6/5);
% psi_qm7n6=interp1(q_initial_profile,psi_scale,7/6);
% psi_qm5n4=interp1(q_initial_profile,psi_scale,5/4);
% 
% 
% % plot([0 2*pi],[psi_qm5n4 psi_qm5n4],'g--','linewidth',3)
% plot([0 2*pi],[psi_qm2n1 psi_qm2n1],'r','linewidth',3)
% plot([0 2*pi],[psi_qm3n2 psi_qm3n2],'r--','linewidth',3)
% plot([0 2*pi],[psi_qm5n2 psi_qm5n2],'r--','linewidth',3)
% plot([0 2*pi],[psi_qm3n1 psi_qm3n1],'r','linewidth',3)
% plot([0 2*pi],[psi_qm1n1 psi_qm1n1],'r','linewidth',3)
% plot([0 2*pi],[psi_qm7n2 psi_qm7n2],'r--','linewidth',3)
% plot([0 2*pi],[psi_qm4n1 psi_qm4n1],'r','linewidth',3)

if POLOIDAL_PLOT==1
    figure(2)
    set(gca,'fontsize',20)
    
    hold on
    contour(scale_X+R0,scale_Z,psi_norm_XZsmall_map',32)
    xlabel('R')
    ylabel('Z')
end
% 
% INITIME=16000
% ENDTIME=size(phipos_outputG,2)-1000
INITIME=1
ENDTIME=size(phipos_outputG,2)-1
% 

for n=1:128
    % find the particles at phi = 0
    PARTPOP=(wrap2pi(phipos_outputG(n,INITIME:ENDTIME))<PHI_SMALL_ANGLE)+(abs(wrap2pi(phipos_outputG(n,INITIME:ENDTIME))-2*pi)<PHI_SMALL_ANGLE);
    PARTPOP(PARTPOP>1)=1;
    TORPOSLOOP=find(PARTPOP)+INITIME-1;
    figure(1);
%     plot(theta_outputG(n,TORPOSLOOP),psi_value_outputG(n,TORPOSLOOP),'.','color',[0 0 0.6])
%     plot(theta_outputG(n,TORPOSLOOP),pphi_recalc_outputG(n,TORPOSLOOP),'b.')
    plot(theta_outputG(n,TORPOSLOOP),pphi_recalc_outputG(n,TORPOSLOOP),'.','color',[0 0 0.6])
%     plot(theta_outputG(n,TORPOSLOOP),pphi_recalc_outputG(n,TORPOSLOOP),'r.')
    if POLOIDAL_PLOT==1
        figure(2)
        Rvalues=interp2(theta_values,psi_scale,X_PR_map',theta_outputG(n,TORPOSLOOP),psi_value_outputG(n,TORPOSLOOP))+R0;
        Zvalues=interp2(theta_values,psi_scale,Z_PR_map',theta_outputG(n,TORPOSLOOP),psi_value_outputG(n,TORPOSLOOP));
        plot(Rvalues,Zvalues,'k.')
    end
end

figure(1);



xlim([0 2*pi])



% 
% for n=1:128
%     % find the particles at phi = 0
%     PARTPOP=(wrap2pi(phipos_outputG(n,INITIME:ENDTIME))<PHI_SMALL_ANGLE)+(abs(wrap2pi(phipos_outputG(n,INITIME:ENDTIME))-2*pi)<PHI_SMALL_ANGLE);
%     PARTPOP(PARTPOP>1)=1;
%     TORPOSLOOP=find(PARTPOP);
%     figure(1);
%     plot(theta_outputG(n,TORPOSLOOP),psi_value_outputG(n,TORPOSLOOP),'.','color',[0 0 0.6])
%     if POLOIDAL_PLOT==1
%         figure(2)
%         Rvalues=interp2(theta_values,psi_scale,X_PR_map',theta_outputG(n,TORPOSLOOP),psi_value_outputG(n,TORPOSLOOP))+R0;
%         Zvalues=interp2(theta_values,psi_scale,Z_PR_map',theta_outputG(n,TORPOSLOOP),psi_value_outputG(n,TORPOSLOOP));
%         plot(Rvalues,Zvalues,'k.')
%     end
% end
