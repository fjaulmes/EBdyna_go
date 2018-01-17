% load('initialG_alphas_vA_all_pre_collapse.mat')
% load('alphas_vA_all_collapse_Glisa_fc2h2_G100214.mat')

Nalphas_simulated=length(find(~alphas_ejected))
volume_vessel=sum(sum(volume_tor_diff));
avg_density_alphas=0.005*Ne0
Frescale=(1/volume_vessel)/Nalphas_simulated

Btot_PR_map=sqrt(Bpol_PR_map.^2+Btor_PR_map.^2);
Btot_radial_profile=mean(Btot_PR_map(1:NP-1,:),1);

n=8
m=9

qTAE=(m+0.5)/n
kpll=1/(2*R0*qTAE)


PSI_BIN_SIZE=8
PSI_BIN_SIZE_HALF=round(0.5*PSI_BIN_SIZE);

psi_pos_bins=PSI_BIN_SIZE_HALF+1:PSI_BIN_SIZE:Nradial-PSI_BIN_SIZE_HALF
psi_bins=psi_scale(PSI_BIN_SIZE_HALF+1:PSI_BIN_SIZE:end-PSI_BIN_SIZE_HALF);
NB_PSI_BINS=length(psi_bins)

psi_bins_inf=psi_scale(1:PSI_BIN_SIZE:end-PSI_BIN_SIZE+1);
psi_bins_sup=psi_scale(PSI_BIN_SIZE+1:PSI_BIN_SIZE:end);
psi_bins_lim=(1:PSI_BIN_SIZE:257);

q_values=interp1(psi_scale,q_final_profile_diff,psi_bins);
% Ne_values=interp1(psi_scale,Ne_profile,psi_bins);
% B_values=interp1(psi_scale,Ne_profile,psi_bins);

% alphas_vperp=sqrt(2*Bavg*eV*alphas_mm/mDT);
alphas_vperp=sqrt(2*alphas_Ekin*eV/mHe-alphas_vpll.^2);

VBIN_SIZE=10*1e5;
V_HALFBIN_SIZE=0.5*VBIN_SIZE;
vperp_bins_lim=(0:VBIN_SIZE:16*VBIN_SIZE);
vperp_bins=(0:VBIN_SIZE:15*VBIN_SIZE)+V_HALFBIN_SIZE;
VPERP0_BIN_POS=1;
vpll_bins=(-15*VBIN_SIZE:VBIN_SIZE:15*VBIN_SIZE)+V_HALFBIN_SIZE;
VPLL0_BIN_POS=round(interp1(vpll_bins,1:length(vpll_bins),0));

% Ekin_bins_lim1=0.5*(mHe/eV)*(vA^2+vperp_bins_lim.^2);
% Ekin_bins_lim2=0.5*(mHe/eV)*((vA/3)^2+vperp_bins_lim.^2);
% Ekin_bins1=0.5*(mHe/eV)*((vA/3)^2+vperp_bins_lim.^2);
% Ekin_bins2=0.5*(mHe/eV)*((vA/3)^2+vperp_bins_lim.^2);
Ekin_bins_lim=0.5*(mHe/eV)*(vperp_bins_lim.^2);
Ekin_bins=0.5*(mHe/eV)*(vperp_bins.^2);

psi_pos_TAE=interp1(q_values,(1:length(psi_pos_bins)),qTAE)
Ne_value=interp1(1:Nradial,Ne_profile,psi_pos_TAE)
B_value=interp1(1:Nradial,Btot_radial_profile,psi_pos_TAE)
vA=B_value/sqrt(mu0*Ne_value*mDT)
wTAE=kpll*vA
DV=0.1*vA;


Fpsi_pop_inf=(alphas_psi>=(psi_pos_TAE-1.5)*PSI_BIN_SIZE).*(alphas_psi<=(psi_pos_TAE-0.5)*PSI_BIN_SIZE);
Fpsi_pop=(alphas_psi>=(psi_pos_TAE-0.5)*PSI_BIN_SIZE).*(alphas_psi<=(psi_pos_TAE+0.5)*PSI_BIN_SIZE);
Fpsi_pop_sup=(alphas_psi>=(psi_pos_TAE+0.5)*PSI_BIN_SIZE).*(alphas_psi<=(psi_pos_TAE+1.5)*PSI_BIN_SIZE);
vpll_pop1=(alphas_vpll>=vA/3-DV).*(alphas_vpll<=vA/3+DV);
vpll_pop2=(alphas_vpll>=vA-DV).*(alphas_vpll<=vA+DV);

RES_POP=find(Fpsi_pop.*vpll_pop1);
Fvperp1=histc(alphas_vperp(RES_POP),vperp_bins_lim);
FEkin1=histc(alphas_Ekin(RES_POP),Ekin_bins_lim);
dFvperp1=gradient(Fvperp1(1:end-1),vperp_bins);
dFEkin1=gradient(FEkin1(1:end-1),Ekin_bins);

RES_POP=find(Fpsi_pop.*vpll_pop2);
Fvperp2=histc(alphas_vperp(RES_POP),vperp_bins_lim);
FEkin2=histc(alphas_Ekin(RES_POP),Ekin_bins_lim);
dFvperp2=gradient(Fvperp2(1:end-1),vperp_bins);
dFEkin2=gradient(FEkin2(1:end-1),Ekin_bins);

dFpsi1=zeros(length(vperp_bins),length(psi_bins));
dFpsi2=zeros(length(vperp_bins),length(psi_bins));

for (vperp_pos=VPERP0_BIN_POS:length(vperp_bins))
    
    vperp_pop=(alphas_vperp>=(vperp_pos-0.5)*VBIN_SIZE).*(alphas_vperp<=(vperp_pos+0.5)*VBIN_SIZE);
    
    figure(1)
    grid on;hold on;
    RES_POP=find(vpll_pop1.*vperp_pop);
    Fpsi1=histc(alphas_psi(RES_POP),psi_bins_lim);
    dFpsi1(vperp_pos,:)=gradient(Fpsi1(1:end-1),psi_bins);
    RES_POP=find(vpll_pop2.*vperp_pop);
    Fpsi2=histc(alphas_psi(RES_POP),psi_bins_lim);
    dFpsi2(vperp_pos,:)=gradient(Fpsi2(1:end-1),psi_bins);
    plot(q_values,Fpsi1(1:end-1));plot(q_values,Fpsi2(1:end-1),'r')
    % plot(q_values,dFpsi1);hold on;plot(q_values,dFpsi2,'r')
    
end

integ1=0;
integ1_Ekin=0;
integ2=0;
integ2_Ekin=0;

for (vperp_pos=VPERP0_BIN_POS:length(vperp_bins))
%     vpll=vA;
%     vperp=vperp_pos*VBIN_SIZE;
%     integ1_Ekin=integ1_Ekin+(vpll^2+0.5*vperp^2)*((wTAE)*(vpll^2+0.5*vperp^2)*vperp*dFEkin1(vperp_pos))*VBIN_SIZE;
%     integ1=integ1+(vpll^2+0.5*vperp^2)*((wTAE)*(vpll^2+0.5*vperp^2)*vperp*dFEkin1(vperp_pos)-(n/ZHe)*(vpll^2+0.5*vperp^2)*vperp*interp1(1:length(psi_bins),dFpsi1(vperp_pos,:),psi_pos_TAE))*VBIN_SIZE;
%     vpll=vA/3;
%     integ2=integ2+(vpll^2+0.5*vperp^2)*((wTAE)*(vpll^2+0.5*vperp^2)*vperp*dFEkin2(vperp_pos)-(n/ZHe)*(vpll^2+0.5*vperp^2)*vperp*interp1(1:length(psi_bins),dFpsi2(vperp_pos,:),psi_pos_TAE))*VBIN_SIZE;
%     integ2_Ekin=integ2_Ekin+(vpll^2+0.5*vperp^2)*((wTAE)*(vpll^2+0.5*vperp^2)*vperp*dFEkin2(vperp_pos))*VBIN_SIZE;
    vpll=vA;
    vperp=vperp_pos*VBIN_SIZE;
    integ1_Ekin=integ1_Ekin+(vpll^2+0.5*vperp^2)*((wTAE)*(eV/mHe)*(vpll^2+0.5*vperp^2)*dFvperp1(vperp_pos))*VBIN_SIZE;
    integ1=integ1+(vpll^2+0.5*vperp^2)*((wTAE)*(eV/mHe)*(vpll^2+0.5*vperp^2)*dFvperp1(vperp_pos)-(n/ZHe)*(vpll^2+0.5*vperp^2)*vperp*interp1(1:length(psi_bins),dFpsi1(vperp_pos,:),psi_pos_TAE))*VBIN_SIZE;
    vpll=vA/3;
    integ2=integ2+(vpll^2+0.5*vperp^2)*((wTAE)*(eV/mHe)*(vpll^2+0.5*vperp^2)*dFvperp2(vperp_pos)-(n/ZHe)*(vpll^2+0.5*vperp^2)*vperp*interp1(1:length(psi_bins),dFpsi2(vperp_pos,:),psi_pos_TAE))*VBIN_SIZE;
    integ2_Ekin=integ2_Ekin+(vpll^2+0.5*vperp^2)*((wTAE)*(eV/mHe)*(vpll^2+0.5*vperp^2)*dFvperp2(vperp_pos))*VBIN_SIZE;
end

gamma1=Frescale*(1/eV)*(1/Bavg)*kpll*vA*2*(pi^2)*mu0*(mHe^2)*R0*(qTAE^3)*integ1
gamma2=Frescale*(1/eV)*(1/Bavg)*kpll*vA*2*(pi^2)*mu0*(mHe^2)*R0*(qTAE^3)*integ2
gamma_TAE=gamma1+gamma2

figure(2)
grid on
set(gca,'fontsize',20)
hold on
mHist = hist2d ([alphas_vpll, alphas_vperp], vpll_bins, vperp_bins);
contour(vpll_bins(1:end-1)+V_HALFBIN_SIZE,vperp_bins(1:end-1)+V_HALFBIN_SIZE,mHist',12)
plot([vA-DV vA-DV],[0 vperp_bins(end)],'r','linewidth',2)
plot([vA+DV vA+DV],[0 vperp_bins(end)],'r','linewidth',2)
plot([vA/3-DV vA/3-DV],[0 vperp_bins(end)],'r','linewidth',2)
plot([vA/3+DV vA/3+DV],[0 vperp_bins(end)],'r','linewidth',2)
xlabel('v_{||} (m/s)')
ylabel('v_{perp} (m/s)')
% plot2Dhist(mHist', vpll_bins,vperp_bins,  'vpll', 'vperp', 'Nb part'); 


RES_POP=find(Fpsi_pop);
disp('safety factor position')
psi_rank=interp1(1:length(psi_pos_bins),psi_pos_bins,psi_pos_TAE)

close all
figure(3)
grid on
set(gca,'fontsize',20)
title(strcat('q=',num2str(qTAE)))
hold on
mHist = hist2d ([alphas_vpll(RES_POP), alphas_vperp(RES_POP)], vpll_bins, vperp_bins);
contour(vpll_bins(1:end-1)+V_HALFBIN_SIZE,vperp_bins(1:end-1)+V_HALFBIN_SIZE,mHist',50:50:450)
plot([vA-DV vA-DV],[0 vperp_bins(end)],'r','linewidth',2)
plot([vA+DV vA+DV],[0 vperp_bins(end)],'r','linewidth',2)
plot([vA/3-DV vA/3-DV],[0 vperp_bins(end)],'r','linewidth',2)
plot([vA/3+DV vA/3+DV],[0 vperp_bins(end)],'r','linewidth',2)
xlim([-vperp_bins(end-4) vperp_bins(end-4)])
ylim([vpll_bins(17) vpll_bins(end-4)])
xlabel('v_{||} (m/s)')
ylabel('v_{perp} (m/s)')

disp('Number of alphas in radial slice=')
length(RES_POP)
% figure(2)
% RES_POP=find(Fpsi_pop.*vpll_pop1);
% Fvperp1=histc(alphas_vperp(RES_POP),vperp_bins);
% disp('length(RES_POP1)=')
% length(RES_POP)
% 
% RES_POP=find(Fpsi_pop.*vpll_pop2);
% disp('length(RES_POP2)=')
% length(RES_POP)
% Fvperp2=histc(alphas_vperp(RES_POP),vperp_bins);
% plot(Fvperp1);hold on;plot(Fvperp2,'r')

