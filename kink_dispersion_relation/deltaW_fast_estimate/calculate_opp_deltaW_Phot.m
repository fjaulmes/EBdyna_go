
reset_data_analysis_environment;

load('NBIopp_Phot_data.mat')
scale_X_P=scale_X;
scale_Z_P=scale_Z;

initialize_xi_map_calculation_context;
rescaling_to_XZsmall_maps;
% [XX_P ZZ_P]=meshgrid(scale_X_P,scale_Z_P);
NBI_Phot_XZ_map_small=interp2(scale_X_P,scale_Z_P,NBI_Phot_XZ_map',XX_small,ZZ_small);
NBI_Phot_XZ_map_small=NBI_Phot_XZ_map_small';


deltaW_est_phi=zeros(NB_PHI-1,1);
deltaW_XZ_map_phi=zeros(NB_PHI-1,sizeX,sizeZ);

Bphi_PR_map=Btor_PR_map(:,1:size_r);

load('../../E_maps/E0141.mat');
load('../../B_maps/B0141.mat');

frame_rank=15
% phi_index=45
calculate_gradPhot;

psi_star_dot_omega_map(:,:)=psi_star_dot_evol(frame_rank,:,:);

[scXX scZZ]=meshgrid(scale_X,scale_Z);

xi0_core=ksi0_evol_lin(frame_rank)
core_radial_pos=round(interp1(radial_r_value_flux,1:Nradial,xi0_core));
volume_XZ_map=Rpos_XZsmall_map.*2*pi*DX*DX;
tot_vol=sum(sum(volume_XZ_map(psi_XZ_map_mask~=0)));

% should we use the finesse output pressure ??

% P_initial_profile=PTOT_profile_interp_ini;

% we'd rahter not, since it hsa to match against the B field we have ?????


for index=1:NB_THETA
    %         phi_rank=round(2*phi_index-1);
    P_PR_map(index,:)=P_initial_profile;
end
P_data=reshape(P_PR_map(:,1:radial_q1_size),NP*radial_q1_size,1);
P_data=double(P_data);
% P_XZ_map=griddata(finesse_data_X,finesse_data_Z,P_data,XX_small,ZZ_small);
% P_XZ_map(isnan(P_XZ_map))=0;
% P_XZ_map=P_XZ_map';

P_XZ_map=gridfit(finesse_data_X,finesse_data_Z,P_data,scale_X,scale_Z,'smoothness',0.5);
P_XZ_map=P_XZ_map';
   
% calculate_gradP_XZ;

Btot_ini_XZ_map=sqrt(BpolX_initial_XZsmall_map.^2+BpolZ_initial_XZsmall_map.^2+Bphi_XZsmall_map.^2);
% calculate curvature
Etot_XZ_map=P_XZ_map+0.5*Btot_ini_XZ_map.^2/mu0;
calculate_gradE_XZ;
calculate_gradB_XZ;

BxgradB_X=-Bphi_XZsmall_map.*gradB_Z;
BxgradB_Z=Bphi_XZsmall_map.*gradB_X;
BxgradB_phi=BpolX_initial_XZsmall_map.*gradB_Z-BpolZ_initial_XZsmall_map.*gradB_X;

% small toroidal component
kappa_phi_XZ_map=BxgradB_X.*BpolZ_initial_XZsmall_map-BxgradB_Z.*BpolX_initial_XZsmall_map;
kappa_phi_XZ_map=kappa_phi_XZ_map./Btot_ini_XZ_map.^3;


kappa_X_XZ_map=mu0*(gradE_X./Btot_ini_XZ_map.^2);
kappa_Z_XZ_map=mu0*(gradE_Z./Btot_ini_XZ_map.^2);

kappa_X_XZ_map=kappa_X_XZ_map.*psi_XZ_map_mask;
kappa_Z_XZ_map=kappa_Z_XZ_map.*psi_XZ_map_mask;

save quiver_data.mat scXX scZZ kappa_X_XZ_map kappa_Z_XZ_map kappa_phi_XZ_map gradP_X gradP_Z

%%
for phi_index=1:NB_PHI-1
    disp('----------------------------------------------')
    phi_index
    
    psi_star_dot_PR_map_phi=zeros(NB_THETA,NB_THETA,size_r);
    for index=1:NB_THETA
        %         phi_rank=round(2*phi_index-1);
        psi_star_dot_PR_map_phi(index,:,:)=rotate_map_phi(psi_star_dot_omega_map,index)';
    end
    
    Btot_PR_map=squeeze(Btot_map_phi(phi_index,:,:));
    BpolX_PR_map=squeeze(bX_map_phi(phi_index,:,:)).*Btot_PR_map;
    BpolZ_PR_map=squeeze(bZ_map_phi(phi_index,:,:)).*Btot_PR_map;
    
    EX_PR_map=squeeze(Efield_X_map_phi(phi_index,:,:));
    EZ_PR_map=squeeze(Efield_Z_map_phi(phi_index,:,:));
    grad_Phi_tor_PR_map=squeeze(grad_Phi_tor_map_phi(phi_index,:,:));
    psi_star_dot_PR_map=squeeze(psi_star_dot_PR_map_phi(phi_index,:,:));
    
    Ephi_PR_map=psi_star_dot_PR_map./Rpos_PR_map(:,1:size_r)-grad_Phi_tor_PR_map;
    
    Ephi_PR_map_recalc=BpolX_PR_map.*EX_PR_map+BpolZ_PR_map.*EZ_PR_map;
    Ephi_PR_map_recalc=-Ephi_PR_map_recalc./Bphi_PR_map;
    Ephi_PR_map=0.5*(Ephi_PR_map+Ephi_PR_map_recalc);
%     Bphi_PR_map_recalc=BpolX_PR_map.*EX_PR_map+BpolZ_PR_map.*EZ_PR_map;
%     Bphi_PR_map_recalc=-Bphi_PR_map_recalc./Ephi_PR_map;
    
    vE_X_PR_map=EZ_PR_map.*Bphi_PR_map-Ephi_PR_map.*BpolZ_PR_map;
    vE_Z_PR_map=Ephi_PR_map.*BpolX_PR_map-EX_PR_map.*Bphi_PR_map;
    vE_phi_PR_map=EX_PR_map.*BpolZ_PR_map-EZ_PR_map.*BpolX_PR_map;
    
    vE_X_PR_map=vE_X_PR_map./Btot_PR_map.^2;
    vE_Z_PR_map=vE_Z_PR_map./Btot_PR_map.^2;
    vE_phi_PR_map=vE_phi_PR_map./Btot_PR_map.^2;
    
    vE_tot_PR_map=sqrt(vE_X_PR_map.^2+vE_Z_PR_map.^2+vE_phi_PR_map.^2);
    
%     E_data=reshape(BpolX_PR_map(:,1:radial_q1_size),NP*radial_q1_size,1);
%     EX_XZ_map=gridfit(finesse_data_X,finesse_data_Z,E_data,scale_X,scale_Z,'smoothness',0.5);
%     EX_XZ_map=EX_XZ_map';
%     
%     
%     % now get XZ map for B
%     B_data=reshape(BpolX_PR_map(:,1:radial_q1_size),NP*radial_q1_size,1);
%     BX_XZ_map=gridfit(finesse_data_X,finesse_data_Z,B_data,scale_X,scale_Z,'smoothness',0.5);
%     BX_XZ_map=BX_XZ_map';
%     
%     B_data=reshape(BpolZ_PR_map(:,1:radial_q1_size),NP*radial_q1_size,1);
%     BZ_XZ_map=gridfit(finesse_data_X,finesse_data_Z,B_data,scale_X,scale_Z,'smoothness',0.5);
%     BZ_XZ_map=BZ_XZ_map';
%     
%     B_data=reshape(Btot_PR_map(:,1:radial_q1_size),NP*radial_q1_size,1);
%     Btot_XZ_map=gridfit(finesse_data_X,finesse_data_Z,B_data,scale_X,scale_Z,'smoothness',0.5);
%     Btot_XZ_map=Btot_XZ_map';
    
    
    % calculate curvature
%     Etot_XZ_map=P_XZ_map+0.5*Btot_XZ_map.^2/mu0;
%     calculate_gradE_phi_rank;
%     calculate_gradB_phi_rank;
%     
%     BxgradB_X=-Bphi_XZsmall_map.*gradB_Z;
%     BxgradB_Z=Bphi_XZsmall_map.*gradB_X;
%     BxgradB_phi=BX_XZ_map.*gradB_Z-BZ_XZ_map.*gradB_X;
%     
%     kappa_X_XZ_map=BxgradB_Z.*Bphi_XZsmall_map-BxgradB_phi.*BZ_XZ_map;
%     kappa_Z_XZ_map=BxgradB_phi.*BX_XZ_map-BxgradB_X.*Bphi_XZsmall_map;
%     kappa_phi_XZ_map=BxgradB_X.*BZ_XZ_map-BxgradB_Z.*BX_XZ_map;
    
    % small toroidal component
%     kappa_phi_XZ_map=BxgradB_X.*BZ_XZ_map-BxgradB_Z.*BX_XZ_map;
%     kappa_phi_XZ_map=kappa_phi_XZ_map./Btot_XZ_map.^3;
% 
%     kappa_X_XZ_map=mu0*(gradE_X./Btot_XZ_map.^2);
%     kappa_Z_XZ_map=mu0*(gradE_Z./Btot_XZ_map.^2);
      
%     kappa_X_XZ_map=kappa_X_XZ_map./Btot_XZ_map.^3;
%     kappa_Z_XZ_map=kappa_Z_XZ_map./Btot_XZ_map.^3;
%     kappa_phi_XZ_map=kappa_phi_XZ_map./Btot_XZ_map.^3;
%     
%     kappa_X_XZ_map=mu0*gradP_X./Btot_XZ_map.^2+kappa_X_XZ_map;
%     kappa_Z_XZ_map=mu0*gradP_Z./Btot_XZ_map.^2+kappa_Z_XZ_map;
    
%     kappa_X_XZ_map=kappa_X_XZ_map.*psi_XZ_map_mask;
%     kappa_Z_XZ_map=kappa_Z_XZ_map.*psi_XZ_map_mask;
    

    disp('----------------------------------------------')
    %
    
    % now get XZ map for vE
    vE_data=reshape(vE_X_PR_map(:,1:radial_q1_size),NP*radial_q1_size,1);
    vE_X_XZ_map=gridfit(finesse_data_X,finesse_data_Z,vE_data,scale_X,scale_Z,'smoothness',0.5);
    vE_X_XZ_map=vE_X_XZ_map';
    
    vE_data=reshape(vE_Z_PR_map(:,1:radial_q1_size),NP*radial_q1_size,1);
    vE_Z_XZ_map=gridfit(finesse_data_X,finesse_data_Z,vE_data,scale_X,scale_Z,'smoothness',0.5);
    vE_Z_XZ_map=vE_Z_XZ_map';
    
    vE_data=reshape(vE_phi_PR_map(:,1:radial_q1_size),NP*radial_q1_size,1);
    vE_phi_XZ_map=gridfit(finesse_data_X,finesse_data_Z,vE_data,scale_X,scale_Z,'smoothness',0.5);
    vE_phi_XZ_map=vE_phi_XZ_map';
    
    vE_norm_XZ_map=sqrt(vE_X_XZ_map.^2+vE_Z_XZ_map.^2+vE_phi_XZ_map.^2);
    vE_norm_XZ_map=vE_norm_XZ_map.*psi_XZ_map_mask;
    vE_norm_vol_avg=sum(sum(volume_XZ_map(psi_XZ_map_mask~=0).*vE_norm_XZ_map(psi_XZ_map_mask~=0)))/tot_vol
%     vE_norm_axis=vE_tot_PR_map(257-2*(phi_index-1),0.5*(core_radial_pos+psi_rank_q1))
    vE_norm_axis1=max(vE_tot_PR_map(2*(phi_index-1)+1,1:round(0.75*psi_rank_q1)))
    vE_norm_axis2=max(vE_tot_PR_map(round(0.5*NP)-(phi_index-1),1:round(0.75*psi_rank_q1)))
    vE_norm_axis=max(vE_norm_axis1,vE_norm_axis2);
    vE_norm_avg=mean(vE_norm_XZ_map(vE_norm_XZ_map~=0));
    
    
    xi_X_norm_XZ_map= vE_X_XZ_map./vE_norm_axis;
    xi_Z_norm_XZ_map= vE_Z_XZ_map./vE_norm_axis;
    xi_phi_norm_XZ_map= vE_phi_XZ_map./vE_norm_axis;
    
    %
    xi_X_norm_XZ_map=xi_X_norm_XZ_map.*psi_XZ_map_mask;
    xi_Z_norm_XZ_map=xi_Z_norm_XZ_map.*psi_XZ_map_mask;
    xi_phi_norm_XZ_map=xi_phi_norm_XZ_map.*psi_XZ_map_mask;
    
    xi_norm_XZ_map=sqrt(xi_X_norm_XZ_map.^2+xi_Z_norm_XZ_map.^2+xi_phi_norm_XZ_map.^2);
    
    
%     imagesc(scale_X,scale_Z,xi_X_norm_XZ_map');
    
    
    %
    
    % now finally evaluate deltaW from Phot !
    
    deltaW_XZ_map=xi_X_norm_XZ_map.*gradP_X+xi_Z_norm_XZ_map.*gradP_Z;
    curvature_XZ_map=xi_X_norm_XZ_map.*kappa_X_XZ_map+xi_Z_norm_XZ_map.*kappa_Z_XZ_map+xi_phi_norm_XZ_map.*kappa_phi_XZ_map;
    
    deltaW_XZ_map=deltaW_XZ_map.*curvature_XZ_map;
    
    deltaW_XZ_map_phi(phi_index,:,:)=deltaW_XZ_map;
    
    deltaW_est=0;
    
    for x=1:sizeX
        deltaW_est=deltaW_est+DPHI*DX*DX*sum((scale_X(x)+R0).*squeeze(deltaW_XZ_map(x,1:end)));
    end
    deltaW_est=-deltaW_est;
    
    deltaW_est_hat=(mu0*(R0^3)/(6*pi^2*Bphi0^2*r_value_q1_mean^4))*deltaW_est
    %%
    
    deltaW_est_phi(phi_index)=deltaW_est;
    disp('----------------------------------------------')

end


%%
deltaW_fast_hat=(mu0*(R0^3)/(6*pi^2*Bphi0^2*r_value_q1_mean^4))*sum(deltaW_est_phi)


save('NBIopp_Phot_data.mat','-append','deltaW_fast_hat');
disp('----------------------------------------------')
save('deltaW_NBIopp_data.mat','deltaW_XZ_map_phi','deltaW_est_phi','deltaW_fast_hat');

%%
VECSP=32;

figure(1)
set(gca,'fontsize',20)

hold on
contour(scale_X,scale_Z,psi_XZsmall_map',12,'k.')
contour(scale_X,scale_Z,(NBI_Phot_XZ_map_small)',6,'linewidth',3)
% colormap('summer')
% xi_norm_XZ_map_plot=abs(xi_norm_XZ_map.^0.5);
% xi_X_norm_XZ_map_plot=xi_X_norm_XZ_map./xi_norm_XZ_map_plot;
% xi_Z_norm_XZ_map_plot=xi_Z_norm_XZ_map./xi_norm_XZ_map_plot;

hq=quiver(scXX(1:VECSP:end,1:VECSP:end-VECSP),scZZ(1:VECSP:end,1:VECSP:end-VECSP),xi_X_norm_XZ_map(1:VECSP:end-VECSP,1:VECSP:end)',xi_Z_norm_XZ_map(1:VECSP:end-VECSP,1:VECSP:end)',2.9,'k','linewidth',1.2);
adjust_quiver_arrowhead_size(hq,1.5)
xlim([-0.15 0.26])
ylim([-0.25 0.25])

xlabel('X (m)');
ylabel('Z (m)');

%%
figure(1)
hold on
contour(scale_X,scale_Z,-deltaW_XZ_map',20)
% colormap('summer')
quiver(scXX(1:VECSP:end,1:VECSP:end-VECSP),scZZ(1:VECSP:end,1:VECSP:end-VECSP),gradP_X(1:VECSP:end-VECSP,1:VECSP:end)',gradP_Z(1:VECSP:end-VECSP,1:VECSP:end)','k')
xlim([-0.15 0.26])
ylim([-0.25 0.25])

figure(2)
hold on
contour(scale_X,scale_Z,-deltaW_XZ_map',20)
quiver(scXX(1:VECSP:end,1:VECSP:end-VECSP),scZZ(1:VECSP:end,1:VECSP:end-VECSP),xi_X_norm_XZ_map(1:VECSP:end-VECSP,1:VECSP:end)',xi_Z_norm_XZ_map(1:VECSP:end-VECSP,1:VECSP:end)',2,'k')
xlim([-0.15 0.26])
ylim([-0.25 0.25])

figure(3)
hold on
contour(scale_X,scale_Z,-deltaW_XZ_map',20)
quiver(scXX(1:VECSP:end,1:VECSP:end-VECSP),scZZ(1:VECSP:end,1:VECSP:end-VECSP),kappa_X_XZ_map(1:VECSP:end-VECSP,1:VECSP:end)',kappa_Z_XZ_map(1:VECSP:end-VECSP,1:VECSP:end)',0.5,'k')
xlim([-0.15 0.26])
ylim([-0.25 0.25])


