close all;
% load('map_dimensions.mat')
% load('geomtry_bins.mat')
load('speed_bins.mat')
% load('XZsmall_fields_tokamak.mat')
% load('physics_constants.mat')
% 
load('initial_alphas_distributions.mat');
% load('initial_alphas_pre_collapse');

disp('Total number of particles simulates');
disp(Nalphas_simulated);

X0=alphas_pos_x;
Z0=alphas_pos_x;
vpll0=alphas_vpll;
Ekin0=alphas_Ekin;
mm0=alphas_mm;
Bfield0=interp2(scale_X,scale_Z,Btot_XZ_map',alphas_pos_x,alphas_pos_z,'*linear');
psi_pos0=interp2(scale_X,scale_Z,psi_norm_XZsmall_map',alphas_pos_x,alphas_pos_z,'*linear');

clear LOW_ENERGY_POPULATION

%load('alphas_ejected_count.mat');

alphas_lambda=alphas_mm./alphas_Ekin;
ENERGY_POPULATION_0=find(alphas_Ekin>2*1e3);
%LOW_ENERGY_POPULATION=find(~alphas_ejected);
ENERGY_POPULATION=find(alphas_Ekin>2*1e3);
%ENERGY_POPULATION_EJECTED=find(alphas_ejected.*(alphas_Ekin>2000*1e3).*(alphas_lambda<=0.3));

%ejected_fraction=size(ENERGY_POPULATION_EJECTED,1)/size(ENERGY_POPULATION_0,1);
post_dist_rescale=size(ENERGY_POPULATION_0,1)/size(ENERGY_POPULATION,1);

%LOW_ENERGY_POPULATION=find(abs(alphas_mm./alphas_Ekin<=0.3));
%LOW_ENERGY_POPULATION=find(abs(alphas_mm(LOW_ENERGY_POPULATION)./alphas_Ekin(LOW_ENERGY_POPULATION))<=0.3);
%LOW_ENERGY_POPULATION=find(alphas_Ekin<45*1e3);
%LOW_ENERGY_POPULATION=find(alphas_Ekin(LOW_ENERGY_POPULATION)<95*1e3);
disp('ENERGY_POPULATION')
disp(size(ENERGY_POPULATION,1));
%LOW_ENERGY_POPULATION=find(alphas_Ekin(LOW_ENERGY_POPULATION)<120*1e3);

%load('initial_alphas_pre_collapse');
load('alphas_pre_record_10000.mat');
%load('initial_alphas_pre_collapse');

alphas_pos_psi=interp2(scale_X,scale_Z,psi_norm_XZsmall_map',alphas_pos_x,alphas_pos_z,'*linear');

figure(1)
hold on
grid on;

Nvpll0=histc(X0(ENERGY_POPULATION_0),(-20:20)*0.04);xlim([-0.7 0.7]);
Nvpll=post_dist_rescale*histc(alphas_pos_x(ENERGY_POPULATION),(-20:20)*0.04);xlim([-0.7 0.7]);

plot((-20:20)*0.04,Nvpll0,'b');
plot((-20:20)*0.04,Nvpll,'r');

% h = findobj(gca,'Type','patch');
% set(h(2),'FaceColor','b','EdgeColor','k');
% set(h(1),'FaceColor','r');

legend('initial','after mixing')
title('axial position (X)');


figure(2)
hold on;
grid on;

Nvpll0=histc(vpll0(ENERGY_POPULATION_0),vpll_range-0.5*vpll_bin_size);
Nvpll=post_dist_rescale*histc(alphas_vpll(ENERGY_POPULATION),vpll_range-0.5*vpll_bin_size);

% h = findobj(gca,'Type','patch');
% set(h(2),'FaceColor','b','EdgeColor','k');
% set(h(1),'FaceColor','r');

plot(vpll_range,Nvpll0,'b');
plot(vpll_range,Nvpll,'r');

legend('initial','after mixing');
title('v_{||}');


figure(3)
hold on;
grid on;

Npsi0=density_part_ratio*histc(psi_pos0(ENERGY_POPULATION_0),psi_bin_pos-0.5*radial_bin_size);
Npsi=post_dist_rescale*density_part_ratio*histc(alphas_pos_psi(ENERGY_POPULATION),psi_bin_pos-0.5*radial_bin_size);

Npsi0=Npsi0./sum(volume_bin(:,:),2);
Npsi=Npsi./sum(volume_bin(:,:),2);
% hist(psi_pos0,30);
% hist(alphas_pos_psi,30);
% 
% h = findobj(gca,'Type','patch');
% set(h(2),'FaceColor','b','EdgeColor','k');
% set(h(1),'FaceColor','r');

plot(psi_bin_pos,Npsi0,'b');
plot(psi_bin_pos,Npsi,'r');


legend('initial','after mixing');
title('flux surface position (\psi)');



[XX_psi PSI_xx]=meshgrid(scale_X(X_bin_pos),psi_bin_pos);

% alphas_vpll0_map=griddata(psi_pos0,X0,abs(vpll0),PSI_xx,XX_psi);
% alphas_vpll0_map(isnan(alphas_vpll0_map))=0;
%alphas_vpll0_map=alphas_vpll0_map.*volume_bin;

% alphas_vpll_map=griddata(alphas_pos_psi,alphas_pos_x,abs(alphas_vpll),PSI_xx,XX_psi);
% alphas_vpll_map(isnan(alphas_vpll_map))=0;
%alphas_vpll_map=alphas_vpll_map.*volume_bin;

% figure(4)
% grid on; 
% hold on;
% plot(scale_X(X_bin_pos),mean(alphas_vpll0_map,1),'b');
% plot(scale_X(X_bin_pos),mean(alphas_vpll_map,1),'r');
% legend('initial','after mixing')
% title('v_{||}');
% xlim([-0.6 0.6]);

% figure(3)
% imagesc(scale_X,scale_Z,alphas_vpll0_map);
% axis xy
% figure(4)
% imagesc(scale_X,scale_Z,alphas_vpll_map);
% axis xy


% alphas_Ekin0_map=griddata(psi_pos0,X0,Ekin0,PSI_xx,XX_psi);
% alphas_Ekin0_map(isnan(alphas_Ekin0_map))=0;
%alphas_Ekin0_map=alphas_Ekin0_map.*volume_bin;

% alphas_Ekin_map=griddata(alphas_pos_psi,alphas_pos_x,alphas_Ekin,PSI_xx,XX_psi);
% alphas_Ekin_map(isnan(alphas_Ekin_map))=0;
%alphas_Ekin_map=alphas_Ekin_map.*volume_bin;
% 
% figure(5)
% grid on; 
% hold on;
% plot(scale_X(X_bin_pos),mean(alphas_Ekin0_map,1),'b');
% plot(scale_X(X_bin_pos),mean(alphas_Ekin_map,1),'r');
% legend('initial','after mixing')
% title('Ekin');
% xlim([-0.6 0.6]);

% figure(6)
% imagesc(scale_X,scale_Z,alphas_Ekin0_map);
% axis xy
% figure(7)
% imagesc(scale_X,scale_Z,alphas_Ekin_map);
% axis xy

figure(4)
hold on;
grid on;

p_phi0=mHe*(X0(ENERGY_POPULATION_0)+R0).*abs(vpll0(ENERGY_POPULATION_0))-2*eV*(psi_pos0(ENERGY_POPULATION_0)-257)*psi_global/257;
p_phi=mHe*(alphas_pos_x(ENERGY_POPULATION)+R0).*abs(alphas_vpll(ENERGY_POPULATION))-2*eV*(alphas_pos_psi(ENERGY_POPULATION)-257)*psi_global/257;

Npsi0=density_part_ratio*histc(abs(p_phi0),(1:20)*max(abs(p_phi0))/20);
Npsi=post_dist_rescale*density_part_ratio*histc(abs(p_phi),(1:20)*max(abs(p_phi0))/20);

% hist(psi_pos0,30);
% hist(alphas_pos_psi,30);
% 
% h = findobj(gca,'Type','patch');
% set(h(2),'FaceColor','b','EdgeColor','k');
% set(h(1),'FaceColor','r');

plot(1:20,Npsi0,'b');
plot(1:20,Npsi,'r');


legend('initial','after mixing');
title('canonical parallel angular momentum');