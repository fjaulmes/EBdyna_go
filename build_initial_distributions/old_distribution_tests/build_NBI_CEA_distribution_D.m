clear all;
initialize_folder_names;
filename=strcat(DATA_FOLDER,'q_profile.mat');
load(filename);
filename=strcat(DATA_FOLDER,'physics_constants.mat');
load(filename);
filename=strcat(DATA_FOLDER,'pressure_profile.mat');
load(filename);
filename=strcat(DATA_FOLDER,'motions_map_dimensions.mat');
load(filename);
filename=strcat(DATA_FOLDER,'flux_geometry.mat');
load(filename);
filename=strcat(DATA_FOLDER,'XZsmall_fields_tokamak_pre_collapse.mat');
load(filename);

% from eqdsk file
rgeo=6.3397;
zgeo=0.3999;
rmaxis=6.6989;
zmaxis=0.5196;

% for consistency with our shift
% rmaxis=(R0+X_axis)
% zmaxis=0.4*zmaxis

import_NBI_ions_distrib;

PROCESS_NUMBER=15

NB_PROCESS=16

%pause
close all

% Main important values for description of the distribution
mHe=mD
ZHe=1

NBI_RESCALE=0.99
% DIM_RESCALE=0.99

Nalphas_simulated=length(iondistrib(:,1))

alphas_pos_x=zeros(Nalphas_simulated,1);
alphas_pos_z=zeros(Nalphas_simulated,1);
alphas_pos_phi=zeros(Nalphas_simulated,1);
alphas_Ekin=zeros(Nalphas_simulated,1);
alphas_mm=zeros(Nalphas_simulated,1);
alphas_vpll=zeros(Nalphas_simulated,1);


% INIT_DATA_POS=1

d_pos_x=iondistrib(:,1)-rgeo;
d_pos_z=-iondistrib(:,2)+zgeo;
d_Ekin=NBI_RESCALE*iondistrib(:,3);
d_pitch=iondistrib(:,4);
d_vtot=sqrt(2*(eV/mHe)*d_Ekin);
% d_vpll=sign(d_pitch).*sqrt(2*(eV/mHe)*d_Ekin.*(d_pitch.^2)./(1+d_pitch.^2));
d_vpll=d_pitch.*d_vtot;
d_weight=iondistrib(:,5);

alphas_pos_x=NBI_RESCALE*[d_pos_x ];
alphas_pos_z=NBI_RESCALE*[d_pos_z ];
alphas_Ekin=[d_Ekin];
alphas_vpll=[d_vpll];

Epll=0.5*(mHe/eV)*alphas_vpll.^2;
Eperp=max(alphas_Ekin-0.5*(mHe/eV)*alphas_vpll.^2,0);
alphas_vperp=sqrt(2*(eV/mHe)*Eperp);
alphas_vtot=sqrt(2*(eV/mHe)*alphas_Ekin);

particles_weight=2*mean(iondistrib(:,5))/NB_PROCESS



disp('number of particles generated');
Nalphas_simulated

%%
% display particles
figure(1);
set(gca,'FontSize',22);

hold on;
axis xy square


hold on; grid on
plot(alphas_pos_x+R0,alphas_pos_z,'b.');
contour(scale_X+R0,scale_Z,psi_XZsmall_map',psi_scale(2:22:end),'k')
contour(scale_X+R0,scale_Z,psi_XZsmall_map',psi_scale(psi_rank_q1),'r','linewidth',4)
contour(scale_X+R0,scale_Z,psi_XZsmall_map',[psi_scale(end) psi_scale(end)] ,'k','linewidth',4)

xlabel('R')
ylabel('Z')

pause(0.1);

alphas_psi_value=interp2(scale_X,scale_Z,psi_XZsmall_map',alphas_pos_x,alphas_pos_z,'*cubic');
alphas_psi=interp1(psi_scale,1:257,alphas_psi_value);


% filling up the uniformly distributes toroidal position values
for(n=1:Nalphas_simulated)
    alphas_pos_phi(n)=my_rand(1)*2*pi;
end

alphas_mm=(alphas_Ekin-0.5*(mHe*alphas_vpll.^2)/eV)./interp2(scale_X,scale_Z,Btot_XZ_map',alphas_pos_x,alphas_pos_z,'*linear');

NEGATIVE_MM=find(alphas_mm<0);

disp('Beware of negative mm!!!!!!')
disp('length(NEGATIVE_MM)');
disp(length(NEGATIVE_MM));
alphas_mm=max(alphas_mm,0);

EXCESS_PSI=find((alphas_psi>=Nradial-1).*isnan(alphas_psi));

disp('Beware of particles outside of simulation!!!!!!')
disp('length(EXCESS_PSI)');
disp(length(EXCESS_PSI));
alphas_mm=max(alphas_mm,0);

INNER_PSI=find((alphas_psi<=Nradial-1).*(~isnan(alphas_psi)));


alphas_pos_x_global=alphas_pos_x(INNER_PSI);
alphas_pos_z_global=alphas_pos_z(INNER_PSI);
alphas_pos_phi_global=alphas_pos_phi(INNER_PSI);
alphas_Ekin_global=alphas_Ekin(INNER_PSI);
alphas_mm_global=alphas_mm(INNER_PSI);
alphas_vpll_global=alphas_vpll(INNER_PSI);
Nalphas_simulated=length(alphas_vpll_global)


disp('---------------------------------------------')
disp('splitting the data in two subsets for speedup')
disp('---------------------------------------------')
alphas_psi_value=interp2(scale_X,scale_Z,psi_XZsmall_map',alphas_pos_x_global,alphas_pos_z_global,'*cubic');
alphas_psi=interp1(psi_scale,1:257,alphas_psi_value);
INNER_PSI=find((alphas_psi<=Nradial-1).*(~isnan(alphas_psi)));
half_Nalphas=round(0.5*Nalphas_simulated)

alphas_pos_x=alphas_pos_x_global(INNER_PSI(1:half_Nalphas));
alphas_pos_z=alphas_pos_z_global(INNER_PSI(1:half_Nalphas));
alphas_pos_phi=alphas_pos_phi_global(INNER_PSI(1:half_Nalphas));
alphas_Ekin=alphas_Ekin_global(INNER_PSI(1:half_Nalphas));
alphas_mm=alphas_mm_global(INNER_PSI(1:half_Nalphas));
alphas_vpll=alphas_vpll_global(INNER_PSI(1:half_Nalphas));
Nalphas_simulated=length(alphas_vpll)
hold on; grid on
plot(alphas_pos_x+R0,alphas_pos_z,'r.');



FILENAME=strcat('initial_NBI_1MEV_D_distribution',num2str(PROCESS_NUMBER),'.mat')
save (FILENAME,'alphas_pos_x','alphas_pos_z','alphas_pos_phi','alphas_Ekin','alphas_mm','alphas_vpll','Nalphas_simulated','particles_weight','mHe','ZHe');


alphas_pos_x=alphas_pos_x_global(INNER_PSI(half_Nalphas+1:end));
alphas_pos_z=alphas_pos_z_global(INNER_PSI(half_Nalphas+1:end));
alphas_pos_phi=alphas_pos_phi_global(INNER_PSI(half_Nalphas+1:end));
alphas_Ekin=alphas_Ekin_global(INNER_PSI(half_Nalphas+1:end));
alphas_mm=alphas_mm_global(INNER_PSI(half_Nalphas+1:end));
alphas_vpll=alphas_vpll_global(INNER_PSI(half_Nalphas+1:end));
Nalphas_simulated=length(alphas_vpll)
hold on; grid on
plot(alphas_pos_x+R0,alphas_pos_z,'g.');


FILENAME=strcat('initial_NBI_1MEV_D_distribution',num2str(PROCESS_NUMBER+1),'.mat')
save (FILENAME,'alphas_pos_x','alphas_pos_z','alphas_pos_phi','alphas_Ekin','alphas_mm','alphas_vpll','Nalphas_simulated','particles_weight','mHe','ZHe');


figure(2)
hist(alphas_Ekin,20)
xlabel('Ekin (eV)')
