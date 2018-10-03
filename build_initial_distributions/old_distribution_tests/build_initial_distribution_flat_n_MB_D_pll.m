clear all;
initialize_folder_names;
build_radial_MBD_binning_TAE;
energy_D_range=energy_H_range;
f_D_energy_percentage=f_H_energy_percentage;
save f_D_energy_percentage_map.mat  radial_bins energy_D_range f_D_energy_percentage

mHe=mD
ZHe=1
% filename=strcat(DATA_FOLDER,'plasma_current.mat');
% load(filename, 'JPHI_XZ_map')

% rmaxis = 1.6715
% zmaxis = 0.0747


if exist('OMEGA_profile_radial_ini')
    rotation_profile=OMEGA_profile_radial_ini;
else
    rotation_profile=Ne_profile*0;
end
CORE_ROTATION=rotation_profile(1)

% counter field rotation (co-current)
% rotation_profile=-rotation_profile;

%%
filename=strcat(DATA_FOLDER,'tokamak_PR_map.mat');
load(filename);
filename=strcat(DATA_FOLDER,'flux_geometry.mat');
load(filename);
filename=strcat(DATA_FOLDER,'volume_flux_geometry.mat');
load(filename);
filename=strcat(DATA_FOLDER,'XZsmall_fields_tokamak_pre_collapse.mat');
load(filename);

Npart_simulated=5.0e4;
NPART_DENSITY_RATIO=1.69e-13

NB_PROCESS=16

%pause
%

psi_rank_max=210

Ni_avg=mean(fast_density_profile(1:psi_rank_max))

close all


%Normalize the alpha density
Nalpha_max=max(Nalpha_binned);
Nalpha_binned_norm=Nalpha_binned/Nalpha_max;


% this number need to be distributed
% according to the density of particles



% loop through all flux surface bins
delta_psi=50;
psi_bin_pos=1+radial_bin_half_size+radial_bin_size*(0:N_radial_bins-1)
% theta_bin_values=theta_bin_bounds(1:end-1)+0.5*(theta_bin_bounds(2)-theta_bin_bounds(1));

% volume of vertical sclices of toroidal dougnnuts shaped bins
% (yummy)



%
X_min=mid_Xaxis_large-0.5*size(scale_X,2)

Npart_rank_begin=1;


% display options
figure(1);
set(gca,'FontSize',22);

hold on;
axis xy square


%should be an even number
X_bin_size=10;
X_bin_half_std=round(0.5*(X_bin_size));

X_min=X1_Nradial(psi_rank_max+0.3*X_bin_size);
X_max=X2_Nradial(psi_rank_max+2*X_bin_size);

X_offset=-round(0.5*mod(X_max-X_min,X_bin_size))+1+X_bin_half_std;
if mod(X_max-X_min,X_bin_size)==0
    N_X_bins=floor((X_max-X_min)/X_bin_size);
else
    N_X_bins=floor((X_max-X_min)/X_bin_size)+1;
end
N_X_bins=N_X_bins-1

X_bin_pos=X_min+X_offset+X_bin_size*(0:N_X_bins);
X_bin_pos_remap=X_bin_pos+mid_Xaxis_large-mid_X;
X_bin_middle_pos=X_bin_pos(1:end-1)+X_bin_half_std;



X1_Nradial=X1_Nradial-mid_Xaxis_large+mid_X;
X2_Nradial=X2_Nradial-mid_Xaxis_large+mid_X;
Z_psi_fit_up_small=Z_psi_fit_up(:,X_min:X_max);
Z_psi_fit_down_small=Z_psi_fit_down(:,X_min:X_max);

N_radial_bins=1

% volume of vertical sclices of toroidal dougnnuts shaped bins
% (yummy)



Z_sup_max=Z_psi_fit_up_small(psi_rank_max,:);
% Z_sup_min=Z_psi_fit_up_small(1,:);

% Z_inf_min=Z_psi_fit_down_small(1,:);
Z_inf_max=Z_psi_fit_down_small(psi_rank_max,:);



volume_bin_vertical=zeros(1,N_X_bins);
Npart_binned=zeros(1,N_X_bins);


for X_bin_rank=2:N_X_bins-1
    DV=0;
    for xpos=1:X_bin_size
        DV=DV+mean((Z_sup_max(X_bin_size*(X_bin_rank-1)+xpos:X_bin_size*X_bin_rank+xpos+1)-Z_inf_max(X_bin_size*(X_bin_rank-1)+xpos:X_bin_size*X_bin_rank+xpos+1)))...
            *DX*2*pi*(R0+mean(X_scale(X_bin_pos(X_bin_rank)+xpos-1:X_bin_pos(X_bin_rank)+xpos)));
    end
    volume_bin_vertical(X_bin_rank)=DV;
end

Npart_binned(1,:)=round(volume_bin_vertical(1,:)*Nalpha_max*NPART_DENSITY_RATIO);

Npart_bins_boundaries=0.5*(volume_bin_vertical(2:end)+volume_bin_vertical(1:end-1));


Nalphas_simulated=sum(sum(Npart_binned));

disp('number of particles generated');
Nalphas_simulated

alphas_pos_x=zeros(Nalphas_simulated,1);
alphas_pos_z=zeros(Nalphas_simulated,1);
alphas_pos_phi=zeros(Nalphas_simulated,1);
alphas_Ekin=zeros(Nalphas_simulated,1);
alphas_mm=zeros(Nalphas_simulated,1);
alphas_vpll=zeros(Nalphas_simulated,1);

alphas_pos_psi=interp2(scale_X,scale_Z,psi_norm_XZsmall_map',alphas_pos_x,alphas_pos_z);



Z_index=1;

psi_rank=1;




disp('Number of particles distributed in total volume:');
disp(sum(Npart_binned))


% display particles
figure(1);
set(gca,'FontSize',22);

hold on;
axis xy square

for X_bin_rank=2:N_X_bins-2
    
    %%
    
    X0=X_bin_pos(X_bin_rank)-X_min+1;
    X_bin_half=X_bin_half_std;
    
    %     Z_sup_max_value=max(Z_sup_max(X0-X_bin_half:X0+X_bin_half+1));
    %     Z_sup_pos=Z_inf_min_value;
    %
    %     Z_inf_min_value=min((Z_inf_max(X0-X_bin_half:X0+X_bin_half+1)));
    %     Z_inf_pos=Z_inf_min_value;
    
    x_min=interp1(1:length(X_scale),X_scale,X_bin_pos(X_bin_rank)-X_bin_half,'*linear');
    x_max=interp1(1:length(X_scale),X_scale,X_bin_pos(X_bin_rank)+X_bin_half,'*linear');
    
    z_min=1;
    z_max=-1;
    z=0;
    Npart_rank_end=min(Npart_rank_begin+Npart_binned(X_bin_rank)-1,Nalphas_simulated);
    
    if (Npart_rank_begin<Npart_rank_end)
        Npart_binned(X_bin_rank)
        
        for(n=Npart_rank_begin:Npart_rank_end)
            % random coordinate on [x_min ; x_max[
            %while (z<z_min)||(z>z_max)
            % can be improved here by linear dist
            x=my_rand_linear_dist(Npart_bins_boundaries(X_bin_rank-1),Npart_bins_boundaries(X_bin_rank),x_min,x_max);
            %x=my_rand(1)*(x_max-x_min)+x_min;
            % upper part
            z_max=interp1(X_scale(X_bin_pos(X_bin_rank)-X_bin_half:X_bin_pos(X_bin_rank)+X_bin_half),Z_sup_max(X0-X_bin_half:X0+X_bin_half),x,'*linear');
            % lower part
            z_min=interp1(X_scale(X_bin_pos(X_bin_rank)-X_bin_half:X_bin_pos(X_bin_rank)+X_bin_half),Z_inf_max(X0-X_bin_half:X0+X_bin_half),x,'*linear');
            % lower part
            z=(rand(1)*(z_max-z_min)+z_min);
            % end
            alphas_pos_x(n)=x;
            alphas_pos_z(n)=z;
            
            z_min=1;
            z_max=-1;
            z=0;
        end
        %Npart_rank_begin
        %Npart_rank_end
        figure(1);
        plot(alphas_pos_x(Npart_rank_begin:Npart_rank_end),alphas_pos_z(Npart_rank_begin:Npart_rank_end),'.');
        pause(0.05);
        
    end
    %alphas_pos_psi=interp2(scale_X,scale_Z,psi_norm_XZsmall_map',alphas_pos_x(Npart_rank_begin:Npart_rank_end),alphas_pos_z(Npart_rank_begin:Npart_rank_end));
    %alphas_psi=interp2(scale_X,scale_Z,psi_XZsmall_map',alphas_pos_x(Npart_rank_begin:Npart_rank_end),alphas_pos_z(Npart_rank_begin:Npart_rank_end));
    %alphas_pos_psi=interp1(psi_scale,1:Nradial,alphas_psi);
    %disp('find(alphas_pos_psi<=1)')
    %disp(length(find(alphas_pos_psi<=1)))
    
    Npart_rank_begin=Npart_rank_begin+Npart_binned(X_bin_rank);
    
    
end

% happily shuffle everything
permutation_vector=randperm(Nalphas_simulated);
alphas_pos_x=alphas_pos_x(permutation_vector);
alphas_pos_z=alphas_pos_z(permutation_vector);



%alphas_pos_psi=interp2(scale_X,scale_Z,psi_norm_XZsmall_map',alphas_pos_x,alphas_pos_z);
alphas_psi=interp2(scale_X,scale_Z,psi_XZsmall_map',alphas_pos_x,alphas_pos_z);
alphas_pos_psi=interp1(psi_scale,1:Nradial,alphas_psi);

disp('find(alphas_pos_psi<=1)')
disp(length(find(alphas_pos_psi<=1)))

%%
Nrank=1;

alphas_weight=interp1(1:Nradial,fast_density_profile,alphas_pos_psi)/Ni_avg;
disp('average weight according to radial position');
disp(mean(alphas_weight));


%%


%%
Epart_vector=ones(Nalphas_simulated,1);
weight_vector=ones(Nalphas_simulated,1);

Nrank=1;
%         % to improve statistics, we use a uiniform spread in vtot
%         % and then weight the particles according to their energy

vtot_vector=rand(Nalphas_simulated,1)*max(sqrt(2*(eV/mHe)*energy_H));
Epart_vector=(0.5*(mHe/eV)*vtot_vector.^2);

weight_vector=4*pi*((mHe/(2*pi*Ti_value*eV))^(3/2)).*(vtot_vector.^2).*exp(-Epart_vector/Ti_value);
weight_vector=weight_vector*1.5*Ti_value/mean(Epart_vector.*weight_vector);


vpll_sup_values=sqrt(2*(energy_H+0.5*energy_bin_size)*eV/mHe);

n_vpll_count=0;
Nrank=1;

vpll_vector=ones(1,1);
vpll_max=sqrt(2*(eV/mHe)*Epart_vector);
max_vpll_Erank=vpll_max(1:Nalphas_simulated);
vpll_vector=rand(Nalphas_simulated,1).*(2*max_vpll_Erank)-max_vpll_Erank;


permutation_vector=randperm(Nalphas_simulated);

    

alphas_Ekin=Epart_vector(permutation_vector);
alphas_vpll=vpll_vector(permutation_vector);
alphas_weight_Ekin=weight_vector(permutation_vector);


%%

%     alphas_pos_psi=interp2(scale_X,scale_Z,psi_norm_XZsmall_map',alphas_pos_x,alphas_pos_z);
alphas_momentum=alphas_weight*0;
% alphas_current=interp2(X_scale,Z_scale,JPHI_XZ_map',alphas_pos_x,alphas_pos_z);


disp('average weight according to radial position');
disp(mean(alphas_weight_Ekin));

alphas_weight=alphas_weight.*alphas_weight_Ekin;
disp('average weight according to combined density and Ekin values');
disp(mean(alphas_weight));


% filling up the uniformly distributes toroidal position values
for(n=1:Nalphas_simulated)
    alphas_pos_phi(n)=rand(1)*2*pi;
end

%load('XZsmall_fields_tokamak_pre_collapse.mat', 'Btot_XZ_map');
alphas_mm=(alphas_Ekin-0.5*(mHe*alphas_vpll.^2)/eV)./interp2(scale_X,scale_Z,Btot_XZ_map',alphas_pos_x,alphas_pos_z,'*linear');



%     MEAN_Q1_ROTATION=mean(alphas_momentum);

Nalphas_simulated
alphas_x_rank_tot=interp1(X_scale,1:length(X_scale),alphas_pos_x);

PART_POP_TOT=find((alphas_x_rank_tot>X1_Nradial(psi_rank_max)).*(alphas_x_rank_tot<X2_Nradial(psi_rank_max)));
Nalphas_simulated=length(PART_POP_TOT)
alphas_pos_x_tot=alphas_pos_x(PART_POP_TOT);
alphas_pos_z_tot=alphas_pos_z(PART_POP_TOT);
alphas_pos_phi_tot=alphas_pos_phi(PART_POP_TOT);
alphas_Ekin_tot=alphas_Ekin(PART_POP_TOT);
alphas_mm_tot=alphas_mm(PART_POP_TOT);
alphas_vpll_tot=alphas_vpll(PART_POP_TOT);
alphas_momentum_tot=alphas_momentum(PART_POP_TOT);
alphas_weight_tot=alphas_weight(PART_POP_TOT);

PARTICLES_SPLIT=floor(Nalphas_simulated/NB_PROCESS)

Nalphas_simulated=NB_PROCESS*PARTICLES_SPLIT

disp('press a key to record split files....')
%pause

%%
PART_POP_TOT=(1:Nalphas_simulated);
Nalphas_simulated=length(PART_POP_TOT)
permutation_vector=randperm(Nalphas_simulated);

alphas_pos_x_tot=alphas_pos_x(PART_POP_TOT);
alphas_pos_z_tot=alphas_pos_z(PART_POP_TOT);
alphas_pos_phi_tot=alphas_pos_phi(PART_POP_TOT);
alphas_Ekin_tot=alphas_Ekin(PART_POP_TOT);
alphas_mm_tot=alphas_mm(PART_POP_TOT);
alphas_vpll_tot=alphas_vpll(PART_POP_TOT);
alphas_momentum_tot=alphas_momentum(PART_POP_TOT);
alphas_weight_tot=alphas_weight(PART_POP_TOT);


alphas_pos_x_tot=alphas_pos_x_tot(permutation_vector);
alphas_pos_z_tot=alphas_pos_z_tot(permutation_vector);
alphas_pos_phi_tot=alphas_pos_phi_tot(permutation_vector);
alphas_Ekin_tot=alphas_Ekin_tot(permutation_vector);
alphas_mm_tot=alphas_mm_tot(permutation_vector);
alphas_vpll_tot=alphas_vpll_tot(permutation_vector);
alphas_momentum_tot=alphas_momentum(permutation_vector);
alphas_weight_tot=alphas_weight(permutation_vector);




Nalphas_simulated_tot=Nalphas_simulated

disp('press a key to record split files....')

%%
for PROCESS_NUMBER=1:16
    PART_POP_PLL=PART_POP_TOT((PROCESS_NUMBER-1)*PARTICLES_SPLIT+1:(PROCESS_NUMBER)*PARTICLES_SPLIT);
    
    alphas_pos_x=alphas_pos_x_tot(PART_POP_PLL);
    alphas_pos_z=alphas_pos_z_tot(PART_POP_PLL);
    alphas_pos_phi=alphas_pos_phi_tot(PART_POP_PLL);
    alphas_Ekin=alphas_Ekin_tot(PART_POP_PLL);
    alphas_mm=alphas_mm_tot(PART_POP_PLL);
    alphas_vpll=alphas_vpll_tot(PART_POP_PLL);
    alphas_momentum=alphas_momentum_tot(PART_POP_PLL);
    alphas_weight=alphas_weight_tot(PART_POP_PLL);
    
    Nalphas_simulated=length(PART_POP_PLL)
    
    FILENAME=strcat(EQ_FOLDER,'initial_MBD_',num2str(MB_TEMP),'_f_distribution',num2str(PROCESS_NUMBER),'.mat')
    save (FILENAME,'mHe','ZHe','alphas_weight','alphas_pos_x','alphas_pos_z','alphas_pos_phi','alphas_Ekin','alphas_mm','alphas_vpll','alphas_momentum','Nalphas_simulated');
    
    % save 'speed_bins.mat' Evalues N_energy_bins vpll_range energy_bin_size N_vparallel_bins vpll_bin_size
    % save 'geomtry_bins.mat' density_part_ratio Z_bin_vector Z_volume_bin volume_bin Npart_binned psi_bin_pos N_radial_bins N_X_bins X_bin_size radial_bin_size X_bin_pos Z_psi_fit_up_small Z_psi_fit_down_small
    
    % Nalphas_simulated=length(alphas_vpll);
    
end


