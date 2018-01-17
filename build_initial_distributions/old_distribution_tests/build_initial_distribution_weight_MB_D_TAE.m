clear all;
initialize_folder_names;
reset_data_analysis_environment;
format compact

filename=strcat(DATA_FOLDER,'flux_geometry.mat');
load(filename);
filename=strcat(DATA_FOLDER,'physics_constants.mat');
load(filename);
filename=strcat(DATA_FOLDER,'pressure_profile.mat');
load(filename);
filename=strcat(DATA_FOLDER,'motions_map_dimensions.mat');
load(filename);
filename=strcat(DATA_FOLDER,'tokamak_PR_map.mat');
load(filename);

TOROIDAL_FIELD_DIRECTION=sign(mean(mean(Bphi_XZsmall_map)))

rho_tor=(tor_flux_profile/max(tor_flux_profile));
nH0=1.44131e17
fast_density_profile=nH0*0.521298*exp(-(0.198739/0.298228)*tanh((sqrt(rho_tor)-0.49123)/0.198739));
radial_bin_size=8;
radial_bin_half_size=0.5*(radial_bin_size);
radial_bins=[radial_bin_half_size+1:radial_bin_size:Nradial-4]';
N_radial_bins=size(radial_bins,1)
Nalpha_binned=ones(N_radial_bins,1);
mHe=mD
ZHe=1

density_on_axis=max(fast_density_profile)

MB_TEMP=250  % keV
energy_max=4500*1e3; %keV

% energy_D_range=energy_H_range;
% f_D_energy_percentage=f_H_energy_percentage;
% save f_D_energy_percentage_map.mat  radial_bins energy_D_range f_D_energy_percentage



%%
filename=strcat(DATA_FOLDER,'volume_flux_geometry.mat');
load(filename);
filename=strcat(DATA_FOLDER,'XZsmall_fields_tokamak_pre_collapse.mat');
load(filename);

Npart_simulated=2e6;
NPART_DENSITY_RATIO=1.69e-13


psi_rank_max=470

Ni_avg=mean(fast_density_profile(1:psi_rank_max))

close all


%Normalize the alpha density
Nalpha_per_m3=5.0e4


% loop through all flux surface bins
psi_bin_pos=1+radial_bin_half_size+radial_bin_size*(0:N_radial_bins-1)

%
X_min=mid_Xaxis_large-0.5*size(scale_X,2)



% display options
figure(1);
set(gca,'FontSize',22);

hold on;
axis xy square


%should be an even number
X_bin_size=10;
X_bin_half_std=round(0.5*(X_bin_size));

X_min=X1_Nradial(psi_rank_max+1.3*X_bin_size);
X_max=X2_Nradial(psi_rank_max+2*X_bin_size);

if mod(X_max-X_min,X_bin_size)==0
    N_X_bins=floor((X_max-X_min)/X_bin_size);
else
    N_X_bins=floor((X_max-X_min)/X_bin_size)+1;
end
N_X_bins=N_X_bins-1

X_bin_pos=X_min+X_bin_size*(0:N_X_bins);



N_radial_bins=1

% volume of vertical sclices of toroidal dougnnuts shaped bins
% (yummy)

Z_sup_max=Z_psi_fit_up(psi_rank_max,:);
Z_inf_max=Z_psi_fit_down(psi_rank_max,:);

Z_sup_max_shift=Z_psi_fit_up(psi_rank_max,X_min:X_max);
Z_inf_max_shift=Z_psi_fit_down(psi_rank_max,X_min:X_max);


volume_bin_vertical=zeros(1,N_X_bins);


for X_bin_rank=1:N_X_bins-1
    DV=0;
    
    for xpos=1:X_bin_size
        Xscpos=X_bin_size*(X_bin_rank-1)+xpos;
        DeltaZ=Z_sup_max_shift(Xscpos)-Z_inf_max_shift(Xscpos);
        DV=DV+DeltaZ*DX*2*pi*(R0+X_scale(Xscpos));
    end
    volume_bin_vertical(X_bin_rank)=DV;
end

Npart_binned(1,:)=round(volume_bin_vertical(1,:)*Nalpha_per_m3);
Npart_binned=Npart_simulated*Npart_binned/sum(Npart_binned);
Npart_binned=round(Npart_binned);

Npart_bins_boundaries=volume_bin_vertical*0;
Npart_bins_boundaries(2:end)=0.5*(volume_bin_vertical(2:end)+volume_bin_vertical(1:end-1));


Nalphas_simulated=round(sum(sum(Npart_binned)));

disp('number of particles generated');
Nalphas_simulated

alphas_pos_x=zeros(Nalphas_simulated,1);
alphas_pos_z=zeros(Nalphas_simulated,1);
alphas_pos_phi=zeros(Nalphas_simulated,1);
alphas_Ekin=zeros(Nalphas_simulated,1);
alphas_mm=zeros(Nalphas_simulated,1);
alphas_vpll=zeros(Nalphas_simulated,1);

alphas_pos_psi=interp2(scale_X,scale_Z,psi_norm_XZsmall_map',alphas_pos_x,alphas_pos_z);


psi_rank=1;
Npart_rank_begin=1;

disp('Number of particles distributed in total volume:');
disp(sum(Npart_binned))


% display particles
figure(1);
set(gca,'FontSize',22);

hold on;
axis xy square

for X_bin_rank=1:N_X_bins-2
    %%
    
    x_min=interp1(1:length(X_scale),X_scale,X_bin_pos(X_bin_rank),'*PCHIP');
    x_max=interp1(1:length(X_scale),X_scale,X_bin_pos(X_bin_rank+1),'*PCHIP');
    
    z_min=1;
    z_max=-1;
    z=0;
    Npart_rank_end=min(Npart_rank_begin+Npart_binned(X_bin_rank)-1,Nalphas_simulated);
    
    if (Npart_rank_begin<Npart_rank_end)
        Npart_binned(X_bin_rank)
        
        for(n=Npart_rank_begin:Npart_rank_end)
            % random coordinate on [x_min ; x_max[
            % can be improved here by linear dist
            x=my_rand_linear_dist(Npart_bins_boundaries(X_bin_rank),Npart_bins_boundaries(X_bin_rank+1),x_min,x_max);
            % upper part
            z_max=interp1(X_scale(X_bin_pos(X_bin_rank):X_bin_pos(X_bin_rank+1)),Z_sup_max(X_bin_pos(X_bin_rank):X_bin_pos(X_bin_rank+1)),x,'*linear');
            % lower part
            z_min=interp1(X_scale(X_bin_pos(X_bin_rank):X_bin_pos(X_bin_rank+1)),Z_inf_max(X_bin_pos(X_bin_rank):X_bin_pos(X_bin_rank+1)),x,'*linear');
            % position between upper and lower boundary
            z=(rand(1)*(z_max-z_min)+z_min);
            alphas_pos_x(n)=x;
            alphas_pos_z(n)=z;
            
            z_min=1;
            z_max=-1;
            z=0;
        end

        figure(1);
        plot(alphas_pos_x(Npart_rank_begin:Npart_rank_end),alphas_pos_z(Npart_rank_begin:Npart_rank_end),'.');
        pause(0.05);
        
    end
    
    Npart_rank_begin=Npart_rank_begin+Npart_binned(X_bin_rank);
    
    
end
%%
% happily shuffle everything
permutation_vector=randperm(Nalphas_simulated);
alphas_pos_x=alphas_pos_x(permutation_vector);
alphas_pos_z=alphas_pos_z(permutation_vector);



alphas_psi=interp2(scale_X,scale_Z,psi_XZsmall_map',alphas_pos_x,alphas_pos_z);
alphas_pos_psi=interp1(psi_scale,1:Nradial,alphas_psi);

disp('find(alphas_pos_psi<=1)')
disp(length(find(alphas_pos_psi<=1)))

%%

alphas_weight_density=interp1(1:Nradial,fast_density_profile,alphas_pos_psi)/Ni_avg;
disp('average weight according to radial position');
disp(mean(alphas_weight_density));




%%
Epart_vector=ones(Nalphas_simulated,1);
weight_vector=ones(Nalphas_simulated,1);

energy_bin_size=10.0*1e3; %keV
energy_H(1)=0.5*energy_bin_size;
N_energy_bins=round(energy_max/energy_bin_size)
for (x=2:N_energy_bins)
    energy_H(x)= energy_H(x-1)+energy_bin_size;
end

Ti_value=MB_TEMP*1e3

% to improve statistics, we use a uiniform spread in vtot
% and then weight the particles according to their energy

vtot_vector=rand(Nalphas_simulated,1)*max(sqrt(2*(eV/mHe)*energy_H));
Epart_vector=(0.5*(mHe/eV)*vtot_vector.^2);

weight_vector=4*pi*((mHe/(2*pi*Ti_value*eV))^(3/2)).*(vtot_vector.^2).*exp(-Epart_vector/Ti_value);
weight_vector=weight_vector*1.5*Ti_value/mean(Epart_vector.*weight_vector);


vpll_sup_values=sqrt(2*(energy_H+0.5*energy_bin_size)*eV/mHe);


vpll_vector=ones(1,1);
vpll_max=sqrt(2*(eV/mHe)*Epart_vector);
max_vpll_Erank=vpll_max(1:Nalphas_simulated);
vpll_vector=rand(Nalphas_simulated,1).*(2*max_vpll_Erank)-max_vpll_Erank;


permutation_vector=randperm(Nalphas_simulated);

    
alphas_pos_x=alphas_pos_x(permutation_vector);
alphas_pos_z=alphas_pos_z(permutation_vector);
alphas_Ekin=Epart_vector(permutation_vector);
alphas_vpll=vpll_vector(permutation_vector);
alphas_weight_Ekin=weight_vector(permutation_vector);
alphas_weight_density=alphas_weight_density(permutation_vector);

alphas_psi=interp2(scale_X,scale_Z,psi_XZsmall_map',alphas_pos_x,alphas_pos_z);
alphas_pos_psi=interp1(psi_scale,1:Nradial,alphas_psi);


%%

alphas_momentum=alphas_Ekin*0;

disp('average weight according to Energy');
disp(mean(alphas_weight_Ekin));
disp('MB temperature (eV)');
disp((2/3)*mean(alphas_weight_Ekin.*alphas_Ekin));

alphas_weight=alphas_weight_density.*alphas_weight_Ekin;
disp('average weight according to combined density and Ekin values');
disp(mean(alphas_weight));


% filling up the uniformly distributes toroidal position values
for(n=1:Nalphas_simulated)
    alphas_pos_phi(n)=rand(1)*2*pi;
end

alphas_mm=(alphas_Ekin-0.5*(mHe*alphas_vpll.^2)/eV)./interp2(scale_X,scale_Z,Btot_XZ_map',alphas_pos_x,alphas_pos_z,'*cubic');

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
alphas_weight_density_tot=alphas_weight_density(PART_POP_TOT);
alphas_weight_Ekin_tot=alphas_weight_Ekin(PART_POP_TOT);



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
alphas_weight_density_tot=alphas_weight_density(PART_POP_TOT);
alphas_weight_Ekin_tot=alphas_weight_Ekin(PART_POP_TOT);


alphas_pos_x_tot=alphas_pos_x_tot(permutation_vector);
alphas_pos_z_tot=alphas_pos_z_tot(permutation_vector);
alphas_pos_phi_tot=alphas_pos_phi_tot(permutation_vector);
alphas_Ekin_tot=alphas_Ekin_tot(permutation_vector);
alphas_mm_tot=alphas_mm_tot(permutation_vector);
alphas_vpll_tot=alphas_vpll_tot(permutation_vector);
alphas_momentum_tot=alphas_momentum(permutation_vector);
alphas_weight_tot=alphas_weight(permutation_vector);
alphas_weight_density_tot=alphas_weight_density(permutation_vector);
alphas_weight_Ekin_tot=alphas_weight_Ekin(permutation_vector);


Nalphas_simulated_tot=Nalphas_simulated

alphas_pos_x=alphas_pos_x_tot;
alphas_pos_z=alphas_pos_z_tot;
alphas_pos_phi=alphas_pos_phi_tot;
alphas_Ekin=alphas_Ekin_tot;
alphas_mm=alphas_mm_tot;
alphas_vpll=alphas_vpll_tot;
alphas_momentum=alphas_momentum_tot;
alphas_weight=alphas_weight_tot;

%% Dealing with a gyro-center position is better



% direction of the B field at local positions
bX=interp2(scale_X,scale_Z,bX_XZ_map',alphas_pos_x,alphas_pos_z);
bZ=interp2(scale_X,scale_Z,bZ_XZ_map',alphas_pos_x,alphas_pos_z);
bphi=TOROIDAL_FIELD_DIRECTION*sqrt(1-(bX.^2+bZ.^2));

alphas_Bfield=interp2(scale_X,scale_Z,Btot_XZ_map',alphas_pos_x,alphas_pos_z);

alphas_Eperp=(alphas_Ekin-0.5*(mHe*alphas_vpll.^2)/eV);
alphas_vperp=sqrt((eV/mHe)*2*alphas_Eperp);

% attribute gyro phase if not included in initial data
norm_angle=rand(Nalphas_simulated,1)*2*pi;
uX=sqrt(1./(1+(bX./bphi).^2));
uZ=0*uX;
uphi=-(bX./bphi).*uX;
unorm=sqrt(uX.^2+uZ.^2+uphi.^2);
wX=1./sqrt((uX./uphi).^2+1+((1./bZ).*(uX./uphi).*(bphi-bX)).^2);
wZ=(wX./bZ).*((uX./uphi).*bphi-bX);
wphi=-(uX./uphi).*wX;
wnorm=sqrt(wX.^2+wZ.^2+wphi.^2);
wX=wX./wnorm;
wZ=wZ./wnorm;
wphi=wphi./wnorm;
% normal vector N = (cos) u + (sin) w
% that N.b=0 precisely
NX=(cos(norm_angle).*uX+sin(norm_angle).*wX);
NZ=((cos(norm_angle).*uZ+sin(norm_angle).*wZ));
Nphi=(cos(norm_angle).*uphi+sin(norm_angle).*wphi);


v0_X=alphas_vpll.*bX+alphas_vperp.*NX;
v0_Z=alphas_vpll.*bZ+alphas_vperp.*NZ;
v0_phi=alphas_vpll.*bphi+alphas_vperp.*Nphi;

v0_X=double(v0_X);
v0_Z=double(v0_Z);
v0_phi=double(v0_phi);

v_X=v0_X;
v_Z=v0_Z;
v_phi=v0_phi;



pos_X_gc=(mHe/eV)*(1/ZHe)*(v_Z.*bphi-v_phi.*bZ)./alphas_Bfield;
pos_Z_gc=(mHe/eV)*(1/ZHe)*(v_phi.*bX-v_X.*bphi)./alphas_Bfield;


%% Moving the particles on their gyro-orbits
alphas_pos_x_gc=alphas_pos_x;
alphas_pos_z_gc=alphas_pos_z;

alphas_pos_x=alphas_pos_x_gc-pos_X_gc;
alphas_pos_z=alphas_pos_z_gc-pos_Z_gc;



FILENAME=strcat(EQ_FOLDER,'initial_MBD_',num2str(MB_TEMP),'_f_distribution_all.mat')
save (FILENAME,'mHe','ZHe','alphas_weight','alphas_weight_density','alphas_weight_Ekin','alphas_pos_x','alphas_pos_z','alphas_pos_phi','v_X','v_Z','v_phi','alphas_Ekin','alphas_mm','alphas_vpll','alphas_momentum','Nalphas_simulated');


