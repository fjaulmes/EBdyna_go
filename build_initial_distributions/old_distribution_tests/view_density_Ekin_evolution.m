%close all;

% load('map_dimensions.mat')
% load('geomtry_bins.mat')
% load('speed_bins.mat')
% load('XZsmall_fields_tokamak.mat')

load('physics_constants.mat')
load('initial_alphas_distributions.mat');
load('initial_alphas_pre_collapse');

disp('Total number of particles simulated');
disp(Nalphas_simulated);

load volume_flux_geometry.mat

% rescaling X bins sizes
X_min=mid_Xaxis_large-0.5*size(scale_X,2);
X_max=mid_Xaxis_large+0.5*size(scale_X,2)-1;

X1_Nradial=X1_Nradial-mid_Xaxis_large+mid_X;
X2_Nradial=X2_Nradial-mid_Xaxis_large+mid_X;

X_bin_size=8;
X_bin_half=0.5*(X_bin_size);
X_min=X1_Nradial(simulation_size_r);
X_max=X2_Nradial(simulation_size_r);

X_offset=-round(0.5*mod(X_max-X_min,X_bin_size))+1+X_bin_half;
if mod(X_max-X_min,X_bin_size)==0
    N_X_bins=ceil((X_max-X_min)/X_bin_size);
else
    N_X_bins=ceil((X_max-X_min)/X_bin_size)+1;
end

X_bin_pos=X_min+X_offset+X_bin_size*(0:N_X_bins-1)
X_bin_pos_remap=X_bin_pos+mid_Xaxis_large-mid_X;

N_Z_bins=round(1.4*N_X_bins);
Z_bin_pos=mid_Z+(-round(0.5*N_Z_bins):round(0.5*N_Z_bins))*X_bin_size;
N_Z_bins=size(Z_bin_pos,2);

X0=alphas_pos_x;
Z0=alphas_pos_x;
vpll0=alphas_vpll;
Ekin0=alphas_Ekin;
mm0=alphas_mm;
Bfield0=interp2(scale_X,scale_Z,Btot_XZ_map',alphas_pos_x,alphas_pos_z,'*linear');
psi_pos0=interp2(scale_X,scale_Z,radial_XZsmall_map',alphas_pos_x,alphas_pos_z,'*nearest');


DX_bin_size=X_bin_size*DX;
DX_bin_size=round(0.5*X_bin_size)*DX;


for (f=2:2:24)
     
    filename='.\';
    frame_name=strcat('alphas_record_',num2str(f));
    filename=strcat(filename,frame_name,'000.mat');
    load(filename);
    
    alphas_pos_psi=interp2(scale_X,scale_Z,radial_XZsmall_map',alphas_pos_x,alphas_pos_z,'*nearest');
    
    clear vparallel_vector density_part_vector Xpos_vector Zpos_vector
    
    bin_index=1;
    
    % toroidal square bins density values
    for(z=1:2:N_Z_bins)
        %     z_pos=mid_Z-mid_X+X_bin_pos(z);
        z_pos=Z_bin_pos(z);
        for(x=1:2:N_X_bins)
            indexes_volume_bin=(alphas_pos_z>=scale_Z(z_pos)-DX_bin_size);
            indexes_volume_bin=indexes_volume_bin.*(alphas_pos_z<scale_Z(z_pos)+DX_bin_size);
            indexes_volume_bin=indexes_volume_bin.*(alphas_pos_x>=scale_X(X_bin_pos(x))-DX_bin_size);
            indexes_volume_bin=indexes_volume_bin.*(alphas_pos_x<scale_X(X_bin_pos(x))+DX_bin_size);
            Npart_axis_bin=size(find(indexes_volume_bin),1);
            Rpos=Rpos_XZsmall_map(X_bin_pos(x),z_pos);
            volume_axis_bin=2*pi*Rpos*(2*DX_bin_size)*(2*DX_bin_size);
            %Ekin_vector(bin_index)=mean(vExB(find(indexes_volume_bin)));
            lambda_vector(bin_index)=mean(alphas_mm(find(indexes_volume_bin))./alphas_Ekin(find(indexes_volume_bin)));
            density_part_vector(bin_index)=density_part_ratio*Npart_axis_bin/volume_axis_bin;
            Xpos_vector(bin_index)=scale_X(X_bin_pos(x));
            Zpos_vector(bin_index)=scale_Z(z_pos);
            bin_index=bin_index+1;
        end
    end
    
    
    
    [XX_psi ZZ_psi]=meshgrid(scale_X(X_bin_pos(1:2:end)),scale_Z(Z_bin_pos(1:2:end)));
    alphas_density_map=griddata(Xpos_vector,Zpos_vector,density_part_vector,XX_psi,ZZ_psi);
    alphas_density_map=alphas_density_map';
    alphas_density_map(isnan(alphas_density_map))=0;
    alphas_density_map=smooth_small_map(alphas_density_map);
    
    alphas_lambda_map=griddata(Xpos_vector,Zpos_vector,lambda_vector,XX_psi,ZZ_psi);
    alphas_lambda_map=alphas_lambda_map';
    alphas_lambda_map(isnan(alphas_lambda_map))=0;
    alphas_lambda_map=smooth_small_map(alphas_lambda_map);
    
%     alphas_Ekin_map=griddata(Xpos_vector,Zpos_vector,Ekin_vector,XX_psi,ZZ_psi);
%     alphas_Ekin_map=alphas_Ekin_map';
%     alphas_Ekin_map(isnan(alphas_Ekin_map))=0;
%     alphas_Ekin_map=smooth_small_map(alphas_Ekin_map);
    
    
    figure(1);
    imagesc(scale_X(X_bin_pos(1:2:end)),scale_Z(Z_bin_pos(1:2:end)),alphas_density_map',[0.1 0.5]*1e16);
    axis xy square; colorbar;
    hold on;
    titre=strcat('density of alphas at time=',num2str(time))
    title(titre);
    xlabel('X (m)');
    ylabel('Z (m)');
    

    figure(2);
    imagesc(scale_X(X_bin_pos(1:2:end)),scale_Z(Z_bin_pos(1:2:end)),alphas_lambda_map',[0.05 0.5]);
    axis xy square; colorbar;
    hold on;
    titre=strcat('averaged \lambda of alphas at time=',num2str(time));
    title(titre);
    xlabel('X (m)');
    ylabel('Z (m)');
    
    pause(0.1);
    
end
