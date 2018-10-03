
    clear all;
    clc
    format compact
    initialize_folder_names_test;
    
    filename=strcat(DATA_FOLDER,'physics_constants.mat');
    load(filename);
    filename=strcat(DATA_FOLDER,'XZsmall_fields_tokamak_pre_collapse.mat');
    load(filename);
    filename=strcat(DATA_FOLDER,'motions_map_dimensions.mat');
    load(filename);
    filename=strcat(DATA_FOLDER,'psi_profiles_kadomtsev.mat');
    load(filename,'psi_star_initial');

    filename=strcat(DATA_FOLDER,'tokamak_PR_map.mat');
    load(filename, 'Rpos_PR_map');
    load(filename, 'radial_r_value_flux');
    filename=strcat(DATA_FOLDER,'B_fields.mat');
    load(filename,'Bpol_PR_map');
    load(filename,'Btor_PR_map');

    filename=strcat(DATA_FOLDER,'TAE_data.mat');
    load(filename);
    
    Btot_PR_map_ini=sqrt(Bpol_PR_map.^2+Btor_PR_map.^2);
%     TAE_angle=2*pi/nTAE
    DPHI=TAE_angle/(NB_PHI-1);
    size_r=size_r_TAE;
    scale_phi=TAE_angle*(0:NB_PHI-1)/(NB_PHI-1);
    scale_psi=pTAE_inf:pTAE_sup;
%     TAE_amplitude=100.0;



%%
% [XX,ZZ] = meshgrid(scale_X,scale_Z);

% F2_PR_map=(Btor_PR_map.*Rpos_PR_map).^2;
% F2_profile=mean(F2_PR_map(1:NP-1,:),1);
% P_prime=gradient(P_initial_profile,psi_scale);
% F2_prime=gradient(F2_profile,psi_scale);

% for(p=1:NP)
%     for(r=1:Nradial)
%         Pprime_PR_map(p,r)=P_prime(r);
%         FFprime_PR_map(p,r)=0.5*F2_prime(r);
%     end
% end

% NX=400;
% NZ=4*NX;
% finesse_data_X=reshape((Rpos_PR_map(:,:)-R0),NP*Nradial,1);
% finesse_data_Z=reshape(Z_PR_map(:,:),NP*Nradial,1);
Rpos_map=zeros(size_X,size_Z);

for (x=1:size_X)
    for (z=1:size_Z)
%        Rpos_map(x,z)=R0+(x-mid_X)*DX;
         Rpos_map(:,z)=scale_X+R0;
    end
end



gpsi_X=zeros(size_X,size_Z);
gpsi_Z=zeros(size_X,size_Z);

for (x=3:size_X-2)
    for (z=3:size_Z-2)
        gpsi_X(x,z)=(1/12)*(-psi_XZsmall_map(x+2,z)+psi_XZsmall_map(x-2,z))+(2/3)*(psi_XZsmall_map(x+1,z)-psi_XZsmall_map(x-1,z));
        gpsi_Z(x,z)=(1/12)*(-psi_XZsmall_map(x,z+2)+psi_XZsmall_map(x,z-2))+(2/3)*(psi_XZsmall_map(x,z+1)-psi_XZsmall_map(x,z-1));
    end
end
gpsi_X=gpsi_X/DX;
gpsi_Z=gpsi_Z/DX;

% here you can verify the good consistency
% of the psi XZ map and the BZ map
BZ_recalc=gpsi_X./Rpos_map;
BX_recalc=-gpsi_Z./Rpos_map;




figure(6);
imagesc(scale_X,scale_Z,((BpolX_initial_XZsmall_map-(BX_recalc))/mean(mean(abs(BX_recalc))))',[-0.01 0.01]);
axis xy;
colorbar;
title('B_X FINESSE error')

figure(5);
imagesc(scale_X,scale_Z,((BpolZ_initial_XZsmall_map-BZ_recalc)/mean(mean(abs(BZ_recalc))))',[-0.01 0.01]);
axis xy;
colorbar;
title('B_Z FINESSE error')

