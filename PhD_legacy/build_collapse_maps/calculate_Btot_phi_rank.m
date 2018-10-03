

BpolX_XZ_map=zeros(sizeX,sizeZ);
BpolZ_XZ_map=zeros(sizeX,sizeZ);

psi_PR_map(:,1:psi_max)=psiH_PR_map(:,1:psi_max)+psi_star_PR_map(:,1:psi_max);
%psi_data=reshape(psi_PR_map(:,1:NR),NB_THETA*NR,1);
% psi_XZ_map=griddata(finesse_data_X_extended,finesse_data_Z_extended,psi_data,XX_small,ZZ_small,'linear');
%psi_XZ_map=dtinterp(finesse_mesh_extended,finesse_mesh_extended_dtri,psi_PR_map(IDTRI_EXT),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
%psi_XZ_map(isnan(psi_XZ_map))=0;
%psi_XZ_map=psi_XZ_map';
%psi2D=psi_XZ_map;

% calculate Bstar THEN add up BH
% much better precision

psi_data=reshape(psi_star_PR_map(:,1:size_r),NB_THETA*size_r,1);
psi_star_XZ_map=gridfit(finesse_data_X,finesse_data_Z,psi_data,scale_X,scale_Z,'smoothness',INTERP_SMOOTHNESS);
psi_star_XZ_map=psi_star_XZ_map';
psi2D=psi_star_XZ_map;


% psi_data=reshape(psi_PR_map(:,1:NR),NB_THETA*NR,1);
% % psi_XZ_map=gridfit(finesse_data_X_extended,finesse_data_X_extended,psi_data,scale_X,scale_Z,'smoothness',0.5);
% psi_XZ_map=gridfit(finesse_data_X_extended,finesse_data_Z_extended,psi_data,scale_X,scale_Z);
% psi_XZ_map=psi_XZ_map'; 
% psi2D=psi_XZ_map;

% if CALCULATE_BSTAR_DATA==1
%     psi_star_omega_map_ini(:,:)=rotate_map_phi(psi_star_omega_map_rank_ini,phi_rank);
%     psi_star_PR_map_ini=psi_star_omega_map_ini';
%     psi_star_PR_map_ext=0*psi_PR_map;
%     psi_star_PR_map_ext(:,1:psi_max)=psi_star_PR_map(:,1:psi_max)-psi_star_PR_map_ini(:,1:psi_max);
% %    psi_data=reshape(psi_star_PR_map(:,1:psi_max),NB_THETA*psi_max,1);
% %      psi_star_XZ_map=dtinterp(finesse_mesh_extended,finesse_mesh_extended_dtri,psi_star_PR_map_ext(IDTRI_EXT),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
% %      psi_star_XZ_map(isnan(psi_star_XZ_map))=0;
% %      psi_star_XZ_map=psi_star_XZ_map';
%     %psi_star_XZ_map=psi_XZ_map-psi_XZ_map_ini;
%     
%     psi_data=reshape(psi_star_PR_map_ext(:,1:NR),NB_THETA*NR,1);
%     psi_star_XZ_map=gridfit(finesse_data_X_extended,finesse_data_Z_extended,psi_data,scale_X,scale_Z);
%     psi_star_XZ_map=psi_star_XZ_map'; 
% end

gpsi_R=zeros(sizeX,sizeZ);
gpsi_Z=zeros(sizeX,sizeZ);
if CALCULATE_BSTAR_DATA==1
     gpsistar_R=zeros(sizeX,sizeZ);
     gpsistar_Z=zeros(sizeX,sizeZ);
end

VZ=(3:sizeZ-2);

if CALCULATE_BSTAR_DATA==1
   for (x=3:sizeX-2)

%     gpsi_R(x,3:sizeZ-2)=(1/12)*(-psi2D(x+2,3:sizeZ-2)+psi2D(x-2,3:sizeZ-2))+(2/3)*(psi2D(x+1,3:sizeZ-2)-psi2D(x-1,3:sizeZ-2));
%     gpsi_Z(x,3:sizeZ-2)=(1/12)*(-psi2D(x,(3:sizeZ-2)-2)+psi2D(x,(3:sizeZ-2)-2))+(2/3)*(psi2D(x,(3:sizeZ-2)+1)-psi2D(x,(3:sizeZ-2)-1));
       gpsistar_R(x,VZ)=(1/12)*(-psi_star_XZ_map(x+2,VZ)+psi_star_XZ_map(x-2,VZ))+(2/3)*(psi_star_XZ_map(x+1,VZ)-psi_star_XZ_map(x-1,VZ));
       gpsistar_Z(x,VZ)=(1/12)*(-psi_star_XZ_map(x,VZ+2)+psi_star_XZ_map(x,VZ-2))+(2/3)*(psi_star_XZ_map(x,VZ+1)-psi_star_XZ_map(x,VZ-1));

%     %         gpsi_Z(x,z)=0.5*(psi2D(x,z+1)-psi2D(x,z-1));
   end
end

for (x=3:sizeX-2)
    gpsi_R(x,VZ)=(1/12)*(-psi2D(x+2,VZ)+psi2D(x-2,VZ))+(2/3)*(psi2D(x+1,VZ)-psi2D(x-1,VZ));
    gpsi_Z(x,VZ)=(1/12)*(-psi2D(x,VZ+2)+psi2D(x,VZ-2))+(2/3)*(psi2D(x,VZ+1)-psi2D(x,VZ-1));
end
% for (x=3:sizeX-2)
%     for (z=3:sizeZ-2)
%         gpsi_R(x,z)=(1/12)*(-psi2D(x+2,z)+psi2D(x-2,z))+(2/3)*(psi2D(x+1,z)-psi2D(x-1,z));
%         gpsi_Z(x,z)=(1/12)*(-psi2D(x,z+2)+psi2D(x,z-2))+(2/3)*(psi2D(x,z+1)-psi2D(x,z-1));
% %         gpsi_Z(x,z)=0.5*(psi2D(x,z+1)-psi2D(x,z-1));
%     end
% end
gpsi_R=gpsi_R/DX;
gpsi_Z=gpsi_Z/DX;
if CALCULATE_BSTAR_DATA==1
     gpsistar_R=gpsistar_R/DX;
     gpsistar_Z=gpsistar_Z/DX;
end

BpolX_XZ_map=-gpsi_Z./Rpos_XZsmall_map;
BpolZ_XZ_map=gpsi_R./Rpos_XZsmall_map;
if CALCULATE_BSTAR_DATA==1
     BstarX_XZ_map=-gpsistar_R./Rpos_XZsmall_map;
     BstarZ_XZ_map=gpsistar_Z./Rpos_XZsmall_map;
end

%%

BpolX_XZ_map=BHpolX_initial_XZsmall_map+BpolX_XZ_map;
BpolZ_XZ_map=BHpolZ_initial_XZsmall_map+BpolZ_XZ_map;
%removing wrong boundaries
BpolX_XZ_map=psi_XZ_map_mask.*BpolX_XZ_map;
BpolZ_XZ_map=psi_XZ_map_mask.*BpolZ_XZ_map;
if CALCULATE_BSTAR_DATA==1
   BstarX_XZ_map=psi_XZ_map_mask.*BstarX_XZ_map;
   BstarZ_XZ_map=psi_XZ_map_mask.*BstarZ_XZ_map;
end

if CALCULATE_PR_DATA_FILE==1
    
    Bstar_data=reshape(BpolX_XZ_map(:,:)',sizeX*sizeZ,1);
    BpolX_PR_map=griddata(X_scale_data,Z_scale_data,Bstar_data,RR,Z_PR_map(:,1:size_r),'cubic');
    
    Bstar_data=reshape(BpolZ_XZ_map(:,:)',sizeX*sizeZ,1);
    BpolZ_PR_map=griddata(X_scale_data,Z_scale_data,Bstar_data,RR,Z_PR_map(:,1:size_r),'cubic');
    
    Btot_PR_map=sqrt(BpolX_PR_map.^2+BpolZ_PR_map.^2+Btor_PR_map(:,1:size_r).^2);
    
end

% BpolX_PR_map=BHpol_X_PR_map_size_r+BstarX_PR_map;
% BpolZ_PR_map=BHpol_Z_PR_map_size_r+BstarZ_PR_map;

% B_data=reshape(BpolX_PR_map(:,1:size_r),NP*size_r,1);
% BpolX_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,B_data(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
% BpolX_XZ_map=BpolX_XZ_map';
% BpolX_XZ_map(isnan(BpolX_XZ_map))=0;
% 
% B_data=reshape(BpolZ_PR_map(:,1:size_r),NP*size_r,1);
% BpolZ_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,B_data(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
% BpolZ_XZ_map=BpolZ_XZ_map';
% BpolZ_XZ_map(isnan(BpolZ_XZ_map))=0;
% 
% BpolX_XZ_map=BpolX_XZ_map.*(psi_XZsmall_map<xpos_sep)+(psi_XZsmall_map>=xpos_sep).*BpolX_initial_XZsmall_map;
% BpolZ_XZ_map=BpolZ_XZ_map.*(psi_XZsmall_map<xpos_sep)+(psi_XZsmall_map>=xpos_sep).*BpolZ_initial_XZsmall_map;
%
% if CALCULATE_VD_DATA_FILE==1
    Btot_XZ_map=sqrt(BpolX_XZ_map.^2+BpolZ_XZ_map.^2+Bphi_XZsmall_map.^2);
    if CALCULATE_BSTAR_DATA==1
        Bstartot_XZ_map=sqrt(BstarX_XZ_map.^2+BstarZ_XZ_map.^2);
    end
% end

% Btot_PR_map(:,1)=Btot_PR_map(:,2);

energy_mag=0;
energy_star_mag=0;
for (x=Xinf+1:Xsup-1)
    Z_inf=round(interp1(scale_Z,1:length(scale_Z),Z_psi_fit_down(psi_max,x)));
    Z_sup=round(interp1(scale_Z,1:length(scale_Z),Z_psi_fit_up(psi_max,x)));
    for (z=Z_inf:Z_sup)
        energy_mag=energy_mag+0.5*(Btot_XZ_map(x-Xinf,z)^2)*(scale_X(x-Xinf)+R0)*(DPHI*DX*DX/mu0);
       if CALCULATE_BSTAR_DATA==1
            energy_star_mag=energy_star_mag+0.5*(Bstartot_XZ_map(x-Xinf,z)^2)*(scale_X(x-Xinf)+R0)*(DPHI*DX*DX/mu0);
       end
    end
end

energy_phi_evol(round(frame_rank),phi_index)=energy_mag;
energy_star_phi_evol(round(frame_rank),phi_index)=energy_star_mag;