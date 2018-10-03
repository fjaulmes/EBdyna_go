


% psi_data=reshape(psiH_PR_map(:,1:NR),NB_THETA*NR,1);
% % psi_XZ_map=griddata(finesse_data_X_extended,finesse_data_Z_extended,psi_data,XX_small,ZZ_small,'cubic');
% psiH_XZ_map=dtinterp(finesse_mesh_extended,finesse_mesh_extended_dtri,psiH_PR_map(IDTRI_EXT),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
% psiH_XZ_map(isnan(psiH_XZ_map))=0;
% psiH_XZ_map=psiH_XZ_map';
psi2D=psiH_XZ_map;


gpsi_R=zeros(NZ,NZ);
gpsi_Z=zeros(NZ,NZ);

for (x=3:NZ-2)
    for (z=3:NZ-2)
        gpsi_R(x,z)=(1/12)*(-psi2D(x+2,z)+psi2D(x-2,z))+(2/3)*(psi2D(x+1,z)-psi2D(x-1,z));
        gpsi_Z(x,z)=(1/12)*(-psi2D(x,z+2)+psi2D(x,z-2))+(2/3)*(psi2D(x,z+1)-psi2D(x,z-1));
%         gpsi_Z(x,z)=0.5*(psi2D(x,z+1)-psi2D(x,z-1));
    end
end
gpsi_R=gpsi_R/DX;
gpsi_Z=gpsi_Z/DX;

BHpolX_XZ_map_recalc=-gpsi_Z./Rpos_map;
BHpolZ_XZ_map_recalc=gpsi_R./Rpos_map;

BH_data=reshape(BHpolX_XZ_map_recalc(:,:)',NZ*NZ,1);
BHpolX_PR_map_recalc=griddata(X_scale_data,Z_scale_data,BH_data,Rpos_PR_map-R0,Z_PR_map,'cubic');

BH_data=reshape(BHpolZ_XZ_map_recalc(:,:)',NZ*NZ,1);
BHpolZ_PR_map_recalc=griddata(X_scale_data,Z_scale_data,BH_data,Rpos_PR_map-R0,Z_PR_map,'cubic');
        

figure(8);
imagesc((BHpol_X_map-BHpolX_XZ_map_recalc)',[-0.01 0.01])
axis xy
