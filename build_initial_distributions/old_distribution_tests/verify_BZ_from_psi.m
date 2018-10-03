% load('motions_map_dimensions.mat')
% load('XZsmall_fields_tokamak_post_collapse.mat')
% load('post_collapse_XZsmall_maps.mat')
BX_XZ_map=BpolX_final_XZsmall_map;
BZ_XZ_map=BpolZ_final_XZsmall_map;

% BZ_XZ_map=BZ_XZ_map.*mask_XZ_map;
psi2D=psi_XZsmall_map;
Rpos_map=Rpos_XZsmall_map;

BZ=zeros(NZ,NZ);
gpsi_R=zeros(size_X,size_Z);
gpsi_Z=zeros(size_X,size_Z);

for (x=3:size_X-2)
    for (z=3:size_Z-2)
        gpsi_R(x,z)=(1/12)*(-psi2D(x+2,z)+psi2D(x-2,z))+(2/3)*(psi2D(x+1,z)-psi2D(x-1,z));
        gpsi_Z(x,z)=(1/12)*(-psi2D(x,z+2)+psi2D(x,z-2))+(2/3)*(psi2D(x,z+1)-psi2D(x,z-1));
%         gpsi_Z(x,z)=0.5*(psi2D(x,z+1)-psi2D(x,z-1));
    end
end
gpsi_R=gpsi_R/DX;
gpsi_Z=gpsi_Z/DX;

BX=-gpsi_Z./Rpos_map;
BZ=gpsi_R./Rpos_map;

% [BX BZ]=gradient(psi2D,X_scale,Z_scale);
% BX=-BX./Rpos_map;
% BZ=BZ./Rpos_map;
% imagesc((BX-BR_XZ_map)',[-0.05 0.05]);

figure(4);
imagesc(abs(BX-BX_XZ_map)',[-0.05 0.05]);
axis xy;
colorbar;
% BZ=BZ.*mask_XZ_map;
% BZ=BZ+BZ_XZ_map.*(1-mask_XZ_map);
% 
% for(z=1:NZ)
%     for(x=2:NZ)
%         integ_BZ(x,z)=0.5*(BZ_XZ_map(x,z)+BZ_XZ_map(x-1,z))*((R0+(x-1)*DX-(mid_X-1)*DX)+0.5*DX);
%     end
% end
% integ_BZ=integ_BZ.*mask_XZ_map;
% psi_recalc=zeros(NZ,NZ);
% for(z=1:NZ)
%     for(x=2:NZ)
%         psi_recalc(x,z)=psi_recalc(x-1,z)+integ_BZ(x,z)*DX;
%     end
% end
% 
% figure(4);
% imagesc((psi_XZ_map-psi_recalc)',[-0.05 0.05]);
% axis xy;
% colorbar;
% % psi_recalc=psi_recalc.*mask_XZ_map;
% 
% psi2D=psi_recalc+psi2D.*(1-mask_XZ_map);
% psi2D=min(psi2D,0.1);
% 
% 
% for (x=3:NZ-2)
%     for (z=3:NZ-2)
%         gpsi_R(x,z)=(1/12)*(-psi2D(x+2,z)+psi2D(x-2,z))+(2/3)*(psi2D(x+1,z)-psi2D(x-1,z));
%         gpsi_Z(x,z)=(1/12)*(-psi2D(x,z+2)+psi2D(x,z-2))+(2/3)*(psi2D(x,z+1)-psi2D(x,z-1));
% %         gpsi_Z(x,z)=0.5*(psi2D(x,z+1)-psi2D(x,z-1));
%     end
% end
% gpsi_R=gpsi_R/DX;
% gpsi_Z=gpsi_Z/DX;
% 
% BX=-gpsi_Z./Rpos_map;
% BZ=gpsi_R./Rpos_map;