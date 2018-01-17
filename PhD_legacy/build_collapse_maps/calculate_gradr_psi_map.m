
psi_data=reshape(psiH_PR_map(:,1:size_r)+psi_star_PR_map(:,1:size_r),NP*size_r,1);
psi_map_tiny=dtinterp(finesse_mesh,finesse_mesh_dtri,psi_data(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
psi_map_tiny=psi_map_tiny';
psi_map_tiny(isnan(psi_map_tiny))=0;

Z_data=reshape(Z_PR_map(:,1:size_r),NP*size_r,1);
Z_map=dtinterp(finesse_mesh,finesse_mesh_dtri,Z_data(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
Z_map=Z_map';
Z_map(isnan(Z_map))=0;



GFAC=0.01;

gradr_psi_map=psi_map_tiny*0;
dr_size_map=2*GFAC*sqrt(drX_XZ_map.^2+drZ_XZ_map.^2);
x_plus_map=(Rpos_XZsmall_map-R0)+GFAC*drX_XZ_map;
z_plus_map=(Z_map)+GFAC*drZ_XZ_map;
psi_plus_map=interp2(scale_X,scale_Z,psi_map_tiny',x_plus_map,z_plus_map);
x_minus_map=(Rpos_XZsmall_map-R0)-GFAC*drX_XZ_map;
z_minus_map=(Z_map)-GFAC*drZ_XZ_map;
psi_minus_map=interp2(scale_X,scale_Z,psi_map_tiny',x_minus_map,z_minus_map);
gradr_psi_map=(psi_plus_map-psi_minus_map)./dr_size_map;
gradr_psi_map(isnan(gradr_psi_map))=0;

% for (x=3:size_X-2)
% for (z=3:size_Z-2)
%     drX=drX_XZ_map(x,z);
%     drZ=drZ_XZ_map(x,z);
%     dr_size=2*sqrt(drX^2+drZ^2);
%     x_plus=scale_X(x)+drX;
%     z_plus=scale_Z(x)+drZ;
%     psi_plus=interp2(scale_X,scale_Z,psi_map',x_plus,z_plus);
%     x_minus=scale_X(x)-drX;
%     z_minus=scale_Z(x)-drZ;
%     psi_minus=interp2(scale_X,scale_Z,psi_map',x_minus,z_minus);
%     gradr_psi_map(x,z)=(psi_plus-psi_minus)/dr_size;
% end
% end



% figure(1)
% set(gca,'FontSize',22);
% hold on
% imagesc(scale_X,scale_Z , gradr_psi_map',[-2 2]);
% axis xy
% xlabel('X')
% ylabel('Z')
% xlim([-0.3 0.4])
% ylim([-0.5 0.5])
% t_h=title('$$\nabla_r\psi (\varphi=0)$$','interpreter','latex');
% set(t_h,'Interpreter','latex');
% colorbar;
