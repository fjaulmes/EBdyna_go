
close all

E_data=reshape(psi_star_dot_PR_map(:,1:size_r),NP*size_r,1);
Efield_3_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,E_data(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
Efield_3_XZ_map=Efield_3_XZ_map';
Efield_3_XZ_map(isnan(Efield_3_XZ_map))=0;

Bfield_3_XZmap=Bphi_XZsmall_map./Rpos_XZsmall_map;

B_data=reshape(BstarX_PR_map(:,1:size_r),NP*size_r,1);
BstarX_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,B_data(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
BstarX_XZ_map=BstarX_XZ_map';
BstarX_XZ_map(isnan(BstarX_XZ_map))=0;

B_data=reshape(BstarZ_PR_map(:,1:size_r),NP*size_r,1);
BstarZ_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,B_data(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
BstarZ_XZ_map=BstarZ_XZ_map';
BstarZ_XZ_map(isnan(BstarZ_XZ_map))=0;

verif_dot_product=EX_XZ_map.*BstarX_XZ_map+EZ_XZ_map.*BstarZ_XZ_map+Efield_3_XZ_map.*Bfield_3_XZmap;
%verif_dot_product=abs(EX_XZ_map.*Bstar_X_XZ_map+EZ_XZ_map.*Bstar_Z_XZ_map+Efield_3_RZmap.*Bfield_3_RZmap);
%verif_dot_product=abs(ER_XZ_zoom_map.*Bpol_X_zoom_map+EZ_XZ_zoom_map.*Bpol_Z_zoom_map+Ephi_XZ_zoom_map_num.*Bphi_XZ_zoom_map);

figure(1)
set(gca,'FontSize',22);

imagesc(scale_X,scale_Z,verif_dot_product',[-1 1]*40);axis xy;
xlabel('X (m)');
ylabel('Z (m)');
xlim([-0.4 0.5]);
ylim([-0.6 0.6]);
colorbar;
axis xy square
hold on;

disp('Total cumulated errors');
disp(sum(sum(verif_dot_product)));




