
phi_rank



BsX_PR_map=squeeze(BsX_map_phi(phi_rank,:,:));
BsZ_PR_map=squeeze(BsZ_map_phi(phi_rank,:,:));
Bsphi_PR_map=squeeze(Bsphi_map_phi(phi_rank,:,:));

EX_PR_map=squeeze(Efield_X_map_phi(phi_rank,:,:));
EZ_PR_map=squeeze(Efield_Z_map_phi(phi_rank,:,:));





Bdata=reshape((BsX_PR_map(:,1:size_r_TAE)),NP*size_r_TAE,1);
BsX_XZ_map=griddata(finesse_data_X,finesse_data_Z,Bdata,XX_small,ZZ_small,'cubic');
BsX_XZ_map(isnan(BsX_XZ_map)) = 0;
BsX_XZ_map=BsX_XZ_map';

Bdata=reshape((BsZ_PR_map(:,1:size_r_TAE)),NP*size_r_TAE,1);
BsZ_XZ_map=griddata(finesse_data_X,finesse_data_Z,Bdata,XX_small,ZZ_small,'cubic');
BsZ_XZ_map(isnan(BsZ_XZ_map)) = 0;
BsZ_XZ_map=BsZ_XZ_map';

Bdata=reshape((Bsphi_PR_map(:,1:size_r_TAE)),NP*size_r_TAE,1);
Bsphi_XZ_map=griddata(finesse_data_X,finesse_data_Z,Bdata,XX_small,ZZ_small,'cubic');
Bsphi_XZ_map(isnan(Bsphi_XZ_map)) = 0;
Bsphi_XZ_map=Bsphi_XZ_map';


Edata=reshape((EX_PR_map(:,1:size_r_TAE)),NP*size_r_TAE,1);
EX_XZ_map=griddata(finesse_data_X,finesse_data_Z,Edata,XX_small,ZZ_small,'cubic');
EX_XZ_map(isnan(EX_XZ_map)) = 0;
EX_XZ_map=EX_XZ_map';

Edata=reshape((EZ_PR_map(:,1:size_r_TAE)),NP*size_r_TAE,1);
EZ_XZ_map=griddata(finesse_data_X,finesse_data_Z,Edata,XX_small,ZZ_small,'cubic');
EZ_XZ_map(isnan(EZ_XZ_map)) = 0;
EZ_XZ_map=EZ_XZ_map';


BX_tot_XZ_map=BpolX_initial_XZsmall_map+BsX_XZ_map;
BZ_tot_XZ_map=BpolZ_initial_XZsmall_map+BsZ_XZ_map;
Bphi_tot_XZ_map=Bphi_XZsmall_map+Bsphi_XZ_map;

Btot_XZ_map_rank=sqrt(BX_tot_XZ_map.^2+BZ_tot_XZ_map.^2+Bphi_tot_XZ_map.^2);




Ephi_XZ_map=-(EX_XZ_map.*BX_tot_XZ_map+EZ_XZ_map.*BZ_tot_XZ_map)./(Bphi_tot_XZ_map);
Etot_XZ_map_rank=sqrt(EX_XZ_map.^2+EZ_XZ_map.^2+Ephi_XZ_map.^2);

%estimate the displacement
ksi_dot_XZ_map_frame_rank=Etot_XZ_map_rank./Btot_XZ_map_rank;

exb_X=EZ_XZ_map.*Bphi_tot_XZ_map-Ephi_XZ_map.*BZ_tot_XZ_map;
exb_Z=-EX_XZ_map.*Bphi_tot_XZ_map+Ephi_XZ_map.*BX_tot_XZ_map;
exb_phi=EX_XZ_map.*BZ_tot_XZ_map-EZ_XZ_map.*BX_tot_XZ_map;
norm_exb=sqrt(exb_X.^2+exb_Z.^2+exb_phi.^2);
% exb_X=exb_X./norm_exb;
% exb_Z=exb_Z./norm_exb;
% exb_phi=exb_phi./norm_exb;
exb_X(isnan(exb_X))=0;
exb_Z(isnan(exb_Z))=0;
exb_phi(isnan(exb_phi))=0;
