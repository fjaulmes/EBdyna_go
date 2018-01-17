RR=Rpos_PR_map(:,1:pTAE_sup)-R0;

finesse_data_X=reshape((Rpos_PR_map(:,1:pTAE_sup)-R0),NP*pTAE_sup,1);
finesse_data_Z=reshape(Z_PR_map(:,1:pTAE_sup),NP*pTAE_sup,1);

% NR=size_r_TAE+8
% finesse_data_X_extended=reshape((Rpos_PR_map(:,1:NR)-R0),NP*NR,1);
% finesse_data_Z_extended=reshape(Z_PR_map(:,1:NR),NP*NR,1);

PR_DT_INTERPOLATION_METHOD='linear';    % for the cast to PR maps
% getting the reference field in PR coords
BHpol_X_PR_map_size_r_TAE(:,:)=BHpol_X_PR_map(:,1:pTAE_sup);
BHpol_Z_PR_map_size_r_TAE(:,:)=BHpol_Z_PR_map(:,1:pTAE_sup);
BHpol_X_PR_map_size_r_TAE(:,1)=BHpol_X_PR_map_size_r_TAE(:,2);
BHpol_Z_PR_map_size_r_TAE(:,1)=BHpol_Z_PR_map_size_r_TAE(:,2);



% [XX_small ZZ_small]=meshgrid(scale_X,scale_Z);
%     X_scale_data=reshape(XX_small,sizeX*sizeZ,1);
%     Z_scale_data=reshape(ZZ_small,sizeX*sizeZ,1);
%     XZ_mesh=[X_scale_data Z_scale_data];
%     XZ_mesh_dtri=DelaunayTri(XZ_mesh);
%     psiH_XZsmall_map=zeros(sizeX,sizeZ);
%     psiH_XZsmall_map(:,:)=psiH_XZ_map(Xinf:Xsup,Zinf:Zsup);

gd_data=reshape((psi_XZ_map(Xinf:Xsup,Zinf:Zsup))',(Xsup-Xinf+1)*(Zsup-Zinf+1),1);


psi_XZsmall_map=griddata(gdX_data,gdZ_data,gd_data,XX_small,ZZ_small);
psi_XZsmall_map(isnan(psi_XZsmall_map))=0;
psi_XZsmall_map=psi_XZsmall_map';

%     Bphi_XZsmall_map=zeros(sizeX,sizeZ);
%     Bphi_XZsmall_map(:,:)=Bphi_XZ_map(Xinf:Xsup,Zinf:Zsup);
gd_data=reshape((Bphi_XZ_map(Xinf:Xsup,Zinf:Zsup))',(Xsup-Xinf+1)*(Zsup-Zinf+1),1);
Bphi_XZsmall_map=griddata(gdX_data,gdZ_data,gd_data,XX_small,ZZ_small);
Bphi_XZsmall_map(isnan(Bphi_XZsmall_map))=0;
Bphi_XZsmall_map=Bphi_XZsmall_map';

%     BpolX_initial_XZsmall_map=zeros(sizeX,sizeZ);
%     BpolX_initial_XZsmall_map(:,:)=BpolX_initial_map(Xinf:Xsup,Zinf:Zsup);
gd_data=reshape((BpolX_initial_map(Xinf:Xsup,Zinf:Zsup))',(Xsup-Xinf+1)*(Zsup-Zinf+1),1);
BpolX_initial_XZsmall_map=griddata(gdX_data,gdZ_data,gd_data,XX_small,ZZ_small);
BpolX_initial_XZsmall_map(isnan(BpolX_initial_XZsmall_map))=0;
BpolX_initial_XZsmall_map=BpolX_initial_XZsmall_map';

%     BpolZ_initial_XZsmall_map=zeros(sizeX,sizeZ);
%     BpolZ_initial_XZsmall_map(:,:)=BpolZ_initial_map(Xinf:Xsup,Zinf:Zsup);
gd_data=reshape((BpolZ_initial_map(Xinf:Xsup,Zinf:Zsup))',(Xsup-Xinf+1)*(Zsup-Zinf+1),1);
BpolZ_initial_XZsmall_map=griddata(gdX_data,gdZ_data,gd_data,XX_small,ZZ_small);
BpolZ_initial_XZsmall_map(isnan(BpolZ_initial_XZsmall_map))=0;
BpolZ_initial_XZsmall_map=BpolZ_initial_XZsmall_map';

Btot_XZ_small_map=sqrt(BpolX_initial_XZsmall_map.^2+BpolZ_initial_XZsmall_map.^2+Bphi_XZsmall_map.^2);

bX_XZ_small_map=BpolX_initial_XZsmall_map./Btot_XZ_small_map;
bZ_XZ_small_map=BpolZ_initial_XZsmall_map./Btot_XZ_small_map;

%     Rpos_XZsmall_map=zeros(sizeX,sizeZ);
%     Rpos_XZsmall_map(:,:)=Rpos_map(Xinf:Xsup,Zinf:Zsup);
gd_data=reshape((Rpos_map(Xinf:Xsup,Zinf:Zsup))',(Xsup-Xinf+1)*(Zsup-Zinf+1),1);
Rpos_XZsmall_map=griddata(gdX_data,gdZ_data,gd_data,XX_small,ZZ_small);
Rpos_XZsmall_map(isnan(Rpos_XZsmall_map))=0;
Rpos_XZsmall_map=Rpos_XZsmall_map';



%     radial_XZsmall_map=zeros(sizeX,sizeZ);
%     radial_XZsmall_map(:,:)=radial_XZ_map(Xinf:Xsup,Zinf:Zsup);
 

vA_data=reshape((vA_PR_map(:,1:pTAE_sup)),NP*pTAE_sup,1);
vA_XZ_map=griddata(finesse_data_X,finesse_data_Z,vA_data,XX_small,ZZ_small);
vA_XZ_map(isnan(vA_XZ_map)) = 0; 
vA_XZ_map=vA_XZ_map';
vA_max=max(max(vA_XZ_map));
vA_inf=1*1e6;
vA_XZ_map(vA_XZ_map==0)=vA_max;
vA_XZ_map(vA_XZ_map<=vA_inf)=vA_inf;
vA_XZ_map=vA_XZ_map*1e-6;
