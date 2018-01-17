    [XX_small ZZ_small]=meshgrid(scale_X,scale_Z);
    X_scale_data=reshape(XX_small,sizeX*sizeZ,1);
    Z_scale_data=reshape(ZZ_small,sizeX*sizeZ,1);
    XZ_mesh=[X_scale_data Z_scale_data];
    XZ_mesh_dtri=DelaunayTri(XZ_mesh);

[residue Xinf]=(min(abs(X_scale-scale_X(1))));
[residue Xsup]=(min(abs(X_scale-scale_X(end))));
[residue Zinf]=(min(abs(Z_scale-scale_Z(1))));
[residue Zsup]=(min(abs(Z_scale-scale_Z(end))));
%     psiH_XZsmall_map=zeros(sizeX,sizeZ);
%     psiH_XZsmall_map(:,:)=psiH_XZ_map(Xinf:Xsup,Zinf:Zsup);

%     psi_norm_XZsmall_map=zeros(sizeX,sizeZ);
%     psi_norm_XZsmall_map(:,:)=radial_XZ_map(Xinf:Xsup,Zinf:Zsup);

    psi_XZsmall_map=zeros(sizeX,sizeZ);
    psi_XZsmall_map(:,:)=psi_XZ_map(Xinf:Xsup,Zinf:Zsup);
    psi_norm_XZsmall_map=interp1(psi_scale,1:Nradial,psi_XZsmall_map);

    Bphi_XZsmall_map=zeros(sizeX,sizeZ);
    Bphi_XZsmall_map(:,:)=Bphi_XZ_map(Xinf:Xsup,Zinf:Zsup);

    BpolX_initial_XZsmall_map=zeros(sizeX,sizeZ);
    BpolX_initial_XZsmall_map(:,:)=BpolX_initial_map(Xinf:Xsup,Zinf:Zsup);

    BpolZ_initial_XZsmall_map=zeros(sizeX,sizeZ);
    BpolZ_initial_XZsmall_map(:,:)=BpolZ_initial_map(Xinf:Xsup,Zinf:Zsup);

    Rpos_XZsmall_map=zeros(sizeX,sizeZ);
    Rpos_XZsmall_map(:,:)=Rpos_map(Xinf:Xsup,Zinf:Zsup);
    
BHpolX_initial_XZsmall_map=zeros(sizeX,sizeZ);
BHpolX_initial_XZsmall_map(:,:)=BHpol_X_map(Xinf:Xsup,Zinf:Zsup);

BHpolZ_initial_XZsmall_map=zeros(sizeX,sizeZ);
BHpolZ_initial_XZsmall_map(:,:)=BHpol_Z_map(Xinf:Xsup,Zinf:Zsup);
%     radial_XZsmall_map=zeros(sizeX,sizeZ);
%     radial_XZsmall_map(:,:)=radial_XZ_map(Xinf:Xsup,Zinf:Zsup);
 
 
DT_INTERPOLATION_METHOD='quadratic'     % by default
 
% psi_max_value=interp1(1:Nradial,psi_scale,psi_max+2)
psi_max_value=interp1(1:Nradial,psi_scale,psi_rank_q1+1)

psi_data=reshape(psi_PR_map(:,1:NR),NB_THETA*NR,1);
% psi_XZ_map=griddata(finesse_data_X_extended,finesse_data_Z_extended,psi_data,XX_small,ZZ_small,'cubic');
psi_XZ_map_ini=dtinterp(finesse_mesh_extended,finesse_mesh_extended_dtri,psi_PR_map(IDTRI_EXT),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
psi_XZ_map_ini(isnan(psi_XZ_map_ini))=0;
psi_XZ_map_ini=psi_XZ_map_ini';
psi_XZ_map_mask=psi_XZ_map_ini*0;

for (x=3:sizeX-2)
    for (z=3:sizeZ-2)
        if psi_max_value<0
            if psi_XZ_map_ini(x,z)<psi_max_value
                psi_XZ_map_mask(x,z)=1;
            end
        else
            if psi_XZ_map_ini(x,z)>psi_max_value
                psi_XZ_map_mask(x,z)=1;
            end
        end
    end
end
