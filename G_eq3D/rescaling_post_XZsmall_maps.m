    [XX_small ZZ_small]=meshgrid(scale_X,scale_Z);
    X_scale_data=reshape(XX_small,size_X*size_Z,1);
    Z_scale_data=reshape(ZZ_small,size_X*size_Z,1);
    XZ_mesh=[X_scale_data Z_scale_data];
    XZ_mesh_dtri=DelaunayTri(XZ_mesh);

[residue Xinf]=(min(abs(X_scale-scale_X(1))));
[residue Xsup]=(min(abs(X_scale-scale_X(end))));
[residue Zinf]=(min(abs(Z_scale-scale_Z(1))));
[residue Zsup]=(min(abs(Z_scale-scale_Z(end))));

size_X=length(scale_X);
size_Z=length(scale_Z);

    Bphi_XZsmall_map=zeros(size_X,size_Z);
    Bphi_XZsmall_map(:,:)=Bphi_XZ_map(Xinf:Xsup,Zinf:Zsup);

    psiH_XZsmall_map=zeros(size_X,size_Z);
    psiH_XZsmall_map(:,:)=psiH_XZ_map(Xinf:Xsup,Zinf:Zsup);

    Bphi_XZsmall_map=zeros(size_X,size_Z);
    Bphi_XZsmall_map(:,:)=Bphi_XZ_map(Xinf:Xsup,Zinf:Zsup);

    psi_XZsmall_map=zeros(size_X,size_Z);
    psi_XZsmall_map(:,:)=psi_XZ_map(Xinf:Xsup,Zinf:Zsup);

    BpolX_final_XZsmall_map=zeros(size_X,size_Z);
    BpolX_final_XZsmall_map(:,:)=BX_XZ_map(Xinf:Xsup,Zinf:Zsup);

    BpolZ_final_XZsmall_map=zeros(size_X,size_Z);
    BpolZ_final_XZsmall_map(:,:)=BZ_XZ_map(Xinf:Xsup,Zinf:Zsup);

    Rpos_XZsmall_map=zeros(size_X,size_Z);
    Rpos_XZsmall_map(:,:)=Rpos_map(Xinf:Xsup,Zinf:Zsup);

    radial_XZsmall_map=zeros(size_X,size_Z);
    radial_XZsmall_map(:,:)=radial_XZ_map(Xinf:Xsup,Zinf:Zsup);

    

    
    theta_XZsmall_map=zeros(size_X,size_Z);
    theta_XZsmall_map(:,:)=theta_XZ_map(Xinf:Xsup,Zinf:Zsup);
