    [XX_small ZZ_small]=meshgrid(scale_X,scale_Z);
    X_scale_data=reshape(XX_small,size_X*size_Z,1);
    Z_scale_data=reshape(ZZ_small,size_X*size_Z,1);
    XZ_mesh=[X_scale_data Z_scale_data];
    XZ_mesh_dtri=DelaunayTri(XZ_mesh);
    
    Bphi_XZsmall_map=interp2(X_scale,Z_scale,Bphi_XZ_map',scale_X,scale_Z')';
    Rpos_XZsmall_map=interp2(X_scale,Z_scale,Rpos_map',scale_X,scale_Z')';
    radial_XZsmall_map=interp2(X_scale,Z_scale,radial_XZ_map',scale_X,scale_Z')';
    theta_XZsmall_map=interp2(X_scale,Z_scale,theta_XZ_map',scale_X,scale_Z')';
    psi_XZsmall_map=interp2(X_scale,Z_scale,psi_XZ_map',scale_X,scale_Z')';

    BpolX_initial_XZsmall_map=interp2(X_scale,Z_scale,BpolX_initial_map',scale_X,scale_Z')';
    BpolZ_initial_XZsmall_map=interp2(X_scale,Z_scale,BpolZ_initial_map',scale_X,scale_Z')';

    BHpolX_initial_XZsmall_map=interp2(X_scale,Z_scale,BHpol_X_map',scale_X,scale_Z')';
    BHpolZ_initial_XZsmall_map=interp2(X_scale,Z_scale,BHpol_Z_map',scale_X,scale_Z')';

size_X=size(scale_X,2);
size_Z=size(scale_Z,2);
