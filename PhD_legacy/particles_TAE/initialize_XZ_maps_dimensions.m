

    rx_precise=interp1(1:257,radial_r_value_flux, xPsih_zero);
    rx_rank=round(rx_precise/DX)+1;
    minX=max(mid_X-round(1.1*rx_rank),3);
    maxX=min(mid_X+round(1.1*rx_rank),NZ-2);
    minZ=max(mid_Z-round(0.9*elongation*rx_rank),3);
    maxZ=min(mid_Z+round(0.9*elongation*rx_rank),NZ-2);
    
    sizeX=2*ceil(0.5*(maxX-minX));
    sizeZ=2*ceil(0.5*(maxZ-minZ));
    
    mid_X=ceil(0.5*(maxX-minX))+1;
    mid_Z=ceil(0.5*(maxZ-minZ))+1;
    scale_X=DX*((1:sizeX)-mid_X)
    scale_Z=DX*((1:sizeZ)-mid_Z);
    
%     scale_X_data=reshape(XX(:,:),NZ*NZ,1);
%     scale_Z_data=reshape(ZZ(:,:),NZ*NZ,1);
