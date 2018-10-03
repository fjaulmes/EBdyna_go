
%     NZ=NX;

Z_axis_shift=round(Z_axis/DX)

    rx_precise=rx_evol(end);
    rx_rank=round(rx_precise/DX)+1;
    minX=max(mid_X-round(1.1*rx_rank),3);
    maxX=min(mid_X+round(1.4*rx_rank),NZ-2);
    minZ=max(mid_Z+Z_axis_shift-round(0.9*elongation*rx_rank),3);
    maxZ=min(mid_Z+Z_axis_shift+round(0.9*elongation*rx_rank),NZ-2);
    
    sizeX=2*ceil(0.5*(maxX-minX));
    sizeZ=2*ceil(0.5*(maxZ-minZ));
    
    mid_X=ceil(0.5*(maxX-minX))+1;
    mid_Z=ceil(0.5*(maxZ-minZ))+1;
    scale_X=DX*((1:sizeX)-mid_X);
    scale_Z=DX*((1:sizeZ)-(mid_Z-Z_axis_shift));
    
%     scale_X_data=reshape(XX(:,:),NZ*NZ,1);
%     scale_Z_data=reshape(ZZ(:,:),NZ*NZ,1);

    NB_PHI=129;
    NB_PHI_DATA_HALF=round(0.5*NB_PHI-1);
    DPHI=(2*pi)/(NB_PHI-1);