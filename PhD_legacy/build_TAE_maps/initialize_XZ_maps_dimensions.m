ZOOM_FACTOR=1
DX=DX/ZOOM_FACTOR

%     NZ=NX;

Z_axis_shift=round(Z_axis/DX);

%     rx_precise=rx_evol(end);
%     rx_rank=round(rx_precise/DX)+1;
    minX=max(mid_X-round(4.6*pTAE_sup),3);
    maxX=min(mid_X+round(4.6*pTAE_sup),NZ-2);
    minZ=max(mid_Z+Z_axis_shift-round(4.6*elongation*pTAE_sup),3);
    maxZ=min(mid_Z+Z_axis_shift+round(4.6*elongation*pTAE_sup),NZ-2);
    
    sizeX=ZOOM_FACTOR*2*ceil(0.5*(maxX-minX));
    sizeZ=ZOOM_FACTOR*2*ceil(0.5*(maxZ-minZ));
    
    mid_X=ZOOM_FACTOR*ceil(0.5*(maxX-minX))+1;
    mid_Z=ZOOM_FACTOR*ceil(0.5*(maxZ-minZ))+1;
    scale_X=DX*((1:sizeX)-mid_X);
    scale_Z=DX*((1:sizeZ)-(mid_Z-Z_axis_shift));
    
%     scale_X_data=reshape(XX(:,:),NZ*NZ,1);
%     scale_Z_data=reshape(ZZ(:,:),NZ*NZ,1);




[XX_small ZZ_small]=meshgrid(scale_X,scale_Z);
X_scale_data=reshape(XX_small,sizeX*sizeZ,1);
Z_scale_data=reshape(ZZ_small,sizeX*sizeZ,1);
% XZ_mesh=[X_scale_data Z_scale_data];
% XZ_mesh_dtri=DelaunayTri(XZ_mesh);


    
    
[residue Xinf]=(min(abs(X_scale-scale_X(1))));
[residue Xsup]=(min(abs(X_scale-scale_X(end))));
[residue Zinf]=(min(abs(Z_scale-scale_Z(1))));
[residue Zsup]=(min(abs(Z_scale-scale_Z(end))));

    Xinf=max(Xinf-1,3);
    Xsup=min(Xsup+1,NZ-2);
    Zinf=max(Zinf-1,3);
    Zsup=min(Zsup+1,NZ-2);

[XX ZZ]=meshgrid(X_scale(Xinf:Xsup),Z_scale(Zinf:Zsup));
gdX_data=reshape(XX,(Xsup-Xinf+1)*(Zsup-Zinf+1),1);
gdZ_data=reshape(ZZ,(Xsup-Xinf+1)*(Zsup-Zinf+1),1);

% 
% filename=strcat(DATA_FOLDER,'Epot_psi_star_dot_evol.mat');
% load(filename);


