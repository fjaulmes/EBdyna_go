    [XX_small ZZ_small]=meshgrid(scale_X,scale_Z);
    X_scale_data=reshape(XX_small,size_X*size_Z,1);
    Z_scale_data=reshape(ZZ_small,size_X*size_Z,1);
    XZ_mesh=[X_scale_data Z_scale_data];
    XZ_mesh_dtri=DelaunayTri(XZ_mesh);

[residue Xinf]=(min(abs(X_scale-scale_X(1))));
[residue Xsup]=(min(abs(X_scale-scale_X(end))));
[residue Zinf]=(min(abs(Z_scale-scale_Z(1))));
[residue Zsup]=(min(abs(Z_scale-scale_Z(end))));
    psiH_XZsmall_map=zeros(size_X,size_Z);
    psiH_XZsmall_map(:,:)=psiH_XZ_map(Xinf:Xsup,Zinf:Zsup);

    psi_XZsmall_map=zeros(size_X,size_Z);
    psi_XZsmall_map(:,:)=psi_XZ_map(Xinf:Xsup,Zinf:Zsup);

%     lap_psi_XZsmall_map=zeros(size_X,size_Z);
%     lap_psi_XZsmall_map(:,:)=lap_psi_XZ_map(Xinf:Xsup,Zinf:Zsup);
% 
%     lap_psiH_XZsmall_map=zeros(size_X,size_Z);
%     lap_psiH_XZsmall_map(:,:)=lap_psiH_XZ_map(Xinf:Xsup,Zinf:Zsup);
% 
%     lap_B_XZsmall_map=zeros(size_X,size_Z);
%     lap_B_XZsmall_map(:,:)=lap_B(Xinf:Xsup,Zinf:Zsup);


    Bphi_XZsmall_map=zeros(size_X,size_Z);
    Bphi_XZsmall_map(:,:)=Bphi_XZ_map(Xinf:Xsup,Zinf:Zsup);

    BpolX_final_XZsmall_map=zeros(size_X,size_Z);
    BpolX_final_XZsmall_map(:,:)=BX_XZ_map(Xinf:Xsup,Zinf:Zsup);

    BpolZ_final_XZsmall_map=zeros(size_X,size_Z);
    BpolZ_final_XZsmall_map(:,:)=BZ_XZ_map(Xinf:Xsup,Zinf:Zsup);

    Rpos_XZsmall_map=zeros(size_X,size_Z);
    Rpos_XZsmall_map(:,:)=Rpos_map(Xinf:Xsup,Zinf:Zsup);

    radial_XZsmall_map=zeros(size_X,size_Z);
    radial_XZsmall_map(:,:)=radial_XZ_map(Xinf:Xsup,Zinf:Zsup);
    

%     Bphi_XZsmall_map=interp2(X_scale,Z_scale,Bphi_XZ_map',scale_X,scale_Z')';
%     Rpos_XZsmall_map=interp2(X_scale,Z_scale,Rpos_map',scale_X,scale_Z')';
%     radial_XZsmall_map=interp2(X_scale,Z_scale,radial_XZ_map',scale_X,scale_Z')';
%     psi_XZsmall_map=interp2(X_scale,Z_scale,psi_XZ_map',scale_X,scale_Z')';
%     psiH_XZsmall_map=interp2(X_scale,Z_scale,psiH_XZ_map',scale_X,scale_Z')';
    
    
%     lap_psi_XZsmall_map=interp2(X_scale,Z_scale,lap_psi_XZ_map',scale_X,scale_Z')';
%     lap_psiH_XZsmall_map=interp2(X_scale,Z_scale,lap_psiH_XZ_map',scale_X,scale_Z')';
%     lap_B_XZsmall_map=interp2(X_scale,Z_scale,lap_B',scale_X,scale_Z')';
% 
%     BpolX_initial_XZsmall_map=interp2(X_scale,Z_scale,BpolX_initial_map',scale_X,scale_Z')';
%     BpolZ_initial_XZsmall_map=interp2(X_scale,Z_scale,BpolZ_initial_map',scale_X,scale_Z')';
    
size_X=length(scale_X);
size_Z=length(scale_Z);
    
    theta_XZsmall_map=zeros(size_X,size_Z);
    theta_XZsmall_map(:,:)=theta_XZ_map(Xinf:Xsup,Zinf:Zsup);

%     theta_XZsmall_map=interp2(X_scale,Z_scale,theta_XZ_map',scale_X,scale_Z')';
% 
%     theta_XZsmall_map(mid_X:size_X,mid_Z+round(Z_axis/DX))=0.25*theta_XZsmall_map(mid_X:size_X,mid_Z+round(Z_axis/DX)+1);
%     theta_XZsmall_map(mid_X:size_X,mid_Z+round(Z_axis/DX)-1)=0.75*2*pi+0.25*theta_XZsmall_map(mid_X:size_X,mid_Z+round(Z_axis/DX)-2);