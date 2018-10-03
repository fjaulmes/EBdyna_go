[XX_small, ZZ_small]=meshgrid(scale_X,scale_Z);
X_scale_data=reshape(XX_small,size_X*size_Z,1);
Z_scale_data=reshape(ZZ_small,size_X*size_Z,1);
XZ_mesh=[X_scale_data Z_scale_data];
XZ_mesh_dtri=DelaunayTri(XZ_mesh);

[~, Xinf]=(min(abs(X_scale-scale_X(1))));
[~, Xsup]=(min(abs(X_scale-scale_X(end))));
[~, Zinf]=(min(abs(Z_scale-scale_Z(1))));
[~, Zsup]=(min(abs(Z_scale-scale_Z(end))));
psiH_XZsmall_map=zeros(size_X,size_Z);

psi_XZsmall_map=zeros(size_X,size_Z);

Bphi_XZsmall_map=zeros(size_X,size_Z);
BpolX_initial_XZsmall_map=zeros(size_X,size_Z);
BpolZ_initial_XZsmall_map=zeros(size_X,size_Z);

Rpos_XZsmall_map=zeros(size_X,size_Z);
radial_XZsmall_map=zeros(size_X,size_Z);
q_initial_XZsmall_map=zeros(size_X,size_Z);

q_initial_XZsmall_map=interp2(X_scale,Z_scale,q_XZ_map',scale_X,scale_Z')';

Bphi_XZsmall_map=interp2(X_scale,Z_scale,Bphi_XZ_map',scale_X,scale_Z')';
Rpos_XZsmall_map=interp2(X_scale,Z_scale,Rpos_map',scale_X,scale_Z')';
radial_XZsmall_map=interp2(X_scale,Z_scale,radial_XZ_map',scale_X,scale_Z')';
psi_XZsmall_map=interp2(X_scale,Z_scale,psi_XZ_map',scale_X,scale_Z')';
psiH_XZsmall_map=interp2(X_scale,Z_scale,psiH_XZ_map',scale_X,scale_Z')';

BpolX_initial_XZsmall_map=interp2(X_scale,Z_scale,BpolX_initial_map',scale_X,scale_Z')';
BpolZ_initial_XZsmall_map=interp2(X_scale,Z_scale,BpolZ_initial_map',scale_X,scale_Z')';

size_X=length(scale_X);
size_Z=length(scale_Z);

theta_XZsmall_map=zeros(size_X,size_Z);
theta_XZsmall_map=interp2(X_scale,Z_scale,theta_XZ_map',scale_X,scale_Z')';


psi_XZsmall_map(isnan(psi_XZsmall_map))=0;
psiH_XZsmall_map(isnan(psiH_XZsmall_map))=0;