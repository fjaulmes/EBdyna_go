function [INDEX_LIST3D_1 INDEX_LIST3D_2 INDEX_LIST3D_3 INDEX_LIST3D_4 INDEX_LIST3D_5 INDEX_LIST3D_6 INDEX_LIST3D_7 INDEX_LIST3D_8 slopex slopey slopez] = build_3Dinterp_indexarrays(scX, scY, scZ, DX,DY,DZ, x, y, z)

% should use identical increments for each scale

pindexx=zeros(length(x),1);
indexx=zeros(length(x),1);
slopex=zeros(length(x),1);

pindexy=zeros(length(x),1);
indexy=zeros(length(x),1);
slopey=zeros(length(x),1);

pindexz=zeros(length(x),1);
indexz=zeros(length(x),1);
slopez=zeros(length(x),1);

INDEX_LIST3D_1=zeros(length(x),1);
INDEX_LIST3D_2=zeros(length(x),1);
INDEX_LIST3D_3=zeros(length(x),1);
INDEX_LIST3D_4=zeros(length(x),1);
INDEX_LIST3D_5=zeros(length(x),1);
INDEX_LIST3D_6=zeros(length(x),1);
INDEX_LIST3D_7=zeros(length(x),1);
INDEX_LIST3D_8=zeros(length(x),1);

Xindex_precise=(x/DX)+1; % phi
Yindex_precise=(y/DY)+1; % theta
% z=min(z,scZ(end)-1); % psi (excluding separatrix)
Zindex_precise=((z-1)/DZ)+1;
% Zindex_precise=min(Zindex_precise,length(Z));

lnX=length(scX);
lnY=length(scY);
lnZ=length(scZ);

pindexx=floor(Xindex_precise);
pindexx=min(pindexx,lnX-1); % 2pi is excluded
pindexy=floor(Yindex_precise);
pindexy=min(pindexy,lnY-1); % 2pi is excluded
pindexz=floor(Zindex_precise);
pindexz=min(pindexz,lnZ-1); % max radial index is excluded
indexx=pindexx+1; 
indexy=pindexy+1; 
indexz=pindexz+1; 

slopex = (Xindex_precise - pindexx) ;
slopey = (Yindex_precise - pindexy) ;
slopez = (Zindex_precise - pindexz) ;



INDEX_LIST3D_1 = sub2ind([lnX, lnY,  lnZ] ,pindexx, pindexy, pindexz );
INDEX_LIST3D_2 = sub2ind([lnX, lnY,  lnZ] ,indexx, pindexy, pindexz );
INDEX_LIST3D_3 = sub2ind([lnX, lnY,  lnZ] ,pindexx, indexy, pindexz );
INDEX_LIST3D_4 = sub2ind([lnX, lnY,  lnZ] ,indexx, indexy, pindexz );
INDEX_LIST3D_5 = sub2ind([lnX, lnY,  lnZ] ,pindexx, pindexy, indexz );
INDEX_LIST3D_6 = sub2ind([lnX, lnY,  lnZ] ,indexx, pindexy, indexz );
INDEX_LIST3D_7 = sub2ind([lnX, lnY,  lnZ] ,pindexx, indexy, indexz );
INDEX_LIST3D_8 = sub2ind([lnX, lnY,  lnZ] ,indexx, indexy, indexz );

