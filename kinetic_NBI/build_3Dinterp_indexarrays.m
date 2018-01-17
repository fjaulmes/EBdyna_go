function [INDEX_LIST3D_1 INDEX_LIST3D_2 INDEX_LIST3D_3 INDEX_LIST3D_4 INDEX_LIST3D_5 INDEX_LIST3D_6 INDEX_LIST3D_7 INDEX_LIST3D_8 slopex slopey slopez] = build_3Dinterp_indexarrays(X, Y, Z, DX,DY,DZ, x, y, z)

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
Zindex_precise=z;
% Zindex_precise=(z/DZ);  % no need since psi is already normalized
pindexx=floor(Xindex_precise);
pindexy=floor(Yindex_precise);
pindexz=floor(Zindex_precise);
indexx=pindexx+1; % 2pi is excluded
indexy=pindexy+1; % 2pi is excluded
indexz=pindexz+1; % size_r is excluded

slopex = (Xindex_precise - pindexx) ;
slopey = (Yindex_precise - pindexy) ;
slopez = (Zindex_precise - pindexz) ;


% 
% for n=1:length(x)
%     
%     pindexx(n) = find((x(n) >= X), 1, 'last');
%     indexx(n) = find((x(n) <= X), 1, 'first');
%     
%     if isempty(pindexx(n))
%         warning('interpolating x value before beginning');
%         pindexx(n) = indexx(n);
%         slopex(n) = 0;
%     elseif isempty(indexx(n))
%         warning('interpolating x value after end');
%         indexx(n) = pindexx(n);
%         slopex(n) = 0;
%     elseif pindexx(n) == indexx(n)
%         slopex(n) = 0
%     else
%         Xp = X(pindexx(n));
%         slopex(n) = (x(n) - Xp) / (X(indexx(n)) - Xp);
%     end
%     
%     pindexy(n) = find((y(n) >=Y), 1, 'last');
%     indexy(n) = find((y(n) <= Y), 1, 'first');
%     
%     if isempty(pindexy(n))
%         warning('interpolating y value before beginning');
%         pindexy(n) = indexy(n);
%         slopey(n) = 0;
%     elseif isempty(indexy(n))
%         warning('interpolating y value after end');
%         indexy(n) = pindexy(n);
%         slopey(n) = 0;
%     elseif pindexy(n) == indexy(n)
%         slopey(n) = 0
%     else
%         Yp = Y(pindexy(n));
%         slopey(n) = (y(n) - Yp) / (Y(indexy(n)) - Yp);
%     end
%     
%     pindexz(n) = find((z(n) >= Z), 1, 'last');
%     indexz(n) = find((z(n)  <= Z), 1, 'first');
%     
%     if isempty(pindexz(n))
%         warning('interpolating z value before beginning');
%         pindexz(n) = indexz(n);
%         slopez(n) = 0;
%     elseif isempty(indexz(n))
%         warning('interpolating z value after end');
%         indexz(n) = pindexz(n);
%         slopez(n) = 0;
%     elseif pindexz(n) == indexz(n)
%         slopez(n) = 0
%     else
%         Zp = Z(pindexz(n));
%         slopez(n) = (z(n)  - Zp) / (Z(indexz(n)) - Zp);
%     end
%     
% end

    INDEX_LIST3D_1 = sub2ind([length(X),length(Y) length(Z)] ,pindexx, pindexy, pindexz );
    INDEX_LIST3D_2 = sub2ind([length(X),length(Y) length(Z)] ,indexx, pindexy, pindexz );
    INDEX_LIST3D_3 = sub2ind([length(X),length(Y) length(Z)] ,pindexx, indexy, pindexz );
    INDEX_LIST3D_4 = sub2ind([length(X),length(Y) length(Z)] ,indexx, indexy, pindexz );
    INDEX_LIST3D_5 = sub2ind([length(X),length(Y) length(Z)] ,pindexx, pindexy, indexz );
    INDEX_LIST3D_6 = sub2ind([length(X),length(Y) length(Z)] ,indexx, pindexy, indexz );
    INDEX_LIST3D_7 = sub2ind([length(X),length(Y) length(Z)] ,pindexx, indexy, indexz );
    INDEX_LIST3D_8 = sub2ind([length(X),length(Y) length(Z)] ,indexx, indexy, indexz );

