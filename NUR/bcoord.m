function [R,Z,phi,dl,X,Y]=bcoord(bmwid,focl,div,nb,app,R_t,d_s,Z_s,Z_a,Rmin,Rmax,nr)

% %parameters
% nb=10000;                        %nb numbear of generated beams
% %R_t=sqrt(3)/2;
% R_t=0.25;
% d_s=4;
% Rmin=0.2;
% Rmax=1;
% nr=1000;
% Z_s=0.2;
% Z_a=0.1;
% 
% div.r = 0.00823688265;                        % ion source horizontal divergence (radians)
% div.z = 0.00823688265;                        % ion source horizontal divergence (radians)
% bmwid.r = 0.1;                      % ion source half-width (m)
% bmwid.z = 0.3;                      % ion source half-width (m)
% bmwid.opt = 1;                      % ion source half-width (m)
% focl.r = 4;                     % ion source horizontal focal length (m)
% focl.z = 4;                     % ion source horizontal focal length (m)
% app.l(1)=3;                 % app aperture to source distance
% app.l(2)=3.2;                 % app aperture to source distance
% app.r(1)=0.07;                  % app(2) apperture radius
% app.r(2)=0.05;                  % app(2) apperture radius
% app.z(1)=100;                  % app(2) apperture radius
% app.z(2)=0.1;                  % app(2) apperture radius
% app.opt(1)=0;                   % option: 0 -circ, 1 - rect
% app.opt(2)=1;                   % option: 0 -circ, 1 - rect
% 


%% the wrong code
% [dR,dZ]=blines(bmwidra,foclra,divra,nb,app);
% 
% d_sa=app(1);
% 
% % x=d_s/d_sa*(dR(2,:)-dR(1,:));
% % y=hypot(x,d_s);
% % Rt=d_s*((-x+R_t)./y);
% % lt=y+x.*Rt/R_t;
% 
% sinA=(Z_s-Z_a)/d_sa;
% cosA=sqrt(1-sinA^2);
% 
% x=(d_s*cosA-dZ(1,:)*sinA)./(d_sa*cosA+(dZ(2,:)-dZ(1,:))*sinA).*(dR(2,:)-dR(1,:));
% y=hypot(x,d_s*cosA-dZ(1,:)*sinA);
% Rt=(R_t-dR(1,:)-x).*(d_s*cosA-dZ(1,:)*sinA)./y;
% lt=y+x.*Rt./(d_s*cosA-dZ(1,:)*sinA);
% 
% lst=real(lt-sqrt(-Rt.^2+Rmax^2));
% 
% lfin_max=real(lt+sqrt(-Rt.^2+Rmax^2));
% lfin_min=real(lt-sqrt(-Rt.^2+Rmin^2));
% 
% ind1=ones(size(Rt));
% ind1(Rt<Rmin)=0;
% 
% lfin=ind1.*lfin_max+(ones(size(ind1))-ind1).*lfin_min;
% 
% lr=(ones(nr,1)*lst)+linspace(0,1,nr)'*(lfin-lst);
% 
% R=hypot(ones(nr,1)*lt-lr,ones(nr,1)*Rt);
% 
% dlr=(lfin-lst)*(1/(nr-1));
% phi=-sign(y-lr).*acos(sqrt(((R_t-dR(1,:)-x).^2+R.^2-(y-lr).^2)./(2.*(R_t-dR(1,:)-x).*R)))/pi*180;
% 
% cosR=(d_s*cosA-dZ(1,:)*sinA)./y;
% %% 
% Zh=Z_s-Z_a+(dZ(1,:)-dZ(2,:))*cosA;
% Xh=d_sa*cosA+(dZ(2,:)-dZ(1,:))*sinA;
% tgB=Zh./Xh;
% 
% dl=sqrt((dlr.*cosR).^2+(dlr.*sqrt(1-cosR.^2)).^2+(dlr.*cosR.*tgB).^2);
% 
% Z=Z_s+dZ(1,:)*cosA-lr.*cosR.*tgB;

%% The right code
[dR,dZ]=blines(bmwid,focl,div,nb,app);

d_sa=app.l(1);

sinA=(Z_s-Z_a)/d_sa;
cosA=sqrt(1-sinA^2);

X_s=-d_s*cosA;
X_a=(-d_s+d_sa)*cosA;
Y_a=R_t;
Y_s=R_t;

Xb(1,:)=X_s+dZ(1,:)*sinA;
Xb(2,:)=X_a+dZ(2,:)*sinA;
Yb(1,:)=Y_s+dR(1,:);
Yb(2,:)=Y_a+dR(2,:);
Zb(1,:)=Z_s+dZ(1,:).*cosA;
Zb(2,:)=Z_a+dZ(2,:).*cosA;

a=diff(Xb);
b=diff(Yb);
c=diff(Zb);

mod=sqrt(a.^2+b.^2+c.^2);
a=a./mod;
b=b./mod;
c=c./mod;

ba=b./a;

A=1+ba.^2;
B=2.*ba.*Yb(1,:)-2.*ba.^2.*Xb(1,:);
C=Yb(1,:).^2-2.*ba.*Yb(1,:).*Xb(1,:)+ba.^2.*Xb(1,:).^2-Rmax^2;
Cmin=Yb(1,:).^2-2.*ba.*Yb(1,:).*Xb(1,:)+ba.^2.*Xb(1,:).^2-Rmin^2;

D=B.^2-4.*A.*C;
Dmin=B.^2-4.*A.*Cmin;

if (D(:)<0)
    error('beamray did not cross the plasma')
end


ind1=ones(size(Dmin));%indexes of beam paths that strikes inner coloumn 
ind1(Dmin>0)=0;

X1=(-B-sqrt(D))./2./A;
X2max=real((-B+sqrt(D))./2./A);
X2min=real((-B-sqrt(Dmin))./2./A);
X2=ind1.*X2max+(ones(size(ind1))-ind1).*X2min;


X1(D<0)=Xb(D<0);
X2(D<0)=Xb(D<0)+1;

Y1=Yb(1,:)+ba.*(X1-Xb(1,:));
Y2=Yb(1,:)+ba.*(X2-Xb(1,:));

Z1=Zb(1,:)+c./a.*(X1-Xb(1,:));
Z2=Zb(1,:)+c./a.*(X2-Xb(1,:));

X=(ones(nr,1)*X1)+linspace(0,1,nr)'*(X2-X1);

Y=(ones(nr,1)*Y1)+linspace(0,1,nr)'*(Y2-Y1);

Z=(ones(nr,1)*Z1)+linspace(0,1,nr)'*(Z2-Z1);

R=hypot(X,Y);
phi=asin(X./R);

dl=sqrt((X2-X1).^2+(Y2-Y1).^2+(Z2-Z1).^2)./(nr-1);

end