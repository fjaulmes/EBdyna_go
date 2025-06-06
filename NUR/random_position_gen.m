function [X_v0,Y_v0,Z_v0,vx,vy,vz,REMOVED_LINES]=random_position_gen(P,X,Y,Z,E)

%E=E/m
mu =    1.66053906660 * 1e-27;
mD =    2.01410177811 * mu;
eV =    1.602176634 * 1e-19;

v_norm=sqrt(2*E*1.60217662e-19/1.660539e-27);

ax=(X(2,:)-X(1,:));
ay=(Y(2,:)-Y(1,:));
az=(Z(2,:)-Z(1,:));

mod=sqrt(ax.^2+ay.^2+az.^2);

ax=ax./mod;
ay=ay./mod;
az=az./mod;

r=rand(1,size(P,2));
P0=P-ones(size(P,1),1)*r;
Pm=P0(1:end-2,:).*P0(2:end-1,:);
Pm(Pm==0)=-1;
Pm(Pm>0)=0;
Pm(Pm<0)=1;

Pm=[zeros(1,size(Pm,2)); Pm; zeros(1,size(Pm,2))];
Pm0=[ Pm(2:end,:); zeros(1,size(Pm,2))];

X_v1=Pm.*X;
Y_v1=Pm.*Y;
Z_v1=Pm.*Z;

X_v0=Pm0.*X;
Y_v0=Pm0.*Y;
Z_v0=Pm0.*Z;

Pm2P0=P0.*Pm;

k=P0.*Pm0./(P0.*Pm0-[Pm2P0(2:end,:);zeros(1,size(Pm,2))]);
k(isnan(k))=0;
X_v0=X_v0+k.*([X_v1(2:end,:);zeros(1,size(Pm,2))]-X_v0);
Y_v0=Y_v0+k.*([Y_v1(2:end,:);zeros(1,size(Pm,2))]-Y_v0);
Z_v0=Z_v0+k.*([Z_v1(2:end,:);zeros(1,size(Pm,2))]-Z_v0);

vx=ax.*v_norm;
vy=ay.*v_norm;
vz=az.*v_norm;

Pm1d=sum(Pm);
Pm1d(Pm1d==0)=NaN;

vx=vx.*Pm1d;
vy=vy.*Pm1d;
vz=vz.*Pm1d;

vx=rmmissing(vx);
vy=rmmissing(vy);
vz=rmmissing(vz);

X_v0_sum=sum(X_v0,1);
REMOVED_LINES=logical(X_v0_sum==0);

X_v0=nonzeros(X_v0)';
Y_v0=nonzeros(Y_v0)';
Z_v0=nonzeros(Z_v0)';


if (size(vx,2)~=size(X_v0,2))
    fprintf("size of vx and X_v0 is different");
    fprintf("size vx %i",size(vx))
    fprintf("size X_v0 %i",size(X_v0))
end


end