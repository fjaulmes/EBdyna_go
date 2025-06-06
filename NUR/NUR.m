function [R_io,Z_io,phi_io,vro,vzo,vphio,E_io,Shth_i,Xwo,Ywo,Zwo,sho,rcoord]=NUR(Beam,Plasma,param)
%% input discription


% Beam.bmwidra source/(last grid) half width [m]
% Beam.foclra distance from the source to optical focual point [m]
% Beam.divra source divergence [rad]
% Beam.div.r = 0.00823688265;                        % ion source horizontal divergence (radians)
% Beam.div.z = 0.00823688265;                        % ion source horizontal divergence (radians)
% Beam.bmwid.r = 0.1;                      % ion source half-width (m)
% Beam.bmwid.z = 0.3;                      % ion source half-width (m)
% Beam.bmwid.opt = 1;                      % ion source half-width (m)
% Beam.focl.r = 4;                     % ion source horizontal focal length (m)
% Beam.focl.z = 4;                     % ion source horizontal focal length (m)
% Beam.app.l(1)=3;                 % app aperture to source distance
% Beam.app.l(2)=3.2;                 % app aperture to source distance
% Beam.app.r(1)=0.07;                  % app(2) apperture radius
% Beam.app.r(2)=0.05;                  % app(2) apperture radius
% Beam.app.z(1)=100;                  % app(2) apperture radius
% Beam.app.z(2)=0.1;                  % app(2) apperture radius
% Beam.app.opt(1)=0;                   % option: 0 -circ, 1 - rect
% Beam.app.opt(2)=1;                   % option: 0 -circ, 1 - rect
% 
% Beam.R_t Tangency radius [m]
% Beam.app_ang last apperture center toroidal angle coordinate.
% Beam.d_s distance from the source to the tangency point [m]
% Beam.Z_s source elevation over midplane [m]
% Beam.Z_a aperture elevation over midplane [m]
% Beam.Ab beam speace atomic number
% Beam.E beam Energy[eV]
% Beam.frac_E1 full energy current fraction
% Beam.frac_E2 half energy current fraction


% param.nmark number of markers
% param.nb number of markers creted during one itaration (have toke it because of shinethrough)
% param.nr number of points along one ray to calculate probability function

% Plasma.ne density [m^-3]
% Plasma.Te [eV]
% Plasma.Zeff 
% Plasma.R [m] 1d array
% Plasma.Z [m] 1d array
% Plasma.Ap plasma speace atomic number
% Plasma.IMP impurity index (for now IMP==1 CARBON impurity)
% Plasma.BT Magnetic field [T]

%% The code
Z_i=[];R_i=[];phi_i=[];
Z_io=[];R_io=[];phi_io=[];

if(Plasma.BT<4.5) 
    BT=1;
else
    BT=5;
end


Shth=0;
i=0;
rcoord=struct();
% rcoord.X=[];
% rcoord.Y=[];
% rcoord.Z=[];


for j=1:3
X0=[];
Y0=[];
Z0=[];
    if (j==1)
        E=Beam.E;
        nmark=round(param.nmark*Beam.frac_E1);
    elseif (j==2)
        E=Beam.E/2;
        nmark=round(param.nmark*Beam.frac_E2);
    else
        E=Beam.E/3;
        nmark=param.nmark-round(param.nmark*Beam.frac_E1)-round(param.nmark*Beam.frac_E2);
    end

sAold=0;
sA=0;
while (sA<nmark)

    i=i+1;
fprintf('Energy number %i  \n\n',j);
fprintf('iteration %i  \n\n',i);

    time1=cputime;
    
tic
% [R,Z,phi,dl,X,Y]=bcoord(Beam.bmwidra,Beam.foclra,Beam.divra,param.nb,Beam.app,Beam.R_t,Beam.d_s,Beam.Z_s,Beam.Z_a,min(Plasma.R),max(Plasma.R),param.nr);
[R,Z,phi,dl,X,Y]=bcoord(Beam.bmwid,Beam.focl,Beam.div,param.nb,Beam.app,Beam.R_t,Beam.d_s,Beam.Z_s,Beam.Z_a,Plasma.Rwall_min,max(Plasma.R),param.nr);


fprintf("bcoord() ");
toc
fprintf("\n");
tic
[ne_b,Te_b,Zeff_b]=bline_profiles(R,Z,Plasma.ne,Plasma.Te/1e3,Plasma.Zeff,Plasma.R,Plasma.Z);
fprintf("bline_profiles() ");
toc
fprintf("\n");

tic
[P, sigma, dP,I]=deposition_pbty(dl,E/1e3,ne_b,Te_b,Zeff_b,Beam.Ab,Plasma.Ap,Plasma.IMP,BT);
fprintf("deposition_pbty() ");
toc
fprintf("\n");


Shth=(i-1)/i*Shth+sum(1-P(end,:))/size(P,2)/i;



fprintf('Shine through of current energy component %7.5f percent ',Shth*100);
fprintf("\n\n");
tic
[X_ii,Y_ii,Z_ii,vx,vy,vz_i,REMOVED_LINES]=random_position_gen(P,X,Y,Z,E/Beam.Ab);



fprintf("random_position_gen() ");
toc
fprintf("\n");

X_reduced=X(:,~REMOVED_LINES);
Y_reduced=Y(:,~REMOVED_LINES);
Z_reduced=Z(:,~REMOVED_LINES);
X_reduced=[X_reduced(1,:) ; X_reduced(end,:)];
Y_reduced=[Y_reduced(1,:) ; Y_reduced(end,:)];
Z_reduced=[Z_reduced(1,:) ; Z_reduced(end,:)];


R_ii=hypot(X_ii,Y_ii);
% phi_ii=atan(X_ii./Y_ii)
phi_ii=0.5*(sign(Y_ii)+1).*asin(X_ii./R_ii)-0.5*(sign(Y_ii)-1).*(-pi/2+asin(Y_ii./R_ii));
% phi_ii=0.5*(sign(Y_ii)+1).*asin(X_ii./R_ii)-0.5*(sign(Y_ii)-1).*(acos(X_ii./R_ii)+pi/2);

ksi=asin(vy./hypot(vx,vy));
% ksi=atan(vy./vx);

vr_i=hypot(vx,vy).*sin(phi_ii+ksi);
vphi_i=hypot(vx,vy).*cos(phi_ii+ksi);

if (i==1)
    R_i=R_ii;
    phi_i=phi_ii;
    Z_i=Z_ii;

    vr=vr_i;
    vphi=vphi_i;

    vz=vz_i;

    sh=1-P(end,:);
    Xw=X(end,:);
    Yw=Y(end,:);
    Zw=Z(end,:);
    X0=X;
    Y0=Y;
    Z0=Z;
    
    X_ray=X_reduced;
    Y_ray=Y_reduced;
    Z_ray=Z_reduced;
else
    R_i=[R_i R_ii];
    phi_i=[phi_i phi_ii];
    Z_i=[Z_i Z_ii];

    vr=[vr vr_i];
    vphi=[vphi vphi_i];
    vz=[vz vz_i];

    sh=[sh nmark/param.nmark*E/Beam.E*(1-P(end,:))];
    Xw=[Xw X(end,:)];
    Yw=[Yw Y(end,:)];
    Zw=[Zw Z(end,:)];

    X0=[X0 X(:,:)];
    Y0=[Y0 Y(:,:)];
    Z0=[Z0 Z(:,:)];

    X_ray=[X_ray X_reduced];
    Y_ray=[Y_ray Y_reduced];
    Z_ray=[Z_ray Z_reduced];

end

sA=size(Z_i,2);

time2=cputime;
trem=ceil((param.nmark-sA-numel(R_io))/(sA-sAold))*(time2-time1)/60;

fprintf('time remaining %8.3f min\n ',trem);
fprintf("\n\n");
sAold=sA;

%%
% sanity check for geometry
% figure;
% hold on
% rann=(1:100:length(phi_ii));
% plot(X_reduced(:,rann),Y_reduced(:,rann))
% plot(R_ii(rann).*sin(phi_ii(rann)),R_ii(rann).*cos(phi_ii(rann)),'r.','Markersize',15)
% plot(X_ii(rann),Y_ii(rann),'b.')
% 
% quiver(R_ii(rann).*sin(phi_ii(rann)),R_ii(rann).*cos(phi_ii(rann))...
%    ,vphi_i(rann).*cos(phi_ii(rann))+vr_i(rann).*sin(phi_ii(rann)),-vphi_i(rann).*sin(phi_ii(rann))+vr_i(rann).*cos(phi_ii(rann)),'k')
% 
% quiver(R_ii(rann).*sin(phi_ii(rann)),R_ii(rann).*cos(phi_ii(rann))...
%    ,vx(rann),vy(rann),'m')
% 
% 
% drawnow

%%
end

if ~isfield(rcoord,'X')
    rcoord.X=X_ray(:,1:nmark);
    rcoord.Y=Y_ray(:,1:nmark);
    rcoord.Z=Z_ray(:,1:nmark);
else
    rcoord.X=[rcoord.X X_ray(:,1:nmark) ];
    rcoord.Y=[rcoord.Y Y_ray(:,1:nmark) ];
    rcoord.Z=[rcoord.Z Z_ray(:,1:nmark) ];
end

if (j==1)
R_io=R_i(1:nmark);
Z_io=Z_i(1:nmark);
phi_io=phi_i(1:nmark);

vro=vr(1:nmark);
vphio=vphi(1:nmark);
vzo=vz(1:nmark);

E_io=E*ones(1,nmark);

Shth_i=nmark/param.nmark*E/Beam.E*Shth;

sho=sh(1:nmark);
Xwo=Xw(1:nmark);
Ywo=Yw(1:nmark);
Zwo=Zw(1:nmark);



else
R_io=[R_io R_i(1:nmark)];
Z_io=[Z_io Z_i(1:nmark)];
phi_io=[phi_io phi_i(1:nmark)];

vro=[vro vr(1:nmark)];
vphio=[vphio vphi(1:nmark)];
vzo=[vzo vz(1:nmark)];

E_io=[E_io E*ones(1,nmark)];  

Shth_i=Shth_i+nmark/param.nmark*E/Beam.E*Shth;

sho=[sho sh(1:nmark)];
Xwo=[Xwo Xw(1:nmark)];
Ywo=[Ywo Yw(1:nmark)];
Zwo=[Zwo Zw(1:nmark)];

% figure;
% hold on
% rann=(1:100:length(phi_ii));
% plot(rcoord.X(:,rann),rcoord.Y(:,rann))
% plot(R_io(rann).*sin(phi_io(rann)),R_io(rann).*cos(phi_io(rann)),'r.','Markersize',15)
% 
% quiver(R_io(rann).*sin(phi_io(rann)),R_io(rann).*cos(phi_io(rann))...
%    ,vphio(rann).*cos(phi_io(rann))+vro(rann).*sin(phi_io(rann)),-vphio(rann).*sin(phi_io(rann))+vro(rann).*cos(phi_io(rann)),'k')
% 
% 

end



% rcoord.X=[X0(:,1:nmark) rcoord.X  ];
% rcoord.Y=[Y0(:,1:nmark) rcoord.Y  ];
% rcoord.Z=[Z0(:,1:nmark) rcoord.Z  ];

Z_i=[];R_i=[];phi_i=[];
Shth=0;
i=0;
end

alpha=asin((Beam.Z_s-Beam.Z_a)/Beam.app.l(1));
phi_io=phi_io+Beam.app_ang+atan2(cos(alpha)*( Beam.d_s-Beam.app.l(1)),Beam.R_t);

fprintf('Shine through TOTAL %7.5f percent ',Shth_i*100);
fprintf("\n\n");


%%
% 
% figure;
% hold on
% rann=(1:1000:length(phi_io));
% plot(rcoord.X(:,rann),rcoord.Y(:,rann))
% plot(R_io(rann).*sin(phi_io(rann)),R_io(rann).*cos(phi_io(rann)),'r.')
% plot(R_io(rann).*sin(phi_io(rann)),R_io(rann).*cos(phi_io(rann)),'r.')

% 
% figure;
% hold on
% rann=(1:100:length(phi_i));
% %  rann=(100000:100012);
% plot(X0(:,rann),Y0(:,rann))
% 
% plot(R_i(rann).*sin(phi_i(rann)),R_i(rann).*cos(phi_i(rann)),'r.')
% 
% % quiver(R_i(rann).*sin(phi_i(rann)),R_i(rann).*cos(phi_i(rann))...
% %    ,vphi(rann).*cos(phi_i(rann))+vR(rann).*sin(phi_i(rann)),-vphi(rann).*sin(phi_i(rann))+vR(rann).*cos(phi_i(rann)),'k')
% % %%
% figure;
% hold on
% 
% rann=(1:2000:length(phi_io));
% % rann=(10000:10100);
% plot(rcoord.X(:,rann),rcoord.Y(:,rann))
% plot(R_io(rann).*sin(phi_io(rann)),R_io(rann).*cos(phi_io(rann)),'r.')


end






