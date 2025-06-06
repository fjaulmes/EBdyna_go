%% read tokamak profiles
%load('/compass/Shared/Common/Theory-Modelling/NBI_FIESTA/EBdyna/data_tokamak/pressure_profile.mat')
%load('/compass/Shared/Common/Theory-Modelling/NBI_FIESTA/EBdyna/data_tokamak/flux_geometry.mat')
%load('/compass/Shared/Common/Theory-Modelling/NBI_FIESTA/EBdyna/data_tokamak/tokamak_map_dimensions.mat')

load('./data_tokamak/pressure_profile.mat')
load('./data_tokamak/flux_geometry.mat')
load('./data_tokamak/tokamak_map_dimensions.mat')

% R0=0.55;%C
% R0=0.894;%CU

psi_XZ_map=(psi_XZ_map-psi_scale(1))/range(psi_scale);
% figure;imagesc(Z_scale,R0+X_scale,psi_XZ_map)

nef=fit([psi_norm ]',[ne_prof ]','cubicinterp');
Tef=fit([psi_norm ]',[Te_prof ]','cubicinterp');
Zefff=fit([psi_norm]',[linspace(0,2,numel(ne_prof))]','cubicinterp');

ne2d=nef(psi_XZ_map);
ne2d=reshape(ne2d,size(psi_XZ_map));
ne2d(isnan(ne2d))=0;
ne2d(psi_XZ_map>1)=0;

Te2d=Tef(psi_XZ_map);
Te2d=reshape(Te2d,size(psi_XZ_map));
Te2d(isnan(Te2d))=0;
Te2d(psi_XZ_map>1)=0;

Zeff2d=Zefff(psi_XZ_map);
Zeff2d=reshape(Zeff2d,size(psi_XZ_map));
Zeff2d(isnan(Zeff2d))=0;
Zeff2d(psi_XZ_map>1)=0;

% make structures

% Beam.bmwidra=0.1;  %source/(last grid) half width [m]
% Beam.foclra=4; %distance from the source to optical focual point [m]
% Beam.divra=2*0.00823688265; %source divergence [rad]
% Beam.app(1)=1; %distance from the source to the apperture [m]
% Beam.app(2)=1;% apperture half width [m]
% Beam.R_t=0.65; %Tangency radius [m]
% Beam.d_s=4.773;% distance from the source to the tangency point [m]
% Beam.Z_s=0.41; %source elevation over midplane [m]
% Beam.Z_a=0.322844; %aperture elevation over midplane [m]
% Beam.Ab=2.01355; %beam speace atomic number
% %Beam.Ab=2;
% Beam.E=80e3;%beam Energy[eV]
% Beam.frac_E1=0.58;
% Beam.frac_E2=0.26;
% 
% 
% param.nmark=1e5; %number of markers
% param.nb=10000; %number of markers creted during one itaration (have toke it because of shinethrough)
% param.nr=1000; %number of points along one ray to calculate probability function


% So particle fractions injected are:
% 0.4650   (80 keV)
% 0.3875   (40 keV)
% 0.3487   (26.7 keV)





Beam.R_t=0.55; %Tangency radius [m]
Beam.d_s=5.74;% distance from the source to the tangency point [m] : sqrt(1.21^2-0.55^2)+0.42+4.24
Beam.Z_s=0.0; %source elevation over midplane [m]
Beam.Z_a=0.0; % beam aperture elevation over midplane [m]
Beam.Ab=2.01355; %beam species atomic number
Beam.E=75e3; %beam Energy[eV]
Beam.frac_E1=0.465;
Beam.frac_E2=0.3875;



Beam.bmwid.r=0.1;  %source/(last grid) half width [m]
Beam.bmwid.z=0.1;  %source/(last grid) half width [m]
Beam.bmwid.opt=0; %option shape 0-circ, 1-rect
Beam.focl.r=3.4; %distance from the source to optical focal point [m]
Beam.focl.z=3.4; %distance from the source to optical focal point [m]
Beam.div.r=8e-3; %source divergence [rad]
Beam.div.z=8e-3; %source divergence [rad]

Beam.app = struct();
Beam.app.l(1)=4.25; %distance from the source to the apperture [m]
Beam.app.r(1)=0.1;% apperture half width [m]
Beam.app.z(1)=0.1;% apperture half width [m]
Beam.app.opt(1)=1;% option shape 0-circ, 1-rect, 2-circ pyram, 3-pyramid

Beam.app_ang=0;
Beam.phi=0;

param.nmark=1e5; %number of markers
param.nb=8e3; %number of markers creted during one itaration (have toke it because of shinethrough)
param.nr=1500; %number of points along one ray to calculate probability function

nea=mean(ne_prof);


nea=mean(ne_prof);
Plasma.ne=ne2d';%/nea*3.5e19; %density [m^-3]
Plasma.Te=Te2d'; %[eV]
Plasma.Zeff=ones(size(Te2d')) * 1.5;
Plasma.R=X_scale+R0; %[m] 1d array
Plasma.Z=Z_scale ;%[m] 1d array
Plasma.Ap=2; %plasma speace atomic number (int)
Plasma.IMP=1; %impurity index (for now IMP==1 CARBON impurity)
Plasma.BT=1; %Magnetic field

rng('shuffle')

[R_i,Z_i,phi_i,vR,vZ,vphi,E_i,Shth,X,Y,Z,sh,rcoord]=NUR(Beam,Plasma,param);
% 
% 

%%
mu =    1.66053906660 * 1e-27;
mD =    2.01410177811 * mu;
eV =    1.602176634 * 1e-19;


v1=sqrt(2*Beam.E*(eV/mD));
v2=sqrt(2*Beam.E/2*(eV/mD));
v3=sqrt(2*Beam.E/3*(eV/mD));

vnorm=sqrt(vR.^2+vZ.^2+vphi.^2);
Enorm=0.5*(mD/eV)*vnorm.^2;

Emain=find(abs(Beam.E-Enorm)<1e3);
E2=find(abs(Beam.E/2-Enorm)<1e3);
E3=find(abs(Beam.E/3-Enorm)<1e3);
vR(Emain)=vR(Emain).*v1./vnorm(Emain);
vR(E2)=vR(E2).*v2./vnorm(E2);
vR(E3)=vR(E3).*v3./vnorm(E3);
vZ(Emain)=vZ(Emain).*v1./vnorm(Emain);
vZ(E2)=vZ(E2).*v2./vnorm(E2);
vZ(E3)=vZ(E3).*v3./vnorm(E3);
vphi(Emain)=vphi(Emain).*v1./vnorm(Emain);
vphi(E2)=vphi(E2).*v2./vnorm(E2);
vphi(E3)=vphi(E3).*v3./vnorm(E3);

NUR_out=struct();
NUR_out.R_i=R_i;
NUR_out.Z_i=Z_i;
NUR_out.phi_i=phi_i;
NUR_out.vR=vR;
NUR_out.vZ=vZ;
NUR_out.vphi=vphi;
NUR_out.Shth=Shth;

INDEXES_NUR=randperm(length(NUR_out.R_i));

NUR_out.R_i=NUR_out.R_i(INDEXES_NUR);
NUR_out.Z_i=NUR_out.Z_i(INDEXES_NUR);
NUR_out.phi_i=NUR_out.phi_i(INDEXES_NUR);
NUR_out.vR=NUR_out.vR(INDEXES_NUR);
NUR_out.vZ=NUR_out.vZ(INDEXES_NUR);
NUR_out.vphi=NUR_out.vphi(INDEXES_NUR);


rcoord.X=rcoord.X(:,INDEXES_NUR);
rcoord.Y=rcoord.Y(:,INDEXES_NUR);
rcoord.Z=rcoord.Z(:,INDEXES_NUR);

save('NUR_output_COMPASS_Rtan55_100k.mat','NUR_out','Plasma','Beam','param','rcoord')


%%


FILENAME='NUR_out_COMPASS.mat'

NUR_out=save_NUR_output(FILENAME,Beam,Plasma,param,R_i,Z_i,phi_i,vR,vZ,vphi,Shth,rcoord)





%plot ions vectors
nv=10000; %how many to plot
%pickup rundom numbers
rann = randi([1 numel(Z_i)],1,nv)

% figure;quiver(R_i(rann),Z_i(rann),vr(rann),vz(rann))
% figure;quiver(R_i(rann),phi_i(rann),vr(rann),vphi(rann))
% figure;quiver(R_i(rann).*sin(phi_i(rann)),R_i(rann).*cos(phi_i(rann)),vphi(rann).*cos(phi_i(rann))+vr(rann).*sin(phi_i(rann)),-vphi(rann).*sin(phi_i(rann))+vr(rann).*cos(phi_i(rann)))

figure;quiver3(R_i(rann).*sin(phi_i(rann)),R_i(rann).*cos(phi_i(rann)),Z_i(rann)...
    ,vphi(rann).*cos(phi_i(rann))+vR(rann).*sin(phi_i(rann)),-vphi(rann).*sin(phi_i(rann))+vR(rann).*cos(phi_i(rann)),vZ(rann))
axis equal
