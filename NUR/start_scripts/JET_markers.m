%% read_ tokamak _profiles

% load('/compass/Shared/Common/COMPASS-UPGRADE/RP1_Design/Scenarios/4.4/EBdyna_minorradius_28cm/data_tokamak/pressure_profile.mat')
% load('/compass/Shared/Common/COMPASS-UPGRADE/RP1_Design/Scenarios/4.4/EBdyna_minorradius_28cm/data_tokamak/flux_geometry.mat')
% load('/compass/Shared/Common/COMPASS-UPGRADE/RP1_Design/Scenarios/4.4/EBdyna_minorradius_28cm/data_tokamak/tokamak_map_dimensions.mat')

addpath '../'
load('/compass/Shared/Common/Theory-Modelling/JET/pulse_92436/EBdyna/FIESTA_equil.mat');
load('/compass/Shared/Common/Theory-Modelling/JET/pulse_92436/EBdyna/METIS_data_JET.mat')

%psi_XZ_map=(psi_XZ_map-psi_scale(1))/range(psi_scale);
psi_XZ_map=FIESTA.Psi_n;
psi_XZ_map(psi_XZ_map>1)=-1;
Z_scale1=FIESTA.scale_Z;
X_scale1=FIESTA.scale_R;
% figure;imagesc(Z_scale1,X_scale1,psi_XZ_map)

psi_norm=METIS.prof.psi_norm;
ne_prof=METIS.prof.ne;
Te_prof=METIS.prof.Te;

nef=fit([psi_norm ]',[ne_prof ]','cubicinterp');
Tef=fit([psi_norm ]',[Te_prof ]','cubicinterp');
Zefff=fit([psi_norm]',[linspace(0,2,numel(ne_prof))]','cubicinterp'); %!!!!!!!!!!!CHECK IT!!!!!!!!!!

ne2d1=nef(psi_XZ_map);
ne2d1=reshape(ne2d1,size(psi_XZ_map));
ne2d1(isnan(ne2d1))=0;
ne2d1(psi_XZ_map>1)=0;
ne2d1(psi_XZ_map<0)=0;

Te2d1=Tef(psi_XZ_map);
Te2d1=reshape(Te2d1,size(psi_XZ_map));
Te2d1(isnan(Te2d1))=0;
Te2d1(psi_XZ_map>1)=0;
Te2d1(psi_XZ_map<0)=0;

Zeff2d1=Zefff(psi_XZ_map);
Zeff2d1=reshape(Zeff2d1,size(psi_XZ_map));
Zeff2d1(isnan(Zeff2d1))=0;
Zeff2d1(psi_XZ_map>1)=0;
Zeff2d1(psi_XZ_map<0)=0;

%figure;imagesc(Z_scale1,X_scale1,Te2d)
%figure;imagesc(Z_scale1,X_scale1,ne2d)

%interpolate on another grid binded by plasma facing components
X_scale=linspace(1.806,3.9,1000);
Z_scale=linspace(-1.45,1.78,1000);

[Z_s1,X_s1]=meshgrid(Z_scale1,X_scale1);
[Z_s,X_s]=meshgrid(Z_scale,X_scale);
ne2d = interp2(Z_s1,X_s1,ne2d1,Z_s,X_s,'cubic');
Te2d = interp2(Z_s1,X_s1,Te2d1,Z_s,X_s,'cubic');
Zeff2d = interp2(Z_s1,X_s1,Zeff2d1,Z_s,X_s,'cubic');
%igure;imagesc(Z_scale,X_scale,ne2d)

%addpath('/compass/home/zadvitskiy/matlab/rw_namelist/')
indata=read_TRANSP_namelist('/compass/home/zadvitskiy/Downloads/92436H28TR.DAT');
% make structures
param.nmark=1.5e5/15; %number of markers
param.nb=1000; %number of markers creted during one itaration (have taken it because of shinethrough)
param.nr=1000; %number of points along one ray to calculate probability function

nea=mean(ne_prof);
Plasma.ne=ne2d';%/nea*3.5e19; %density [m^-3]
Plasma.Te=Te2d'; %[eV]
Plasma.Zeff=ones(size(Te2d')) ;
Plasma.R=X_scale; %[m] 1d array
Plasma.Z=Z_scale ;%[m] 1d array
Plasma.Ap=2; %plasma speace atomic number (int)
Plasma.IMP=1; %impurity index (for now IMP==1 CARBON impurity)
Plasma.BT=1; %Magnetic field



for i=1:15

Beam.bmwid.r=indata.bmwid.r(i);  %source/(last grid) half width [m]
Beam.bmwid.z=indata.bmwid.z(i);  %source/(last grid) half width [m]
Beam.bmwid.opt=indata.bmwid.opt(i); %option shape 0-circ, 1-rect
Beam.focl.r=indata.focl.r(i); %distance from the source to optical focual point [m]
Beam.focl.z=indata.focl.z(i); %distance from the source to optical focual point [m]
Beam.div.r=indata.div.r(i); %source divergence [rad]
Beam.div.z=indata.div.z(i); %source divergence [rad]
Beam.app.l=indata.app.l(i); %distance from the source to the apperture [m]
Beam.app.r=indata.app.r(i);% apperture half width [m]
Beam.app.z=indata.app.z(i);% apperture half width [m]
Beam.app.opt=1;% option shape 0-circ, 1-rect
Beam.R_t=indata.R_t(i); %Tangency radius [m]
Beam.d_s=indata.d_s(i);% distance from the source to the tangency point [m]
Beam.Z_s=indata.Z_s(i); %source elevation over midplane [m]
Beam.Z_a=indata.Z_a(i); %1st aperture elevation over midplane [m]
Beam.Ab=indata.Ab(i); %beam speace atomic number
%Beam.Ab=2;
Beam.E=100e3;
%Beam.E=Beam.Emax; %beam Energy[eV]
Beam.frac_E1=0.58;
Beam.frac_E2=0.26;
Beam.app_ang=0;
Beam.phi=indata.phi(i);

[R_i0,Z_i0,phi_i0,vr0,vz0,vphi0,E_i0,Shth0,X0,Y0,Z0,sh0,rcoord0]=NUR(Beam,Plasma,param);

phi_i0=phi_i0+indata.phi(i)/180*pi;
rcoord0.R=sqrt(rcoord0.X.^2+rcoord0.Y.^2);
rcoord0.phi=0.5*(sign(rcoord0.Y)+1).*asin(rcoord0.X./rcoord0.R)-0.5*(sign(rcoord0.Y)-1).*(acos(rcoord0.X./rcoord0.R)+pi/2);
rcoord0.phi=rcoord0.phi+Beam.phi/180*pi;
rcoord0.phi=mod(rcoord0.phi,2*pi);
rcoord0.X=rcoord0.R.*sin(rcoord0.phi);
rcoord0.Y=rcoord0.R.*cos(rcoord0.phi);

if (i==1)
R_i=[R_i0];
Z_i=[Z_i0];
phi_i=[phi_i0];
vr=[vr0];
vz=[vz0];
vphi=[vphi0];
E_i=[E_i0];
Shth=[Shth0];
X=[X0];
Y=[Y0];
Z=[Z0];
sh=[sh0];
rcoord=[rcoord0];
else
R_i=[R_i R_i0];
Z_i=[Z_i Z_i0];
phi_i=[phi_i phi_i0];
vr=[vr vr0];
vz=[vz vz0];
vphi=[vphi vphi0];
E_i=[E_i E_i0];
Shth=[Shth Shth0];
X=[X X0];
Y=[Y Y0];
Z=[Z Z0];
sh=[sh sh0];
rcoord.X=[rcoord.X rcoord0.X];
rcoord.Y=[rcoord.Y rcoord0.Y];
rcoord.Z=[rcoord.Z rcoord0.Z];

end
end



FILENAME='NUR_out_JET.mat'
NUR_out=save_NUR_output(FILENAME,Beam,Plasma,param,R_i,Z_i,phi_i,vr,vz,vphi,Shth,rcoord)


%plot ions vectors
nv=10000; %how many to plot
%pickup rundom numbers
rann = randi([1 numel(Z_i)],1,nv)

% figure;quiver(R_i(rann),Z_i(rann),vr(rann),vz(rann))
% figure;quiver(R_i(rann),phi_i(rann),vr(rann),vphi(rann))
% figure;quiver(R_i(rann).*sin(phi_i(rann)),R_i(rann).*cos(phi_i(rann)),vphi(rann).*cos(phi_i(rann))+vr(rann).*sin(phi_i(rann)),-vphi(rann).*sin(phi_i(rann))+vr(rann).*cos(phi_i(rann)))

figure;quiver3(R_i(rann).*sin(phi_i(rann)),R_i(rann).*cos(phi_i(rann)),Z_i(rann)...
    ,vphi(rann).*cos(phi_i(rann))+vr(rann).*sin(phi_i(rann)),-vphi(rann).*sin(phi_i(rann))+vr(rann).*cos(phi_i(rann)),vz(rann))
axis equal