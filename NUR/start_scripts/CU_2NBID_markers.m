
format compact
disp('------------------------------------------')
RTAN=40

RTAN_m=RTAN*0.01;

%% === Configuration ===
SHOT_NUMBER = 5400;
TIME_EQUIL = 1.15;         % Equilibrium time in seconds
sim_folder_name='CU_5400_1150_RT42_01'
DATA_INPUT_FOLDER=['./',sim_folder_name,'/input/'];
DATA_PLASMA_FOLDER=['./',sim_folder_name,'/data_plasma/']
DATA_PHYS_FOLDER=['./data_common/physics_data/'];
Zeff_value=1.2;


%% read tokamak profiles

load([DATA_PLASMA_FOLDER,'pressure_profile.mat'])
load([DATA_PLASMA_FOLDER,'XZsmall_fields_tokamak_pre_collapse.mat'])
load([DATA_PLASMA_FOLDER,'motions_map_dimensions.mat'])


%R0=0.55;%C
R0=0.894;%CU

psin_XZ_map=psi_n_map;
% psin_XZ_map=(psi_map-psi_scale(1))/range(psi_scale);
% figure;imagesc(Z_scale,R0+X_scale,psi_XZ_map)

ne_f=fit([psi_norm ],[ne_prof ]','cubicinterp');
Te_f=fit([psi_norm ],[Te_prof ]','cubicinterp');
Zeff_f=fit([psi_norm],Zeff_value*ones(size(ne_prof))','cubicinterp');

ne2d=ne_f(psin_XZ_map);
ne2d=reshape(ne2d,size(psin_XZ_map));
ne2d(isnan(ne2d))=0;
ne2d=ne2d.*(psin_XZ_map<=1);
% ne2d(psin_XZ_map>1)=0;

Te2d=Te_f(psin_XZ_map);
Te2d=reshape(Te2d,size(psin_XZ_map));
Te2d(isnan(Te2d))=0;
Te2d=Te2d.*(psin_XZ_map<=1);
% Te2d(psin_XZ_map>1)=0;

Zeff2d=Zeff_f(psin_XZ_map);
Zeff2d=reshape(Zeff2d,size(psin_XZ_map));
Zeff2d(isnan(Zeff2d))=0;
Zeff2d=Zeff2d.*(psin_XZ_map<=1);
% Zeff2d(psin_XZ_map>1)=0;

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

%beam properties and beamduct
DIST_VALVE_TANK=0.6
DIST_SOURCE_VALVE=2.47+0.47

Beam1.bmwid.r=0.106;  %source/(last grid) half width [m]
Beam1.bmwid.z=0.106;  %source/(last grid) half width [m]
Beam1.bmwid.opt=0; %option shape 0-circ, 1-rect
Beam1.focl.r=4.6; %distance from the source to optical focual point [m]
Beam1.focl.z=4.6; %distance from the source to optical focual point [m]
Beam1.div.r=11e-3; %source divergence [rad]
Beam1.div.z=11e-3; %source divergence [rad]

Beam1.tilt_omp=5.5; % Beam inclination with mid-plane in degrees

Beam1.R_t=0.42; %Tangency radius [m]
Beam1.Z_s=0.505; %source elevation over midplane [m]
Beam1.d_s=(DIST_SOURCE_VALVE+DIST_VALVE_TANK+sqrt(1.67^2-Beam1.R_t^2)+0.62)/cos(Beam1.tilt_omp/180*pi);% distance from the source to the tangency point [m] : 5.88 m
% Beam1.Z_a=sign(Beam1.Z_s)*(abs(Beam1.Z_s)-Beam1.app.l(1)*sin(4.59/180*pi)); % beam aperture elevation over midplane [m]

Beam1.app.l(1)=Beam1.Z_s/tand(Beam1.tilt_omp); %distance from the source to the apperture [m] : range 3.456 - 4.725 m
Beam1.app.r(1)=0.4;% apperture half width [m]
Beam1.app.z(1)=0.4;% apperture half width [m]
Beam1.app.opt(1)=1;% option shape 0-circ, 1-rect, 2-circ pyram, 3-pyramid
Beam1.Z_a=sign(Beam1.Z_s)*(abs(Beam1.Z_s)-Beam1.app.l(1)*sin(Beam1.tilt_omp/180*pi)); % beam aperture elevation over midplane [m] ~ 0


Beam1.Ab=2.01355; %beam species atomic number
Beam1.E=80e3; %beam Energy[eV]

Beam1.app_ang=0; % leave at 0: normal aperture
Beam1.phi=45;   % toroidal position of this NBI source: degrees clockwise wrt. to noon!

Beam1.frac_E1=0.465;
Beam1.frac_E2=0.3875;

Beam2=Beam1;
Beam2.Z_a=-Beam2.Z_a;
Beam2.Z_s=-Beam2.Z_s;


% Beam1.bmwid.r=0.1;  %source/(last grid) half width [m]
% Beam1.bmwid.z=0.1;  %source/(last grid) half width [m]
% Beam1.bmwid.opt=0; %option shape 0-circ, 1-rect
% Beam1.focl.r=3.4; %distance from the source to optical focual point [m]
% Beam1.focl.z=3.4; %distance from the source to optical focual point [m]
% Beam1.div.r=11e-3; %source divergence [rad]
% Beam1.div.z=11e-3; %source divergence [rad]
% 
% Beam1.app.l(1)=5; %distance from the source to the apperture [m]
% Beam1.app.r(1)=0.15;% apperture half width [m]
% Beam1.app.z(1)=0.15;% apperture half width [m]
% Beam1.app.opt(1)=1;% option shape 0-circ, 1-rect, 2-circ pyram, 3-pyramid
% 
% Beam1.R_t=RTAN_m; %Tangency radius [m]
% Beam1.d_s=(sqrt(1.6731^2-Beam1.R_t^2)+4.1275)/cos(4.59/180*pi);% distance from the source to the tangency point [m] : sqrt(1.21^2-0.3^2)+0.42+4.24
% Beam1.Z_s=0.425; %source elevation over midplane [m]
% Beam1.Z_a=sign(Beam1.Z_s)*(abs(Beam1.Z_s)-Beam1.app.l(1)*sin(4.59/180*pi)); % beam aperture elevation over midplane [m]
% Beam1.Ab=2.01355; %beam species atomic number
% Beam1.E=80e3; %beam Energy[eV]
% Beam1.frac_E1=0.465;
% Beam1.frac_E2=0.3875;
% Beam1.app_ang=0;
% Beam1.phi=45;   % degrees clockwise wrt. to noon!
% 
% Beam2=Beam1;
% Beam2.Z_a=-Beam2.Z_a;
% Beam2.Z_s=-Beam2.Z_s;

param.nmark=3e5; %number of markers
param.nb=2e4; %number of markers creted during one iteration (have toke it because of shinethrough)
param.nr=1000; %number of points along one ray to calculate probability function

nea=mean(ne_prof);
Plasma.ne=ne2d';%/nea*3.5e19; %density [m^-3]
Plasma.Te=Te2d'; %[eV]
Plasma.Zeff=ones(size(Te2d')) * 1.3;
Plasma.R=R_grid; %[m] 1d array
Plasma.Z=Z_grid ;%[m] 1d array
Plasma.Ap=2; %plasma speace atomic number (int)
Plasma.IMP=1; %impurity index (for now IMP==1 CARBON impurity)
Plasma.BT=1; %Magnetic field
Plasma.Rwall_min=0.599; %Central column

disp('------------------------------------------')

rng('shuffle')
[R_i0,Z_i0,phi_i0,vr0,vz0,vphi0,E_i0,Shth0,X0,Y0,Z0,sh0,rcoord0]=NUR(Beam1,Plasma,param);

rcoord0.R=sqrt(rcoord0.X.^2+rcoord0.Y.^2);
rcoord0.phi=0.5*(sign(rcoord0.Y)+1).*asin(rcoord0.X./rcoord0.R)-0.5*(sign(rcoord0.Y)-1).*(acos(rcoord0.X./rcoord0.R)+pi/2);
rcoord0.phi=rcoord0.phi+Beam1.phi/180*pi;
rcoord0.phi=mod(rcoord0.phi,2*pi);
rcoord0.X=rcoord0.R.*sin(rcoord0.phi);
rcoord0.Y=rcoord0.R.*cos(rcoord0.phi);

phi_i0=phi_i0+Beam1.phi/180*pi;


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
%         rcoord.X=[rcoord0.X];
%         rcoord.Y=[rcoord0.Y];
%         rcoord.Z=[rcoord0.Z];
rng('shuffle')

[R_i0,Z_i0,phi_i0,vr0,vz0,vphi0,E_i0,Shth0,X0,Y0,Z0,sh0,rcoord0]=NUR(Beam2,Plasma,param);

rcoord0.R=sqrt(rcoord0.X.^2+rcoord0.Y.^2);
rcoord0.phi=0.5*(sign(rcoord0.Y)+1).*asin(rcoord0.X./rcoord0.R)-0.5*(sign(rcoord0.Y)-1).*(acos(rcoord0.X./rcoord0.R)+pi/2);
rcoord0.phi=rcoord0.phi+Beam2.phi/180*pi;
rcoord0.phi=mod(rcoord0.phi,2*pi);
rcoord0.X=rcoord0.R.*sin(rcoord0.phi);
rcoord0.Y=rcoord0.R.*cos(rcoord0.phi);

phi_i0=phi_i0+Beam2.phi/180*pi;



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

%

%

%%
FILENAME=[DATA_INPUT_FOLDER,'NUR_out_CU_2NBI_phi45_Rt' num2str(RTAN) '.mat']

NUR_out=save_NUR_output(FILENAME,Beam1,Plasma,param,R_i,Z_i,phi_i,vr,vz,vphi,Shth,rcoord)

%plot ions vectors
nv=10000; %how many to plot
%pickup rundom numbers

% figure;quiver(R_i(rann),Z_i(rann),vr(rann),vz(rann))
% figure;quiver(R_i(rann),phi_i(rann),vr(rann),vphi(rann))
% figure;quiver(R_i(rann).*sin(phi_i(rann)),R_i(rann).*cos(phi_i(rann)),vphi(rann).*cos(phi_i(rann))+vr(rann).*sin(phi_i(rann)),-vphi(rann).*sin(phi_i(rann))+vr(rann).*cos(phi_i(rann)))
%%
figure
rann = randi([1 numel(Z_i)],1,nv)

% quiver3(R_i(rann).*sin(phi_i(rann)),R_i(rann).*cos(phi_i(rann)),Z_i(rann)...
%     ,vphi(rann).*cos(phi_i(rann))+vr(rann).*sin(phi_i(rann)),-vphi(rann).*sin(phi_i(rann))+vr(rann).*cos(phi_i(rann)),vz(rann),4)
load(FILENAME)
quiver3(NUR_out.R_i(rann).*sin(NUR_out.phi_i(rann)),NUR_out.R_i(rann).*cos(NUR_out.phi_i(rann)),NUR_out.Z_i(rann)...
    ,NUR_out.vphi(rann).*cos(NUR_out.phi_i(rann))+NUR_out.vR(rann).*sin(NUR_out.phi_i(rann)),-NUR_out.vphi(rann).*sin(NUR_out.phi_i(rann))+NUR_out.vR(rann).*cos(NUR_out.phi_i(rann)),NUR_out.vZ(rann),4)

% phipos=phi_i-pi/2;
% quiver(R_i(rann).*sin(phi_i(rann)),R_i(rann).*cos(phi_i(rann)),vr(rann).*cos(phipos(rann))-vphi(rann).*sin(phipos(rann)),-vphi(rann).*cos(phipos(rann))-vr(rann).*sin(phipos(rann)),16)

% figure;quiver3(R_i(rann).*sin(phi_i(rann)),R_i(rann).*cos(phi_i(rann)),Z_i(rann)...
%     ,vphi(rann).*cos(phi_i(rann))+vr(rann).*sin(phi_i(rann)),-vphi(rann).*sin(phi_i(rann))+vr(rann).*cos(phi_i(rann)),vz(rann))
axis equal

read_NUR_markers