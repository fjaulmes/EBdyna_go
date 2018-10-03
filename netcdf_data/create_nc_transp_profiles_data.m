% Test writing and reading NetCDF4 with chunking and compression using
% Native Matlab routines (tested in Matlab 2010b)
% Rich Signell (rsignell@usgs.gov)
load('../data_tokamak/physics_constants.mat')
load('../data_tokamak/flux_geometry.mat')
load('../data_tokamak/pressure_profile.mat')
load('../data_tokamak/q_profile.mat')

rho_tor_scale=sqrt(tor_flux_profile/max(tor_flux_profile));
psi_norm_scale=1-(psi_scale/max(psi_scale));
rho_pol_scale=sqrt(psi_norm_scale);

psi_bar_rho=interp1(rho_tor_scale,psi_norm_scale,(0:256)/256);
r_bar_rho=interp1(rho_tor_scale,radial_r_value_flux,(0:256)/256);
rho_pol_rho_tor=interp1(rho_tor_scale,rho_pol_scale,(0:256)/256);
rho_lin_scale=(0:256)/256;
psi_bar_rho_lin=interp1((0:256)/256,psi_norm_scale,rho_tor_scale);

filename='31557A01.CDF'

% 
% read NetCDF4 file
ncid=netcdf.open(filename,'nowrite');
varid=netcdf.inqVarID(ncid,'ND');
ND_profile=netcdf.getVar(ncid,varid);
varid=netcdf.inqVarID(ncid,'NE');
NE_profile=netcdf.getVar(ncid,varid);
varid=netcdf.inqVarID(ncid,'TE');
TE_profile=netcdf.getVar(ncid,varid);
varid=netcdf.inqVarID(ncid,'TI');
TI_profile=netcdf.getVar(ncid,varid);
varid=netcdf.inqVarID(ncid,'UTOTL');
UTOTL_profile=netcdf.getVar(ncid,varid);
varid=netcdf.inqVarID(ncid,'UPHI');
UPHI_profile=netcdf.getVar(ncid,varid);
varid=netcdf.inqVarID(ncid,'UBPOL');
UBPOL_profile=netcdf.getVar(ncid,varid);
varid=netcdf.inqVarID(ncid,'OMEGA');
OMEGA_profile=netcdf.getVar(ncid,varid);
varid=netcdf.inqVarID(ncid,'TIME');
TIME_profile=netcdf.getVar(ncid,varid);
varid=netcdf.inqVarID(ncid,'Q');
Q_profile=netcdf.getVar(ncid,varid);
varid=netcdf.inqVarID(ncid,'PTOWB');
PTOWB_profile=netcdf.getVar(ncid,varid);
varid=netcdf.inqVarID(ncid,'NIMP');
NIMP_profile=netcdf.getVar(ncid,varid);
% ncinfo(filename,'NIMP')

netcdf.close(ncid);

CRASH_TIME_STAMP=91
DELTA_CRASH=2

TE_profile_interp_ini=interp1((0:39)/39,TE_profile(:,CRASH_TIME_STAMP),rho_lin_scale);
TI_profile_interp_ini=interp1((0:39)/39,TI_profile(:,CRASH_TIME_STAMP),rho_lin_scale);
NE_profile_interp_ini=interp1((0:39)/39,NE_profile(:,CRASH_TIME_STAMP),rho_lin_scale)*1e6;
NI_profile_interp_ini=interp1((0:39)/39,ND_profile(:,CRASH_TIME_STAMP),rho_lin_scale)*1e6;
OMEGA_profile_interp_ini=interp1((0:39)/39,OMEGA_profile(:,CRASH_TIME_STAMP),rho_lin_scale);
Q_profile_interp_ini=interp1((0:39)/39,Q_profile(:,CRASH_TIME_STAMP),rho_lin_scale);
PTOWB_profile_interp_ini=interp1((0:39)/39,PTOWB_profile(:,CRASH_TIME_STAMP),rho_lin_scale);
NIMP_profile_interp_ini=interp1((0:39)/39,NIMP_profile(:,CRASH_TIME_STAMP),rho_lin_scale)*1e6;

TE_profile_interp_end=interp1((0:39)/39,TE_profile(:,CRASH_TIME_STAMP+DELTA_CRASH),rho_lin_scale);
TI_profile_interp_end=interp1((0:39)/39,TI_profile(:,CRASH_TIME_STAMP+DELTA_CRASH),rho_lin_scale);
NE_profile_interp_end=interp1((0:39)/39,NE_profile(:,CRASH_TIME_STAMP+DELTA_CRASH),rho_lin_scale)*1e6;
NI_profile_interp_end=interp1((0:39)/39,ND_profile(:,CRASH_TIME_STAMP+DELTA_CRASH),rho_lin_scale)*1e6;
OMEGA_profile_interp_end=interp1((0:39)/39,OMEGA_profile(:,CRASH_TIME_STAMP+DELTA_CRASH),rho_lin_scale);
Q_profile_interp_end=interp1((0:39)/39,Q_profile(:,CRASH_TIME_STAMP+DELTA_CRASH),rho_lin_scale);
PTOWB_profile_interp_end=interp1((0:39)/39,PTOWB_profile(:,CRASH_TIME_STAMP+DELTA_CRASH),rho_lin_scale);
NIMP_profile_interp_end=interp1((0:39)/39,NIMP_profile(:,CRASH_TIME_STAMP+DELTA_CRASH),rho_lin_scale)*1e6;

N_profile_interp_ini=NE_profile_interp_ini+NI_profile_interp_ini;
T_profile_interp_ini=TE_profile_interp_ini+TI_profile_interp_ini;
N_profile_interp_end=NE_profile_interp_end+NI_profile_interp_end;
T_profile_interp_end=TE_profile_interp_end+TI_profile_interp_end;

PBULK_profile_interp_ini=(NE_profile_interp_ini.*TE_profile_interp_ini)+(NI_profile_interp_ini.*TI_profile_interp_ini);
PBULK_profile_interp_end=(NE_profile_interp_end.*TE_profile_interp_end)+(NI_profile_interp_end.*TI_profile_interp_end);
PFAST_ini=PTOWB_profile_interp_ini-PBULK_profile_interp_ini*eV;
PFAST_end=PTOWB_profile_interp_end-PBULK_profile_interp_end*eV;

PTOT_profile_interp_ini=PTOWB_profile_interp_ini;


PTOT_rho_scale=interp1(rho_lin_scale,PTOWB_profile_interp_ini,rho_tor_scale);


save('../data_tokamak/q_profile.mat','-append','rho_tor_scale','rho_pol_scale','psi_bar_rho','r_bar_rho',...
     'Q_profile_interp_ini','Q_profile_interp_end','CRASH_TIME_STAMP');
save('../data_tokamak/pressure_profile.mat','-append','rho_tor_scale','rho_pol_scale','psi_bar_rho',...
    'TE_profile_interp_ini','TI_profile_interp_ini','NE_profile_interp_ini','NI_profile_interp_ini',...
    'OMEGA_profile_interp_ini','OMEGA_profile_interp_end','PTOT_profile_interp_ini','NIMP_profile_interp_ini','PFAST_ini','PFAST_end');

rho_tor_scale=rho_lin_scale;
KARELTOOLFILE='AUG31557_2p25_F2P_data.mat'
save(KARELTOOLFILE,'-append', 'rho_tor_scale')
save(KARELTOOLFILE,'-append', 'PTOT_profile_interp_ini')
save(KARELTOOLFILE,'-append', 'Q_profile_interp_ini')
rho_tor_scale=sqrt(tor_flux_profile/max(tor_flux_profile));

%%
figure(1)


subplot(2,1,1)
set(gca,'fontsize',20)
hold on
grid on

plot(TIME_profile,NE_profile(1,:)*1e6,'b','linewidth',2)
xlim([2 3.0])
ylabel('n_e')
xlabel('t (s)')

subplot(2,1,2)
set(gca,'fontsize',20)
hold on
grid on

plot(TIME_profile,TI_profile(1,:),'r','linewidth',2)
xlabel('t (s)')
ylabel('T_i')
xlim([2 3.0])
%%
figure(3)
hold on
grid on
plot(rho_lin_scale,Q_profile_interp_ini,'b','linewidth',2)
plot(rho_lin_scale,Q_profile_interp_end,'r','linewidth',2)
xlabel('\rho _{tor}')
ylabel('q')
plot(rho_tor_scale,1.005*q_initial_profile,'g--')
%%

figure(9)
set(gca,'fontsize',20)
pcore=PTOWB_profile_interp_ini(1)
% Prho_ini=interp1(rho_lin_scale,PTOWB_profile_interp_ini,rho_tor_scale);

hold on
grid on
plot(rho_lin_scale,PTOWB_profile_interp_ini,'b','linewidth',2)
plot(rho_lin_scale,PTOWB_profile_interp_end,'r','linewidth',2)

% plot(rho_lin_scale,PTOWB_profile_interp_end/pcore,'r','linewidth',2)
xlabel('\rho _{tor}')
ylabel('Ptot (Pa)')

hold on
% plot(psi_bar_rho,Prho_ini,'b','linewidth',2)
plot(rho_tor_scale,P_initial_profile,'g--','linewidth',2)
% pol_P=[ 1. -4.6 10.8 -12.2 5.2 -0.2  ];
% polytest=fliplr(pol_P);
% polytest=[-0.09 4.6 -11.8 11 -4.7 1.0];
pol_P=[ 1 -2.55 2.37  4.84 -19  22 -8.65 ]
polytest=fliplr(pol_P);
% polytest=[0 -5.0 24 -31 17 -5.4 1.0];
% polytest=[-2.55 4 0.98 -3.2 -0.2 1.0];
polyvalues=polyval(polytest,psi_norm_scale);
% plot(psi_norm_scale,polyvalues,'r--','linewidth',2)
% plot(1-(psi_scale/max(psi_scale)),P_initial_profile/P_initial_profile(1),'g--','linewidth',2)

%%

figure(2)

hold on
grid on
CORE_ROTATION=OMEGA_profile_interp_ini(1)
plot(psi_bar_rho,(OMEGA_profile_interp_ini),'b','linewidth',2)
plot(psi_bar_rho,(OMEGA_profile_interp_end),'r','linewidth',2)

% plot(rho_lin_scale,(OMEGA_profile_interp_ini),'b','linewidth',2)
% plot(rho_lin_scale,(OMEGA_profile_interp_end),'r','linewidth',2)


xlabel('\psi _{pol}')
ylabel('\Omega (rad/s)')
% 
% hold on
% plot(1-psi_scale/max(psi_scale),NE_profile_interp_ini)
% plot(1-psi_scale/max(psi_scale),NE_profile_interp_end,'r')

% plot(1-psi_scale/max(psi_scale),R0*OMEGA_profile_interp_ini*1e-7/(R0*OMEGA_profile_interp_ini(1)*1e-7),'b','linewidth',2)


%%
figure(4)
hold on
grid on
PFAST_ini=PTOWB_profile_interp_ini-PBULK_profile_interp_ini*eV;
PFAST_end=PTOWB_profile_interp_end-PBULK_profile_interp_end*eV;

PFAST_ini_norm=max(PFAST_ini);
% PFAST_ini=PFAST_ini/PFAST_ini_norm;
% PFAST_end=PFAST_end/PFAST_ini_norm;

plot(rho_lin_scale,PFAST_ini,'b','linewidth',2)
plot(rho_lin_scale,PFAST_end,'r','linewidth',2)
xlabel('\rho _{tor}')
ylabel('Ptot (Pa)')



%%

figure(8)
subplot(2,1,1)
set(gca,'fontsize',18)
hold on
grid on
plot(rho_lin_scale,Q_profile_interp_ini,'b','linewidth',3)
plot(rho_lin_scale,Q_profile_interp_end,'r','linewidth',3)
xlabel('\rho _{tor}')
ylabel('q')
plot(rho_tor_scale,1.0*q_initial_profile,'g--','linewidth',3)
xlim([0 1.0])

legend('before TRANSP','after TRANSP','EBdyna')

subplot(2,1,2)
set(gca,'fontsize',18)
hold on
grid on
plot(rho_lin_scale,PTOWB_profile_interp_ini,'b','linewidth',3)
plot(rho_lin_scale,PTOWB_profile_interp_end,'r','linewidth',3)
xlabel('\rho _{tor}')
ylabel('Ptot (Pa)')
hold on
plot(rho_tor_scale,P_initial_profile,'g--','linewidth',3)
xlim([0 1.0])
