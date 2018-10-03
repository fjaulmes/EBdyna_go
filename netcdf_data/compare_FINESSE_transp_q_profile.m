% Test writing and reading NetCDF4 with chunking and compression using
% Native Matlab routines (tested in Matlab 2010b)
% Rich Signell (rsignell@usgs.gov)
load('../data_tokamak/physics_constants.mat')
load('../data_tokamak/flux_geometry.mat')
%load('../data_tokamak/pressure_profile.mat')
rho_tor_scale=sqrt(tor_flux_profile/max(tor_flux_profile));
psi_norm_scale=1-psi_scale/max(psi_scale);
rho_pol_scale=sqrt(1-psi_scale/max(psi_scale));

% psi_bar_rho=interp1(rho_tor_scale,psi_norm_scale,(0:256)/256);
% rho_pol_rho_tor=interp1(rho_tor_scale,rho_pol_scale,(0:256)/256);
rho_lin_scale=(0:256)/256;

filename='30809C02.CDF'

% read NetCDF4 file
ncid=netcdf.open(filename,'nowrite');

varid=netcdf.inqVarID(ncid,'Q');
Q_profile=netcdf.getVar(ncid,varid);
varid=netcdf.inqVarID(ncid,'PTOWB');
PTOWB_profile=netcdf.getVar(ncid,varid);
varid=netcdf.inqVarID(ncid,'TIME');
TIME_profile=netcdf.getVar(ncid,varid);


netcdf.close(ncid);

CRASH_TIME_STAMP=60

Q_profile_interp_ini=interp1((0:39)/39,Q_profile(:,CRASH_TIME_STAMP),rho_lin_scale);
PTOWB_profile_interp_ini=interp1((0:39)/39,PTOWB_profile(:,CRASH_TIME_STAMP),rho_lin_scale);

Q_profile_interp_end=interp1((0:39)/39,Q_profile(:,CRASH_TIME_STAMP+1),rho_lin_scale);
PTOWB_profile_interp_end=interp1((0:39)/39,PTOWB_profile(:,CRASH_TIME_STAMP+1),rho_lin_scale);


%%
figure(1)
hold on
grid on
plot(rho_lin_scale,Q_profile_interp_ini,'b','linewidth',2)
% plot(rho_lin_scale,Q_profile_interp_end,'r','linewidth',2)
xlabel('\rho _{tor}')
ylabel('q')
% plot(sqrt(tor_flux_profile/max(tor_flux_profile)),q_initial_profile,'r')
plot(sqrt(tor_flux_profile/max(tor_flux_profile)),-1.0*finesse_data(1:257,end),'g')

%%
pol_P_reverse=pol_P(end:-1:1);

psi_scale=finesse_data(1:257,end-1);
psi_range=psi_scale/max(psi_scale);


Pinput_profile=polyval(pol_P_reverse,psi_range)*pcore;

figure(2)
pcore=PTOWB_profile_interp_ini(1)
hold on
grid on
plot(rho_lin_scale,PTOWB_profile_interp_ini,'b','linewidth',2)
% plot(sqrt(tor_flux_profile/max(tor_flux_profile)),P_initial_profile,'r')
% plot(rho_lin_scale,PTOWB_profile_interp_end,'r','linewidth',2)
plot(sqrt(tor_flux_profile/max(tor_flux_profile)),Pinput_profile,'r')
% plot(1-psi_scale/max(psi_scale),PTOWB_profile_interp_ini/PTOWB_profile_interp_ini(1),'b','linewidth',2)

xlabel('\rho _{tor}')
ylabel('P NBI')

