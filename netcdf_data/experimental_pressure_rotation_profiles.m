% Test writing and reading NetCDF4 with chunking and compression using
% Native Matlab routines (tested in Matlab 2010b)
% Rich Signell (rsignell@usgs.gov)
load('../data_tokamak/physics_constants.mat')
load('../data_tokamak/flux_geometry.mat')
%load('../data_tokamak/pressure_profile.mat')
rho_tor_scale=(tor_flux_profile/max(tor_flux_profile));
psi_norm_scale=1-psi_scale/max(psi_scale);
rho_pol_scale=1-sqrt(psi_scale/max(psi_scale));

psi_bar_rho=interp1(rho_tor_scale,psi_norm_scale,(0:256)/256);
rho_pol_rho_tor=interp1(rho_tor_scale,rho_pol_scale,(0:256)/256);
rho_lin_scale=(0:256)/256;

filename='30382A07.CDF'
% % if 0
%     % read data from OpenDAP using NJ Toolbox (http://njtbx.sourceforge.net)
%     url='http://coast-enviro.er.usgs.gov/thredds/dodsC/bathy/etopo1_bed_g2';
%     nc=mDataset(url);
%     topo=nc{'topo'}(1:12:end,1:12:end);
%     g=nc{'topo'}(1:12:end,1:12:end).grid;
%     topo(topo<0)=0;
%     lon=g.lon;
%     lat=g.lat;
%     save topo.mat topo lon lat
% % end
% % or load previously save mat file
% load topo.mat
% [ny,nx]=size(topo);
% subplot(211);pcolor(lon,lat,double(topo));shading flat;caxis([0 5000])
% title('Topo from Mat file');
% %%
% % write NetCDF4 with chunking & compression (deflation)
% 
% ncid = netcdf.create(filename,'NETCDF4');
% latdimid = netcdf.defDim(ncid,'lat',ny);
% londimid = netcdf.defDim(ncid,'lon',nx);
% varid = netcdf.defVar(ncid,'topo','short',[latdimid londimid]);
% lonid = netcdf.defVar(ncid,'lon','float',[londimid]);
% latid = netcdf.defVar(ncid,'lat','float',[latdimid]);
% 
% netcdf.defVarChunking(ncid,varid,'CHUNKED',[180 360]);
% netcdf.defVarDeflate(ncid,varid,true,true,5);
% netcdf.putAtt(ncid,latid,'units','degrees_north');
% netcdf.putAtt(ncid,lonid,'units','degrees_east');
% netcdf.putAtt(ncid,varid,'units','m');
% %netcdf.putAtt(ncid,varid,'missing_value',int16(-32767));
% 
% netcdf.putVar(ncid,lonid,[0],[nx],lon(1:nx));
% netcdf.putVar(ncid,latid,[0],[ny],lat(1:ny));
% netcdf.putVar(ncid,varid,[0 0],[ny nx],topo(1:ny,1:nx));
% 
% netcdf.close(ncid);
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


netcdf.close(ncid);
title('Topo from NetCDF4 file');

CRASH_TIME_STAMP=117

TE_profile_interp_ini=interp1((0:39)/39,TE_profile(:,CRASH_TIME_STAMP),rho_lin_scale);
TI_profile_interp_ini=interp1((0:39)/39,TI_profile(:,CRASH_TIME_STAMP),rho_lin_scale);
NE_profile_interp_ini=interp1((0:39)/39,NE_profile(:,CRASH_TIME_STAMP),rho_lin_scale)*1e6;
NI_profile_interp_ini=interp1((0:39)/39,ND_profile(:,CRASH_TIME_STAMP),rho_lin_scale)*1e6;
OMEGA_profile_interp_ini=interp1((0:39)/39,OMEGA_profile(:,CRASH_TIME_STAMP),rho_lin_scale);
Q_profile_interp_ini=interp1((0:39)/39,Q_profile(:,CRASH_TIME_STAMP),rho_lin_scale);
PTOWB_profile_interp_ini=interp1((0:39)/39,PTOWB_profile(:,CRASH_TIME_STAMP),rho_lin_scale);

TE_profile_interp_end=interp1((0:39)/39,TE_profile(:,CRASH_TIME_STAMP+1),rho_lin_scale);
TI_profile_interp_end=interp1((0:39)/39,TI_profile(:,CRASH_TIME_STAMP+1),rho_lin_scale);
NE_profile_interp_end=interp1((0:39)/39,NE_profile(:,CRASH_TIME_STAMP+1),rho_lin_scale)*1e6;
NI_profile_interp_end=interp1((0:39)/39,ND_profile(:,CRASH_TIME_STAMP+1),rho_lin_scale)*1e6;
OMEGA_profile_interp_end=interp1((0:39)/39,OMEGA_profile(:,CRASH_TIME_STAMP+1),rho_lin_scale);
Q_profile_interp_end=interp1((0:39)/39,Q_profile(:,CRASH_TIME_STAMP+1),rho_lin_scale);
PTOWB_profile_interp_end=interp1((0:39)/39,PTOWB_profile(:,CRASH_TIME_STAMP+1),rho_lin_scale);

N_profile_interp_ini=NE_profile_interp_ini+NI_profile_interp_ini;
T_profile_interp_ini=TE_profile_interp_ini+TI_profile_interp_ini;
N_profile_interp_end=NE_profile_interp_end+NI_profile_interp_end;
T_profile_interp_end=TE_profile_interp_end+TI_profile_interp_end;

PBULK_profile_interp_ini=(NE_profile_interp_ini.*TE_profile_interp_ini)+(NI_profile_interp_ini.*TI_profile_interp_ini);
PBULK_profile_interp_end=(NE_profile_interp_end.*TE_profile_interp_end)+(NI_profile_interp_end.*TI_profile_interp_end);

%%
figure(3)
hold on
grid on
plot(rho_lin_scale,Q_profile_interp_ini,'b','linewidth',2)
plot(rho_lin_scale,Q_profile_interp_end,'r','linewidth',2)
xlabel('\psi _{pol}')
ylabel('q')

%%

figure(1)
set(gca,'fontsize',20)
pcore=PTOWB_profile_interp_ini(1)
hold on
grid on
plot(rho_lin_scale,PTOWB_profile_interp_ini,'b','linewidth',3)
plot(rho_lin_scale,PTOWB_profile_interp_end,'r','linewidth',3)
legend('before sawtooth','after sawtooth')
xlabel('\rho _{tor}')
ylabel('Ptot (Pa)')

%%

figure(2)
set(gca,'fontsize',20)

hold on
grid on
CORE_ROTATION=OMEGA_profile_interp_ini(1)
plot(rho_lin_scale,(OMEGA_profile_interp_ini),'b','linewidth',3)
plot(rho_lin_scale,(OMEGA_profile_interp_end),'r','linewidth',3)
legend('before sawtooth','after sawtooth')
xlabel('\rho _{tor}')
ylabel('\Omega (rad/s)')
% 
