% Test writing and reading NetCDF4 with chunking and compression using
% Native Matlab routines (tested in Matlab 2010b)
% Rich Signell (rsignell@usgs.gov)
load('../data_tokamak/physics_constants.mat')
load('../data_tokamak/flux_geometry.mat')
load('../data_tokamak/pressure_profile.mat')
load('../data_tokamak/q_profile.mat')

rho_tor_scale=sqrt(tor_flux_profile/max(tor_flux_profile));
psi_norm_scale=1-(psi_scale/max(psi_scale));
rho_pol_scale=1-sqrt(psi_scale/max(psi_scale));

psi_bar_rho=interp1(rho_tor_scale,psi_norm_scale,(0:256)/256);
rho_pol_rho_tor=interp1(rho_tor_scale,rho_pol_scale,(0:256)/256);
rho_lin_scale=(0:256)/256;
psi_bar_rho_lin=interp1(rho_lin_scale,psi_norm_scale,(0:256)/256);

filename='30382P01.CDF'

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


netcdf.close(ncid);
title('Topo from NetCDF4 file');

CRASH_TIME_STAMP=108

TE_profile_interp_ini=interp1((0:39)/39,TE_profile(:,CRASH_TIME_STAMP),rho_lin_scale);
TI_profile_interp_ini=interp1((0:39)/39,TI_profile(:,CRASH_TIME_STAMP),rho_lin_scale);
NE_profile_interp_ini=interp1((0:39)/39,NE_profile(:,CRASH_TIME_STAMP),rho_lin_scale)*1e6;
NI_profile_interp_ini=interp1((0:39)/39,ND_profile(:,CRASH_TIME_STAMP),rho_lin_scale)*1e6;
OMEGA_profile_interp_ini=interp1((0:39)/39,OMEGA_profile(:,CRASH_TIME_STAMP),rho_lin_scale);
Q_profile_interp_ini=interp1((0:39)/39,Q_profile(:,CRASH_TIME_STAMP),rho_lin_scale);
PTOWB_profile_interp_ini=interp1((0:39)/39,PTOWB_profile(:,CRASH_TIME_STAMP),rho_lin_scale);
NIMP_profile_interp_ini=interp1((0:39)/39,NIMP_profile(:,CRASH_TIME_STAMP),rho_lin_scale)*1e6;

TE_profile_interp_end=interp1((0:39)/39,TE_profile(:,CRASH_TIME_STAMP+1),rho_lin_scale);
TI_profile_interp_end=interp1((0:39)/39,TI_profile(:,CRASH_TIME_STAMP+1),rho_lin_scale);
NE_profile_interp_end=interp1((0:39)/39,NE_profile(:,CRASH_TIME_STAMP+1),rho_lin_scale)*1e6;
NI_profile_interp_end=interp1((0:39)/39,ND_profile(:,CRASH_TIME_STAMP+1),rho_lin_scale)*1e6;
OMEGA_profile_interp_end=interp1((0:39)/39,OMEGA_profile(:,CRASH_TIME_STAMP+1),rho_lin_scale);
Q_profile_interp_end=interp1((0:39)/39,Q_profile(:,CRASH_TIME_STAMP+1),rho_lin_scale);
PTOWB_profile_interp_end=interp1((0:39)/39,PTOWB_profile(:,CRASH_TIME_STAMP+1),rho_lin_scale);
NIMP_profile_interp_end=interp1((0:39)/39,NIMP_profile(:,CRASH_TIME_STAMP+1),rho_lin_scale)*1e6;

N_profile_interp_ini=NE_profile_interp_ini+NI_profile_interp_ini;
T_profile_interp_ini=TE_profile_interp_ini+TI_profile_interp_ini;
N_profile_interp_end=NE_profile_interp_end+NI_profile_interp_end;
T_profile_interp_end=TE_profile_interp_end+TI_profile_interp_end;

PBULK_profile_interp_ini=(NE_profile_interp_ini.*TE_profile_interp_ini)+(NI_profile_interp_ini.*TI_profile_interp_ini);
PBULK_profile_interp_end=(NE_profile_interp_end.*TE_profile_interp_end)+(NI_profile_interp_end.*TI_profile_interp_end);

% load('../data_tokamak/pressure_profile.mat')
% save('../data_tokamak/pressure_profile.mat','-append','rho_tor_scale','rho_pol_scale','TE_profile_interp_ini','TI_profile_interp_ini','NE_profile_interp_ini','NI_profile_interp_ini')
PTOT_profile_interp_ini=PTOWB_profile_interp_ini;
% save('../data_tokamak/pressure_profile.mat','-append','rho_tor_scale','rho_pol_scale','psi_bar_rho','TE_profile_interp_ini','TI_profile_interp_ini','NE_profile_interp_ini','NI_profile_interp_ini','PTOT_profile_interp_ini')
% save('../data_tokamak/pressure_profile.mat','-append','rho_tor_scale','rho_pol_scale','psi_bar_rho',...
%     'TE_profile_interp_ini','TI_profile_interp_ini','NE_profile_interp_ini','NI_profile_interp_ini',...
%     'OMEGA_profile_interp_ini','PTOT_profile_interp_ini','NIMP_profile_interp_ini');


TE_profile_radial_ini=interp1((0:Nradial-1)/(Nradial-1),TE_profile_interp_ini,rho_tor_scale);
TI_profile_radial_ini=interp1((0:Nradial-1)/(Nradial-1),TI_profile_interp_ini,rho_tor_scale);
NE_profile_radial_ini=interp1((0:Nradial-1)/(Nradial-1),NE_profile_interp_ini,rho_tor_scale);
NI_profile_radial_ini=interp1((0:Nradial-1)/(Nradial-1),NI_profile_interp_ini,rho_tor_scale);
PTOT_profile_radial_ini=interp1((0:Nradial-1)/(Nradial-1),PTOT_profile_interp_ini,rho_tor_scale);

for nt=1:size(PTOWB_profile,2)
    PTOT_profile_interp(:,nt)=interp1((0:39)/39,PTOWB_profile(:,nt),rho_lin_scale);
    PTOT_profile_radial(:,nt)=interp1((0:Nradial-1)/(Nradial-1),PTOT_profile_interp(:,nt),rho_tor_scale);
end


psi_mix=size_r-3
psiq1=interp1(Q_profile_interp_ini,1:Nradial,1)
rho_tor_mix=interp1(1:Nradial,rho_tor_scale,psi_mix)
rho_tor_q1=interp1(1:Nradial,rho_lin_scale,psiq1)
rmix=interp1(1:Nradial,radial_r_value_flux,psi_mix)
rho_mix_lin_pos=interp1(rho_lin_scale,1:Nradial,rho_tor_mix)

%%
figure(1)
subplot(2,1,1)
set(gca,'fontsize',16)
hold on
grid on
plot(rho_lin_scale,Q_profile_interp_ini,'b','linewidth',2)
plot(rho_lin_scale,Q_profile_interp_end,'r--','linewidth',2)

plot([rho_tor_q1 rho_tor_q1],[0 4],'y--','linewidth',2)
plot([rho_tor_mix rho_tor_mix],[0 4],'g--','linewidth',2)
xlabel('\rho _{tor}')
ylabel('q')
xlim([0 0.6])
ylim([0.5 2])

hl=legend('before sawtooth','after sawtooth','r_1','r_{mix}')
set(hl,'fontsize',14)


%%

figure(1)
subplot(2,1,2)
set(gca,'fontsize',16)

pcore=PTOWB_profile_interp_ini(1)
% Prho_ini=interp1(rho_lin_scale,PTOWB_profile_interp_ini,rho_tor_scale);

hold on
grid on
plot(rho_lin_scale,PTOWB_profile_interp_ini,'b','linewidth',2)
plot(rho_lin_scale,PTOWB_profile_interp_end,'r--','linewidth',2)
% xlabel('\psi _{pol}')
xlabel('\rho _{tor}')
ylabel('Ptot (Pa)')
plot([rho_tor_q1 rho_tor_q1],[0 1e5],'y--','linewidth',2)
plot([rho_tor_mix rho_tor_mix],[0 1e5],'g--','linewidth',2)

xlim([0 0.6])
ylim([0.5 8]*1e4)


hl=legend('before sawtooth','after sawtooth','r_1','r_{mix}')
set(hl,'fontsize',14)



%%
% 
% figure(1)
% subplot(3,1,3)
% set(gca,'fontsize',16)
% 
% hold on
% grid on
% CORE_ROTATION=OMEGA_profile_interp_ini(1)
% plot(rho_lin_scale,(OMEGA_profile_interp_ini),'b','linewidth',2)
% plot(rho_lin_scale,(OMEGA_profile_interp_end),'r--','linewidth',2)
% % xlabel('\psi _{pol}')
% xlabel('\rho _{tor}')
% ylabel('\Omega (rad/s)')
% plot([rho_tor_mix rho_tor_mix],[0 2e5],'g--','linewidth',2)
% 
% xlim([0 0.6])



%%
figure(2)

psi_mix=size_r-3
psi_rank_q1

subplot(2,1,1)
set(gca,'fontsize',16)
hold on
grid on

plot(TIME_profile,PTOT_profile_radial(1,:)*1e6,'b','linewidth',2)
xlim([2 3.0])
ylabel('P_{tot} (r=0) (Pa)')
xlabel('t (s)')

subplot(2,1,2)
set(gca,'fontsize',16)
hold on
grid on

plot(TIME_profile,PTOT_profile_radial(psi_rank_q1,:),'y--','linewidth',2)
plot(TIME_profile,PTOT_profile_radial(psi_mix,:),'r','linewidth',2)
xlabel('t (s)')
ylabel('P_{tot} (Pa)')
xlim([2 3.0])

legend('r=r_1','r=r_{mix}')