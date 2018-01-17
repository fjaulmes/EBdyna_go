
clc
format compact

METIS_filename='COMPASS_15606';
METIS_filename=strcat(METIS_filename,'.mat')

time_METIS=1.12

SAVENAME='METIS_profiles.mat'

load(METIS_filename);post.profil0d.temps;

%%
B0_mat=repmat(post.z0dinput.geo.b0,1,size(post.profil0d.rmx,2));
rho_tor_metis=post.profil0d.rmx.^2;
rho_tor_metis=rho_tor_metis.*B0_mat*pi;
rho_tor_metis_max=max(rho_tor_metis');
rho_tor_metis_max=repmat(rho_tor_metis_max,size(post.profil0d.rmx,2),1);
rho_tor_metis_max=rho_tor_metis_max';
rho_tor_metis=sqrt(rho_tor_metis./rho_tor_metis_max);

rho_pol_metis=post.profil0d.psi;
rho_pol_metis_max=max(rho_pol_metis');
rho_pol_metis_max=repmat(rho_pol_metis_max,size(post.profil0d.rmx,2),1);
rho_pol_metis_max=rho_pol_metis_max';
psi_metis=post.profil0d.psi-rho_pol_metis_max;
rho_pol_metis_min=min(psi_metis');
rho_pol_metis_min=repmat(rho_pol_metis_min,size(post.profil0d.rmx,2),1);
rho_pol_metis_min=rho_pol_metis_min';
rho_pol_metis=psi_metis./abs(rho_pol_metis_min);
rho_pol_metis=min(rho_pol_metis,0);
rho_pol_metis=max(rho_pol_metis,-1);
psi_norm_metis=-rho_pol_metis;
rho_pol_metis=sqrt(-rho_pol_metis);

R0_mat=repmat(post.z0dinput.geo.R,1,size(post.profil0d.rmx,2));


METISdata=struct;

METISdata.time_METIS=time_METIS;
METISdata.rho_tor_metis=get_profile_at_time(rho_tor_metis,time_METIS,post.profil0d.temps);
METISdata.psi_metis=get_profile_at_time(post.profil0d.psi,time_METIS,post.profil0d.temps);
METISdata.psi_norm_metis=get_profile_at_time(psi_norm_metis,time_METIS,post.profil0d.temps);
METISdata.rho_pol_metis=get_profile_at_time(rho_pol_metis,time_METIS,post.profil0d.temps);
METISdata.q_profile=get_profile_at_time(post.profil0d.qjli,time_METIS,post.profil0d.temps);
METISdata.Te_profile=get_profile_at_time(post.profil0d.tep,time_METIS,post.profil0d.temps);
METISdata.Ti_profile=get_profile_at_time(post.profil0d.tip,time_METIS,post.profil0d.temps);
METISdata.ne_profile=get_profile_at_time(post.profil0d.nep,time_METIS,post.profil0d.temps);
METISdata.ni_profile=get_profile_at_time(post.profil0d.nip,time_METIS,post.profil0d.temps);
METISdata.Omega_profile=get_profile_at_time(post.profil0d.vtor./(R0_mat),time_METIS,post.profil0d.temps);
METISdata.Ptot_profile=get_profile_at_time(post.profil0d.ptot,time_METIS,post.profil0d.temps);
METISdata.dPtot_profile=get_profile_at_time(post.profil0d.dptotdpsi,time_METIS,post.profil0d.temps);
METISdata.nW_profile=get_profile_at_time(post.profil0d.nwp,time_METIS,post.profil0d.temps);
METISdata.pnbi_profile=get_profile_at_time(post.profil0d.pnbi,time_METIS,post.profil0d.temps);
METISdata.pnbi_ion_profile=get_profile_at_time(post.profil0d.pnbi_ion,time_METIS,post.profil0d.temps);
METISdata.Fdia_profile=get_profile_at_time(post.profil0d.fdia,time_METIS,post.profil0d.temps);
METISdata.jphi_profile=get_profile_at_time(post.profil0d.jeff,time_METIS,post.profil0d.temps);
METISdata.iRavg_profile=get_profile_at_time(post.profil0d.ri,time_METIS,post.profil0d.temps);

METISgeo=struct;
METISgeo.time_METIS=time_METIS;
METISgeo.Ip=get_value_at_time(post.z0dinput.cons.ip ,time_METIS,post.profil0d.temps);
METISgeo.amr=get_value_at_time(post.z0dinput.geo.a ,time_METIS,post.profil0d.temps);
METISgeo.kappa=get_value_at_time(post.z0dinput.geo.K ,time_METIS,post.profil0d.temps);
METISgeo.R0=get_value_at_time(post.z0dinput.geo.R ,time_METIS,post.profil0d.temps);
METISgeo.B0=get_value_at_time(post.z0dinput.geo.b0 ,time_METIS,post.profil0d.temps);
METISgeo.Raxis=get_profile_at_time(post.profil0d.Raxe ,time_METIS,post.profil0d.temps);
METISgeo.Raxis=METISgeo.Raxis(1);
METISgeo.Zaxis=get_value_at_time(post.z0dinput.geo.z0 ,time_METIS,post.profil0d.temps);
METISgeo.Rsep=get_profile_at_time(post.z0dinput.exp0d.Rsepa,time_METIS,post.profil0d.temps);
Zs = post.z0dinput.exp0d.Zsepa +  post.z0dinput.geo.z0 * ones(1,size(post.z0dinput.exp0d.Rsepa,2));
METISgeo.Zsep=get_profile_at_time(Zs,time_METIS,post.profil0d.temps);
METISgeo.Rsepc=get_profile_at_time(post.profil0d.Rsepa,time_METIS,post.profil0d.temps);
METISgeo.Zsepc=get_profile_at_time(post.profil0d.Zsepa,time_METIS,post.profil0d.temps);


save(SAVENAME,'METISgeo','METISdata','time_METIS');

disp('---------------------------------------------------------------------')
disp('METIS data ready to use in FINESSE and EBdyna_go !')
disp('---------------------------------------------------------------------')

