% reset_data_analysis_environment;
read_core_data_txt_file;
Raxis=R0+X_axis

coreRini=mean(Rcore_evol(1:150));
coreZini=mean(Zcore_evol(1:150));

Xevol=Rcore_evol-coreRini+X_axis;
Zevol=Zcore_evol-coreZini+Z_axis;

psi_core_evol=interp2(scale_X,scale_Z,psi_norm_XZsmall_map',Xevol,Zevol);
psi_value_evol=interp1(1:Nradial,psi_scale,psi_core_evol);
r_value_evol=interp1(1:Nradial,radial_r_value_flux,psi_core_evol);
rhopol_recalc_core_evol=sqrt(1-psi_value_evol/psi_global)

hold on
plot(time_core_evol-2.51,rhopol_core_evol,'b')
plot(time_core_evol-2.51,rhopol_recalc_core_evol,'r--')

figure(2)
hold on
plot(time_core_evol-2.512,r_value_evol,'k')

R_SXR_evol=Xevol+R0;
Z_SXR_evol=Zevol;

save tomas_SXR_data.mat time_core_evol r_value_evol R_SXR_evol Z_SXR_evol