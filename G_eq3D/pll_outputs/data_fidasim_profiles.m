run('../reset_data_analysis_environment.m');


rhot_exp_fs_ini =[       0.74120785      0.67150947      0.57139373      0.47823533      0.38951596 0.29249917      0.21234468      0.13344137     0.066461420];
exp_fs_ini=[      0.014930318     0.014680675     0.024663160     0.028218698     0.046563364 0.054612653     0.076098148     0.084483462      0.10622026];
error_fs_ini=[      0.0016077158   0.00067935158   0.00046578934   0.00068488629   0.00097915132 0.0017215187    0.0013678425    0.0012021697   0.00091266394];
     
rhot_EBd_fs_ini=[      0.0653424     0.131745     0.212243     0.292027     0.389061     0.477124 0.571105     0.671329     0.741412];
EBd_fs_ini=[     0.074968072     0.063044768     0.049270695     0.037751806     0.025860887 0.017054363     0.011572914    0.0060061099    0.0038760694];


%background offset correction
exp_fs_ini=exp_fs_ini-0.008;

figure(1)
set(gca,'fontsize',16)
hold on
grid on


errorbar(rhot_exp_fs_ini,exp_fs_ini,error_fs_ini,'b','LineWidth',3)
plot(rhot_EBd_fs_ini,EBd_fs_ini,'r--','LineWidth',3)

exp_fs_ini_interp=interp1(rhot_exp_fs_ini,exp_fs_ini,rho_tor_scale);
exp_fs_ini_interp(isnan(exp_fs_ini_interp))=0;
[max_exp_fs_ini_interp max_pos ]=max(exp_fs_ini_interp)
exp_fs_ini_interp(1:max_pos)=max_exp_fs_ini_interp;
[min_exp_fs_ini_interp min_pos ]=min(exp_fs_ini_interp(find(exp_fs_ini_interp>0)))
exp_fs_ini_interp(min_pos:end)=min_exp_fs_ini_interp.*(1-(0:Nradial-min_pos)/(Nradial-min_pos));


EBd_fs_ini_interp=interp1(rhot_EBd_fs_ini,EBd_fs_ini,rho_tor_scale);
EBd_fs_ini_interp(isnan(EBd_fs_ini_interp))=0;
[max_EBd_fs_ini_interp max_pos ]=max(EBd_fs_ini_interp)
EBd_fs_ini_interp(1:max_pos)=max_EBd_fs_ini_interp;
[min_EBd_fs_ini_interp min_pos ]=min(EBd_fs_ini_interp(find(EBd_fs_ini_interp>0)))
EBd_fs_ini_interp(min_pos:end)=min_EBd_fs_ini_interp.*(1-(0:Nradial-min_pos)/(Nradial-min_pos));

figure(2)
set(gca,'fontsize',16)
hold on
grid on


plot(radial_r_value_flux,exp_fs_ini_interp,'b','LineWidth',3)
plot(radial_r_value_flux,EBd_fs_ini_interp,'r--','LineWidth',3)

figure(3)
hold on
grid on
alphas_weight_radial_profile=exp_fs_ini_interp./EBd_fs_ini_interp;
for corr_edge_range=200:Nradial-1
    alphas_weight_radial_profile(corr_edge_range)=alphas_weight_radial_profile(corr_edge_range)+(corr_edge_range-200)^2*0.0002;
end
alphas_weight_radial_profile(end)=alphas_weight_radial_profile(end-1);

plot(radial_r_value_flux,alphas_weight_radial_profile,'b','LineWidth',3)

alphas_weight_radial_profile_poly=polyval(fit_weight.coeff,radial_r_value_flux);
plot(radial_r_value_flux,alphas_weight_radial_profile_poly,'b--','LineWidth',2)


save weight_correction_profile.mat radial_r_value_flux alphas_weight_radial_profile alphas_weight_radial_profile_poly