volume_flux_diff=volume_flux*0;
volume_flux_diff(2:Nradial)=volume_flux(2:end)-volume_flux(1:end-1);
Btot_PR_map=sqrt(Bpol_PR_map.^2+Btor_PR_map.^2);
Btot_radial_profile=mean(Btot_PR_map(1:NP-1,:),1);
beta_e_profile=2*mu0*Ne_profile.*Te_profile./(Btot_radial_profile.^2)

vA_radial_profile=Btot_radial_profile./(sqrt(mu0*Ne_profile));
vthe_radial_profile=sqrt(2*1.5*Te_profile/me);
lambda_e=vA_radial_profile./vthe_radial_profile;
gm_lambdae=(0.5*sqrt(pi))*(1+2*lambda_e.^2+2*lambda_e.^4).*exp(-lambda_e.^2);
gm_lambdae_3=(0.5*sqrt(pi))*(1+2*(lambda_e/3).^2+2*(lambda_e/3).^4).*exp(-(lambda_e/3).^2);
pos_tae=round(interp1(q_initial_profile,1:Nradial,10.5/6))
GT_lambdae=gm_lambdae(pos_tae)+gm_lambdae(pos_tae)
GT_lambdae_radial_profile=gm_lambdae+gm_lambdae;
kpll=-(1./(2*R0*q_initial_profile));
gamma_e=kpll.*vA_radial_profile.*beta_e_profile.*(q_initial_profile.^2).*GT_lambdae_radial_profile;

figure(1)
set(gca,'fontsize',20)
grid on; 
hold on;
plot(q_initial_profile,gamma_e)
% xlim([1 1.5])
xlabel('q')
ylabel('\gamma_e')

plasma_beta_recalc=mean(2*mu0*volume_flux_diff.*P_initial_profile)/mean(volume_flux_diff.*Btot_radial_profile.^2)