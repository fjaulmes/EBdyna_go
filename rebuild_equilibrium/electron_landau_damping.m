volume_flux_diff=volume_flux*0;
volume_flux_diff(2:257)=volume_flux(2:end)-volume_flux(1:end-1);
Bpol_PR_map=sqrt(BX_PR_map.^2+BZ_PR_map.^2);
Btot_PR_map=sqrt(Bpol_PR_map.^2+Btor_PR_map.^2);
Btot_radial_profile=mean(Btot_PR_map(1:NP-1,:),1);

beta_e_profile_initial=2*mu0*Ne_profile.*Te_profile./(Btot_radial_profile.^2);
vA_radial_profile_initial=Btot_radial_profile./(sqrt(mu0*Ne_profile*mDT));
vthe_radial_profile=sqrt(2*1.5*Te_profile/me);
lambda_e=vA_radial_profile_initial./vthe_radial_profile;
gm_lambdae=(0.5*sqrt(pi))*lambda_e.*(1+2*lambda_e.^2+2*lambda_e.^4).*exp(-lambda_e.^2);
gm_lambdae_3=(0.5*sqrt(pi))*(lambda_e/3).*(1+2*(lambda_e/3).^2+2*(lambda_e/3).^4).*exp(-(lambda_e/3).^2);
% pos_tae=round(interp1(q_initial_profile,1:257,1.19));
% GT_lambdae=gm_lambdae(pos_tae)+gm_lambdae(pos_tae);
GT_lambdae_radial_profile=gm_lambdae+gm_lambdae_3;
kpll_ini=(1./(2*R0*q_initial_profile));
gamma_e_initial=-kpll_ini.*vA_radial_profile_initial.*beta_e_profile_initial.*(q_initial_profile.^2).*GT_lambdae_radial_profile;

% these final profiles come from kinetic calculations (outputs folder, kinetic Pprofile script)

beta_e_profile_final=2*mu0*Ne_final_profile.*Te_final_profile./(Btot_radial_profile.^2);
beta_e_profile_final(find(isnan(beta_e_profile_final)))=beta_e_profile_initial(find(isnan(beta_e_profile_final)));
vA_radial_profile_final=Btot_radial_profile./(sqrt(mu0*Ne_final_profile*mDT));
vA_radial_profile_final(find(isnan(vA_radial_profile_final)))=vA_radial_profile_initial(find(isnan(vA_radial_profile_final)));
vthe_radial_profile=sqrt(2*1.5*Te_final_profile/me);
lambda_e_final=vA_radial_profile_final./vthe_radial_profile;
gm_lambdae_final=(0.5*sqrt(pi))*lambda_e_final.*(1+2*lambda_e_final.^2+2*lambda_e_final.^4).*exp(-lambda_e_final.^2);
gm_lambdae_3_final=(0.5*sqrt(pi))*(lambda_e_final/3).*(1+2*(lambda_e_final/3).^2+2*(lambda_e_final/3).^4).*exp(-(lambda_e_final/3).^2);
% pos_tae=round(interp1(q_final_profile_diff,1:257,1.19));
% GT_lambdae=gm_lambdae(pos_tae)+gm_lambdae(pos_tae);
GT_lambdae_radial_profile_final=gm_lambdae_final+gm_lambdae_3_final;
kpll_final=(1./(2*R0*q_final_profile_diff));
gamma_e_final=-kpll_final.*vA_radial_profile_final.*beta_e_profile_final.*(q_final_profile_diff.^2).*GT_lambdae_radial_profile_final;
gamma_e_final(isnan(gamma_e_final))=gamma_e_initial(isnan(gamma_e_final));

figure(1)
set(gca,'fontsize',20)
grid on; 
hold on;
plot(q_initial_profile,gamma_e_initial,'b--','linewidth',3)
plot(q_final_profile_diff,gamma_e_final,'r','linewidth',3)
xlim([1 1.5])
xlabel('q')
ylabel('\gamma_e (s^{-1})')
legend('initial','after sawtooth')

plasma_beta_recalc=mean(2*mu0*volume_flux_diff.*P_initial_profile)/mean(volume_flux_diff.*Btot_radial_profile.^2);