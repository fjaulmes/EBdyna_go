volume_flux_diff=volume_flux*0;
volume_flux_diff(2:257)=volume_flux(2:end)-volume_flux(1:end-1);
Bpol_PR_map=sqrt(BX_PR_map.^2+BZ_PR_map.^2);
Btot_PR_map=sqrt(Bpol_PR_map.^2+Btor_PR_map.^2);
Btot_radial_profile=mean(Btot_PR_map(1:NP-1,:),1);

beta_e_profile_initial=2*mu0*Ne_profile.*Te_profile./(Btot_radial_profile.^2);
vA_radial_profile_initial=Btot_radial_profile./(sqrt(mu0*Ne_profile*mDT));


figure(1)
set(gca,'fontsize',20)
grid on; 
hold on;
% plot(q_final_profile_diff,0.5*(mD/eV)*(vA_radial_profile_initial).^2,'r','linewidth',3)
% plot(q_final_profile_diff,0.5*(mD/eV)*(vA_radial_profile_initial/3).^2,'b--','linewidth',3)
plot(q_final_profile_diff,vA_radial_profile_initial,'r','linewidth',3)
plot(q_final_profile_diff,vA_radial_profile_initial/3,'b--','linewidth',3)
xlim([1 2])
xlabel('q')
ylabel('v_A (m/s)')
legend('vA','vA/3')

plasma_beta_recalc=mean(2*mu0*volume_flux_diff.*P_initial_profile)/mean(volume_flux_diff.*Btot_radial_profile.^2);