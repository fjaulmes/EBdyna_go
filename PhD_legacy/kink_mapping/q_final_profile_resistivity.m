reset_data_analysis_environment

Te1=interp1(1:257,Te_profile/kb,psi_rank_q1)
Ne1=interp1(1:257,Ne_profile,psi_rank_q1)
lambda_coulomb=17;

omega_ce=eV*Bavg/me;

tau_e=(3.44*1e5)*(Te1^1.5)/(Ne1*lambda_coulomb)
sigma0=Ne1*(eV^2)*tau_e/me

sigma_perp=1.96*sigma0;
% sigma_perp=sigma0*(omega_ce*tau_e)^2;

eta=1/sigma_perp

coef_diff=eta/mu0