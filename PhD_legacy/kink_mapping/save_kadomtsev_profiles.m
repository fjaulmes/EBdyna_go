%load('psi_profiles.mat', 'psi_pol_initial_profile')

psi_star_initial=SIGN_PSIH*Psih;
psi_star_final=SIGN_PSIH*Psih_final;
% psi_star_initial=Psih;
% psi_star_final=Psih_final;

% to discriminate between kink and reconnection
% for Epot calculation
ISKINK=1

FILENAME=strcat(DATA_FOLDER,'psi_profiles_kadomtsev.mat')
save (FILENAME,'ISKINK','psi_star_initial','Psih_values','psi_star_final','psi_global','psi_star_max','psi_pol_initial_profile','xPsih_zero','q_initial','q_final','xPsih_max','psi_star_max');
