
% defining an electron temperature on axis
% this is arbitrary and should be tuned according to experimental value
% Te0=(4.4*1e3)
FILENAME=strcat(FINESSE_FOLDER,'finesse_tokamak_parameters.mat');
load(FILENAME)


% alpha_finesse=alpha;
aspect_ratio=1/epsilon;
R0=aspect_ratio*a
psi0_input=(B0*a^2)/alpha_finesse
Baxis=Bphi0
SIGN_CO_CURRENT_FIELD=sign(psi0_input)
% special case where we want to have negative current 
% with respect to field (to keep Z upward)
if (SIGN_CO_CURRENT_FIELD<0) && (SIGN_TOROIDAL_FIELD>0)
   ASDEX_LIXE_EQUILIBRIUM=1
end

plasma_beta_tor=plasma_beta_tot
plasma_beta_tot=1/(1/plasma_beta_tor+1/plasma_beta_pol);


% shafranov_shift=(1/(epsilon))*(sqrt(1+epsilon^2)-1)*a
Raxis=R0+X_axis;

% Overall poloidal flux on magnetic axis
% Iaxis=R0*Bphi0/mu0;


FILENAME=strcat(DATA_FOLDER,'tokamak_parameters.mat')
save (FILENAME,'SIGN_CO_CURRENT_FIELD','plasma_beta_tot','plasma_beta_pol','plasma_beta_tor','aspect_ratio','Baxis');
