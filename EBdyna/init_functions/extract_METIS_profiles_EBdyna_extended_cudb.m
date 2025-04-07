function extract_METIS_profiles_EBdyna_extended_cudb(SHOT_NUMBER,time)

SAVEFILE=1
wall_pos=1.2

load('./data_tokamak/physics_constants.mat')


if nargin()==0
error('give a shot number!')
end

if nargin()<=1
warning('time point of simulation should be specified')
time=1.55
end
cdb=cdb_client();
% 
signal = cdb.get_signal(['Te/OMP_PROFILES:' num2str(SHOT_NUMBER)]);
Te_prof=interp1(signal.time_axis.data,signal.data,time);
signal = cdb.get_signal(['Ti/OMP_PROFILES:' num2str(SHOT_NUMBER)]);
Ti_prof=interp1(signal.time_axis.data,signal.data,time);
signal = cdb.get_signal(['ne/OMP_PROFILES:' num2str(SHOT_NUMBER)]);
ne_prof=interp1(signal.time_axis.data,signal.data,time);
signal = cdb.get_signal(['ni/OMP_PROFILES:' num2str(SHOT_NUMBER)]);
ni_prof=interp1(signal.time_axis.data,signal.data,time);
signal = cdb.get_signal(['vtor/OMP_PROFILES:' num2str(SHOT_NUMBER)]);
vtor_prof=interp1(signal.time_axis.data,signal.data,time);
% signal = cdb.get_signal(['n0/OMP_PROFILES:' num2str(SHOT_NUMBER)]);
% n0_prof=interp1(signal.time_axis.data,signal.data,time);
% signal = cdb.get_signal(['n0_D2/OMP_PROFILES:' num2str(SHOT_NUMBER)]);
% n0_D2_prof=interp1(signal.time_axis.data,signal.data,time);


vthe_prof=sqrt(2*eV*Te_prof/me);
vthi_prof=sqrt(2*eV*Ti_prof/mD);

%%
VFAST_80keV=sqrt(2*80*1e3*(eV/mD));
VFAST_40keV=sqrt(2*40*1e3*(eV/mD));

lambda_De = sqrt((epsilon0.*Te_prof)./(ne_prof.*eV));
b90_D=(eV^2)./(4*pi*epsilon0*0.5*mD.*(VFAST_40keV).^2);
b90_e=(eV^2)./(4*pi*epsilon0*0.5*me.*(max(vthe_prof,VFAST_40keV)).^2);
b90_D_nov=(eV^2)./(4*pi*epsilon0*0.5*mD);
% b90_e_nov=(eV^2)./(4*pi*epsilon0*0.5*me);

log_lambda_ie=log(lambda_De./b90_e);
log_lambda_ii=log(lambda_De./b90_D);
% log_lambda_ie=log(lambda_De./b90_e);
log_lambda_ii_nov=log(lambda_De./b90_D_nov);

nuNBie=(1/(3*(2*pi)^1.5)).*(sqrt(me)*eV^4/epsilon0^2)./(mD*eV^1.5).*(ne_prof)./((Te_prof).^1.5).*log_lambda_ie;
nuNBii_v3=(1/(4*pi)).*(eV^4/epsilon0^2/(mD^2)).*(ni_prof).*log_lambda_ii;

Ecrit=Te_prof.*(mD./me).^(1/3).*(3*sqrt(pi)/4).^(2/3);

% taus=6.32*1e14*(mD/mH)*(Te_prof.^1.5)./(ne_prof.*log_lambda_ie);
% taus=1./nuNBie;

% INJECTION_ENERGY=80*1e3
% taus2=(1./(1.5.*nuNBie)).*log(1+(INJECTION_ENERGY./Ecrit).^1.5);
% taus3=(1./(1.5.*nuNBie)).*log(1+(0.5*INJECTION_ENERGY./Ecrit).^1.5);



%%

INJECTION_ENERGY=80*1e3
taus_80keV=(1./(1.5.*nuNBie)).*log(1+(0.5*INJECTION_ENERGY./Ecrit).^1.5);
THERMALIZATION_TIME=mean(taus_80keV(1:170))

%%
figure
signal = cdb.get_signal(['R/OMP_PROFILES:' num2str(SHOT_NUMBER)]);
R_prof=interp1(signal.time_axis.data,signal.data,time);
Rsep=R_prof(signal.axis1.data==1)

ylim([0 2*taus_80keV(1)])
wall_pos_psi=interp1(R_prof,signal.axis1.data,wall_pos,'nearest')

hold;grid on;
plot(signal.axis1.data,taus_80keV);

plot([wall_pos_psi wall_pos_psi],[0 2*taus_80keV(1)],'k--')
plot([1 1],[0 2*taus_80keV(1)],'b--')
xlabel('\psi_N')
ylabel('\tau_s [s]')

%%

figure;
hold on
grid on
plot(R_prof,ne_prof*1e-20)
plot(R_prof,Te_prof*1e-3)

plot([wall_pos wall_pos],[0 max(Te_prof*1e-3)],'k--')
plot([Rsep Rsep],[0 max(Te_prof*1e-3)],'b--')

xlabel('R [m]')

legend('ne','Te')



%%
if SAVEFILE==1
        save('./data_tokamak/pressure_profile.mat','-append','nuNBie','nuNBii_v3','vthi_prof','vthe_prof','taus_80keV','vtor_prof')
end
