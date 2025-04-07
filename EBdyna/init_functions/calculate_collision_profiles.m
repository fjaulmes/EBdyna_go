% lambdaD=sqrt(epsilon0.*Te_prof./((eV).*ne_prof));
% log_lambda=log(12.*pi.*ne_prof.*lambdaD.^3);
% log_lambda=log(12*pi*sqrt((epsilon0.*Te_prof).^3./(ne_prof.*eV^3)))
% 
% log_lambda2=23-log(2*sqrt(ne_prof*1e-6).*Te_prof.^-1.5);
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


%%
% log_lambda_ie_2=15.2-0.5*log(ne_prof*1e-20)+log(Te_prof.*1e-3);
% log_lambda_ii_2=17.3-0.5*log(ne_prof*1e-20)+1.5*log(0.5*(Te_prof+Ti_prof).*1e-3);

% figure
% plot(log_lambda_ie);hold on;plot(log_lambda_ii)
% legend('log lambda ie','log lambda ii')


% nuNBie2=(1/(3*(pi)^1.5)).*(eV^4/epsilon0^2)./(me*mD).*(ne_prof)./((vthe_prof).^3).*log_lambda;
nuNBie=(1/(3*(2*pi)^1.5)).*(sqrt(me)*eV^4/epsilon0^2)./(mD*eV^1.5).*(ne_prof)./((Te_prof).^1.5).*log_lambda_ie;
nuNBii_v3=(1/(4*pi)).*(eV^4/epsilon0^2/(mD^2)).*(ni_prof).*log_lambda_ii;

% now collision profiles without log(Lambda)
% nuNBie=(1/(3*(2*pi)^1.5)).*(sqrt(me)*eV^4/epsilon0^2)./(mD*eV^1.5).*(ne_prof)./((Te_prof).^1.5);
% nuNBii_v3=(1/(4*pi)).*(eV^4/epsilon0^2/(mD^2)).*(ni_prof);

Ecrit=Te_prof.*(mD./me).^(1/3).*(3*sqrt(pi)/4).^(2/3);

taus=6.32*1e14*(mD/mH)*(Te_prof.^1.5)./(ne_prof.*log_lambda_ie);
taus=1./nuNBie;

INJECTION_ENERGY=80*1e3
taus2=(1./(1.5.*nuNBie)).*log(1+(INJECTION_ENERGY./Ecrit).^1.5);
taus3=(1./(1.5.*nuNBie)).*log(1+(0.5*INJECTION_ENERGY./Ecrit).^1.5);

taus=taus3;

THERMALIZATION_TIME=mean(taus3)

% vthe_prof=sqrt(2*eV*Te_prof/me);
