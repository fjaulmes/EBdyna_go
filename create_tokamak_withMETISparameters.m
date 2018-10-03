% needs to be taken as an input parameter
format compact

disp('----------------------------------------------------------------------------------------')
disp('The following parameters are critical for the')
disp('proper representation of your tokamak. Please make sure they have been correctly defined ')
disp('BEFORE going any further!')
disp('----------------------------------------------------------------------------------------')
load('METIS_profiles.mat')
SIGN_TOROIDAL_FIELD=sign(mean(finesse_data(1:end,end-13)))

a=METISgeo.amr;
B0=SIGN_TOROIDAL_FIELD*METISgeo.B0;

disp('Polynomial from FINESSE.inp')

F2poly_coefs=[0.0    -0.152    -0.004     0.196    -0.104    -0.004];F2oly_coefs=F2poly_coefs(end:-1:1)
Ppoly_coefs=[ 1.0   -0.91   1.515  -3.95   2.75 -0.5 0.1	];Ppoly_coefs=Ppoly_coefs(end:-1:1)

pol_F2=F2poly_coefs;
pol_P =Ppoly_coefs;
  
   
load finesse_data.mat
disp('FINESSE has dimensionless (X,Z) data')
finesse_data(:,1)=a*finesse_data(:,1);
finesse_data(:,2)=a*finesse_data(:,2);

%obvisouly there is often a small mismatch between 
%the theoretical value of a and the one given
%by the solver but normalization were done according to
%the theoretical one, so we need to cope with this small inconsistence
a_recalc=0.5*(max(finesse_data(:,1))+abs(min(finesse_data(:,1))));
b_recalc=0.5*(max(finesse_data(:,2))+abs(min(finesse_data(:,2))));
elongation=b_recalc/a_recalc

elongation_METIS=METISgeo.kappa

% filename='finesse20140123_155923.917_+0100.dat'

% messy file with no clear delimiter
try
    first_line = dlmread(filename,' ',[0 0 0 18]);
catch
    first_line = dlmread(filename,' ',[0 0 0 16]);
end
finesse_tokamak_data=first_line(find(first_line~=0));

% tokamak parameters from the code
epsilon=finesse_tokamak_data(1);
alpha_finesse=finesse_tokamak_data(2);
gamma_thermo=finesse_tokamak_data(3);
X_axis=a*finesse_tokamak_data(4);
Z_axis=a*finesse_tokamak_data(5);

plasma_beta_tot=finesse_tokamak_data(8)
plasma_beta_pol=finesse_tokamak_data(9)
Diamagnetism_on_axis=finesse_tokamak_data(7)

Bphi0=B0*Diamagnetism_on_axis

save('finesse_tokamak_parameters.mat','SIGN_TOROIDAL_FIELD','a','elongation','B0','Bphi0','epsilon','alpha_finesse','X_axis','Z_axis',...
    'Diamagnetism_on_axis','plasma_beta_tot','plasma_beta_pol','pol_F2','pol_P');