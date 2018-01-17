% needs to be taken as an input parameter
format compact

disp('----------------------------------------------------------------------------------------')
disp('The following parameters are critical for the')
disp('proper representation of your tokamak. Please make sure they have been correctly defined ')
disp('BEFORE going any further!')
disp('----------------------------------------------------------------------------------------')

a=0.92
Bphi0=2.95
Te0=(8.0*1e3)
disp('--------press a key to confirm these values---------------------------------------------')

pause

pol_F2=[0, -7.5, -3.3, 14.3, -7.7, 0.1 ]
pol_P=[1.0, -2.1, 2.16088709677, 2.26, -1, -0.12]

pol_N=[ 1. -0.52  -1.85  5.4 -3.7 ]
pol_T=[ 1. -0.8   8  -26.8  32.4 -13.5 ]

BETA_ALPHAS=0.0
   
   
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

% filename='finesse20140123_155923.917_+0100.dat'

% messy file with no clear delimiter
try
    first_line = dlmread(filename,' ',[0 0 0 18]);
catch
    first_line = dlmread(filename,' ',[0 0 0 17]);
end
finesse_tokamak_data=first_line(find(first_line~=0));

% tokamak parameters from the code
epsilon=finesse_tokamak_data(1);
alpha_finesse=finesse_tokamak_data(2);
gamma_thermo=finesse_tokamak_data(3);
X_axis=a*finesse_tokamak_data(4);
Z_axis=a*finesse_tokamak_data(5);

plasma_beta_tot=finesse_tokamak_data(8); 
plasma_beta_pol=finesse_tokamak_data(9);
Diamagnetism_on_axis=finesse_tokamak_data(7)

B0=Bphi0/Diamagnetism_on_axis

save('finesse_tokamak_parameters.mat','a','elongation','B0','Bphi0','epsilon','alpha_finesse','X_axis','Z_axis',...
    'Diamagnetism_on_axis','plasma_beta_tot','plasma_beta_pol','Te0','pol_F2','pol_P','pol_N','pol_T','BETA_ALPHAS');