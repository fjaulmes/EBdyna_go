% Physics constants
kb=1.381e-23;
epsilon0=8.8542e-12;
mu0=4*pi*1e-7;
gamma_thermo=1.6667;
eV=1.6021766208*1e-19;

me=9.1094e-31;
u_mass=1.66054e-27;
mD=2.014*u_mass;
mT=3.016*u_mass;
mDT=(mD+mT)*0.5;
mHe=4.0026*u_mass;
mH=1.00794*u_mass;  

ZHe=2;

FILENAME=strcat(DATA_FOLDER,'physics_constants.mat')
save (FILENAME,'kb','epsilon0','mu0','gamma_thermo','eV','me','mD','mT','mDT','mHe','mH','ZHe');