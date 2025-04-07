
% defining an electron temperature on axis
% this is arbitrary and should be tuned according to experimental value
% Te0=(4.4*1e3)
FILENAME=strcat(FIESTA_FOLDER,FIESTA_FILENAME);
try
    load(FILENAME)
catch
    run([FIESTA_FOLDER 'extract_FIESTA_CU.m']);
    load(FILENAME)
end

% flip_sign_Bt_Ip_FIESTA

try
    FILENAME=strcat(METIS_FOLDER,METIS_FILENAME);
    load(FILENAME)
catch
    error('Could not read METIS input data!')
    return
end
TMETIS=1.55
Psi_METIS=interp1(post.profil0d.temps,post.profil0d.psi,TMETIS);
DPSI_METIS=Psi_METIS(end)-Psi_METIS(1)

disp('Updating flux values in FIESTA data to 2pi Wb')
FIESTA.psi_boundary=FIESTA.psi_boundary/(2*pi);
FIESTA.psi_axis=FIESTA.psi_axis/(2*pi);
FIESTA.Psi=FIESTA.Psi/(2*pi);
FIESTA.params.psia=FIESTA.params.psia/(2*pi);
FIESTA.prof.psi=FIESTA.prof.psi/(2*pi);

FIESTA.psi_Xpoint=interp2(FIESTA.scale_R,FIESTA.scale_Z,FIESTA.Psi',FIESTA.R_xpoint,FIESTA.Z_xpoint)

aspect_ratio=FIESTA.params.aspectratio;

a=FIESTA.params.r0_geom/aspect_ratio
R0=FIESTA.params.r0_geom
epsilon=1/aspect_ratio;

Raxis=FIESTA.params.r0;
X_axis=Raxis-R0
Z_axis=FIESTA.params.z0

DPSI=FIESTA.psi_boundary-FIESTA.psi_axis
psi_global=DPSI;

disp('Toroidal sign convention of direct <> >0 incompatible with (R,Z,phi) tokamak coordinates')
disp('Fortunately Btor>0 always in FIESTA (also in clockwise configuration!)')

Bphi0=FIESTA.params.b0
Fdia_METIS=interp1(post.profil0d.temps,post.profil0d.fdia,TIME_EQUIL);
B0_METIS=Fdia_METIS(end)/0.894
Bvac0=sign(Bphi0).*B0_METIS
FVAC=Bvac0*0.894

SIGN_TOROIDAL_FIELD=sign(Bphi0);

SIGN_CURRENT=sign(DPSI) % current clockwise from top means psi more and more negative (DPSI>0)
SIGN_CO_CURRENT_FIELD=SIGN_CURRENT*SIGN_TOROIDAL_FIELD;
% special case where we want to have negative current 
% with respect to field (to keep Z upward)
if (SIGN_CO_CURRENT_FIELD<0) && (SIGN_TOROIDAL_FIELD>0)
   ASDEX_LIKE_EQUILIBRIUM=1
end

plasma_beta_tor=FIESTA.params.betat
plasma_beta_pol=FIESTA.params.betap

plasma_beta_tot=1/(1/plasma_beta_tor+1/plasma_beta_pol);


FILENAME=strcat(DATA_FOLDER,'tokamak_parameters.mat')
save (FILENAME,'SIGN_CO_CURRENT_FIELD','plasma_beta_pol','plasma_beta_tor','aspect_ratio','DPSI','Bphi0','Bvac0');
