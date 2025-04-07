function [METIS]=extract_METIS_profiles_EBdyna_extended(METISFILE,OUTPUTFILE,time)

load('FIESTA_equil.mat')

if nargin()==0
METISFILE='./CU_4_4_dsol.mat'
end
if nargin()<=1
OUTPUTFILE='METIS_data_CU'
time=1.55
end
SAVEFILE=1

% 
load(METISFILE)
format compact
rext         = (post.profil0d.Raxe + post.z0dinput.geo.a * post.profil0d.xli);
bulk_rot     =  post.profil0d.omega .* 	rext;

load('./data_tokamak/physics_constants.mat')
load('./data_tokamak/motions_map_dimensions.mat')
load('./data_tokamak/flux_geometry.mat','psi_global')

NB_DATA=19;

t_scale=post.profil0d.temps;

output_array=zeros(NB_DATA,length(post.profil0d.xli));

output_array(1,:)=interp1(t_scale,post.profil0d.psi(:,:),time);
output_array(2,:)=interp1(t_scale,post.profil0d.df2dpsi(:,:),time);
output_array(3,:)=interp1(t_scale,post.profil0d.dptotdpsi(:,:),time);
output_array(4,:)=interp1(t_scale,post.profil0d.ptot(:,:),time);
output_array(5,:)=interp1(t_scale,post.profil0d.fdia(:,:),time);
output_array(6,:)=interp1(t_scale,post.profil0d.jli(:,:),time);
output_array(7,:)=interp1(t_scale,post.profil0d.nep(:,:),time);
output_array(8,:)=interp1(t_scale,post.profil0d.tep(:,:),time);
output_array(9,:)=interp1(t_scale,post.profil0d.tip(:,:),time);
output_array(10,:)=interp1(t_scale,post.profil0d.er(:,:),time);
output_array(11,:)=interp1(t_scale,post.profil0d.pnbi(:,:),time);
output_array(12,:)=interp1(t_scale,post.profil0d.pnbi_ion(:,:),time);
output_array(13,:)=interp1(t_scale,post.profil0d.bpol(:,:),time);
output_array(14,:)=interp1(t_scale,post.profil0d.vpr(:,:),time);
output_array(15,:)=interp1(t_scale,post.profil0d.phi(:,:),time);
output_array(16,:)=interp1(t_scale,post.profil0d.vtor(:,:),time);
output_array(17,:)=interp1(t_scale,bulk_rot(:,:),time);
output_array(18,:)=interp1(t_scale,post.profil0d.n0m(:,:),time);
output_array(19,:)=interp1(t_scale,post.profil0d.qjli(:,:),time);
output_array(20,:)=interp1(t_scale,post.profil0d.nip(:,:),time);
output_array(21,:)=interp1(t_scale,post.profil0d.Raxe(:,:),time);
output_array(22,:)=interp1(t_scale,post.profil0d.kx(:,:),time);
output_array(23,:)=interp1(t_scale,post.profil0d.dx(:,:),time);

% output_array(2,:)=post.profil0d.nep(TIME_CU,:)*1e-20;
% output_array(3,:)=post.profil0d.tep(TIME_CU,:)*1e-3;
% output_array(4,:)=post.profil0d.bpol(TIME_CU,:);


t_scale2=post.z0dinput.cons.temps;
t_scale3=post.zerod.temps;

Ip=interp1(t_scale2,post.z0dinput.cons.ip(:),time);
R=interp1(t_scale2,post.z0dinput.geo.R(:),time);
a=interp1(t_scale2,post.z0dinput.geo.a(:),time);
K=interp1(t_scale2,post.z0dinput.geo.K(:),time);
d=interp1(t_scale2,post.z0dinput.geo.d(:),time);
z0=interp1(t_scale2,post.z0dinput.geo.z0(:),time);
b0=interp1(t_scale2,post.z0dinput.geo.b0(:),time);

% flux=interp1(t_scale3,flux_cs,time);
d0=interp1(t_scale3,post.zerod.d0(:),time);
vplasma=interp1(t_scale3,post.zerod.vp(:),time);
vloop=interp1(t_scale3,post.zerod.vloop(:),time);
teped=interp1(t_scale3,post.zerod.teped(:),time);
neped=interp1(t_scale3,post.zerod.neped(:),time);
te0=interp1(t_scale3,post.zerod.te0(:),time);
ne0=interp1(t_scale3,post.zerod.ne0(:),time);
nem=interp1(t_scale3,post.zerod.nem(:),time);
dsol=interp1(t_scale3,post.zerod.dsol(:),time);
ndd=interp1(t_scale3,post.zerod.ndd(:),time);
ndd_th=interp1(t_scale3,post.zerod.ndd_th(:),time);
ndd_nbi_th=interp1(t_scale3,post.zerod.ndd_nbi_th(:),time);

betan=interp1(t_scale3,post.zerod.betan(:),time);

METIS=struct();

METIS.time=time;
METIS.Ip=Ip;
METIS.Raxis=R+d0;
METIS.a=a;
METIS.K=K;
METIS.triangularity=d;
METIS.z0=z0;
METIS.vplasma=vplasma;
METIS.vloop=vloop;
METIS.d0=d0;
METIS.teped=teped;
METIS.neped=neped;
METIS.te0=te0;
METIS.ne0=ne0;
METIS.b0=b0;
METIS.nem=nem;
METIS.dsol=dsol;
METIS.ndd=ndd;

METIS.prof.psi=squeeze(output_array(1,:));
METIS.prof.df2dpsi=squeeze(output_array(2,:));
METIS.prof.dptotdpsi=squeeze(output_array(3,:));
METIS.prof.fdia=squeeze(output_array(5,:));
METIS.prof.ptot=squeeze(output_array(4,:));
METIS.prof.jli=squeeze(output_array(6,:));
METIS.prof.ne=squeeze(output_array(7,:));
METIS.prof.ni=squeeze(output_array(20,:));
METIS.prof.Te=squeeze(output_array(8,:));
METIS.prof.Ti=squeeze(output_array(9,:));
METIS.prof.Rmp=linspace(0,1,21)*(a-d0)+METIS.Raxis;
METIS.prof.Er=squeeze(output_array(10,:));
METIS.prof.pnbi=squeeze(output_array(11,:));
METIS.prof.pnbi_ion=squeeze(output_array(12,:));
METIS.prof.bpol=squeeze(output_array(13,:));
METIS.prof.vprprof=squeeze(output_array(14,:));
METIS.prof.phitor=squeeze(output_array(15,:));
METIS.prof.vtor=squeeze(output_array(16,:));
METIS.prof.bulk_rot=squeeze(output_array(17,:));
METIS.prof.n0m=squeeze(output_array(18,:));
METIS.prof.qp=squeeze(output_array(19,:));
% #20 next to ne
METIS.prof.Raxis=squeeze(output_array(21,:));
METIS.prof.Kx=squeeze(output_array(22,:));
METIS.prof.dx=squeeze(output_array(23,:));


n0_prof=METIS.prof.n0m;

ne_prof=METIS.prof.ne;
ni_prof=METIS.prof.ni;

Te_prof=METIS.prof.Te;
Ti_prof=METIS.prof.Ti;
vtor_prof=METIS.prof.vtor;


xprof=squeeze(post.profil0d.xli);
% bug in METIS beta
beta_plasma=2*mu0*trapz(xprof,METIS.prof.vprprof.*METIS.prof.ptot)/vplasma/(b0^2+(trapz(xprof,METIS.prof.vprprof.*METIS.prof.bpol)/vplasma)^2);
METIS.beta_plasma=beta_plasma;

vA=sqrt(METIS.prof.bpol.^2+b0^2)./sqrt(mu0.*mD.*METIS.prof.ne);

INJECTION_ENERGY=80*1e3
ZEFF=1.3;

calculate_collision_profiles

METIS.prof.Ecrit=Ecrit;
METIS.prof.taus=taus;
% METIS.prof.tauS2=taus2;
% METIS.prof.tauS3=taus3;
METIS.prof.vA=vA;
METIS.prof.nuNBie=nuNBie;
METIS.prof.nuNBii_v3=nuNBii_v3;

psi_norm=1-(METIS.prof.psi-METIS.prof.psi(end))./(METIS.prof.psi(1)-METIS.prof.psi(end));
METIS.prof.psi_norm=psi_norm;
METIS.prof.vthi_prof=vthi_prof;
METIS.prof.vthe_prof=vthe_prof;

if SAVEFILE==1
    save(OUTPUTFILE,'METIS');
end

% now extending the profiles into SOL

% EXP_LINSPACE=linspace(0.5*0.025,0.2,81);
DELTA_PSIN=0.1*0.025
EXP_LINSPACE=DELTA_PSIN:DELTA_PSIN:0.25

psi_norm = [psi_norm 1+EXP_LINSPACE];

scale_psiN_mp=(scale_psi_mp-scale_psi_mp(1))/psi_global;


% Plasma Phys. Control. Fusion 57 (2015) 125011
qcyl=FIESTA.params.b0./FIESTA.params.aspectratio./METIS.prof.bpol(end)

%lambda_T =  2.63  *  FIESTA.params.b0^-0.5  *   FIESTA.params.q95^0.974  *  6^0.05 *1e-3
% a bit innacurate about qcyl!
lambda_T =  2.63  *  FIESTA.params.b0^-0.5  *   qcyl^0.974  *  6^0.05 *1e-3

Rsep=interp1(scale_psiN_mp,scale_R_mp,1)
psi_lambdaT=interp1(scale_R_mp,scale_psiN_mp,Rsep+lambda_T);
delta_psi_lambdaT=psi_lambdaT-1
% plot(scale_R_mp,scale_psiN_mp)


% about 3.2 mm decay length in temperature 
% lambda_T =  2.63  *  5.0^-0.5  *   2.5^0.974  *  6^0.05
% (https://iopscience.iop.org/article/10.1088/0741-3335/57/12/125011/pdf)
% (Plasma Phys. Control. Fusion 57 (2015) 12501)
% 3.2mm here and 
% 10 mm distance in physical space for this extra normalized flux of 0.12
% exp(-0.1/0.1*0.01/0.0032) = 0.0067

a_minor_radius=FIESTA.params.r0/FIESTA.params.aspectratio;

% EXP_LINSPACE=a_minor_radius*EXP_LINSPACE*1/(lambda_T);
EXP_LINSPACE_R=interp1(scale_psiN_mp,scale_R_mp,1+EXP_LINSPACE)-Rsep
EXP_LINSPACE_R=EXP_LINSPACE_R/(lambda_T);

FAC_LN_LT=1/1.35

ne_sep=ne_prof(end)
Te_sep=Te_prof(end)

n0_prof  = [n0_prof EXP_LINSPACE_R*0+n0_prof(end)];
ne_prof  = [ne_prof exp(-FAC_LN_LT*EXP_LINSPACE_R)*ne_prof(end)];
ni_prof  = [ni_prof exp(-FAC_LN_LT*EXP_LINSPACE_R)*ni_prof(end)];
Te_prof  = [Te_prof exp(-EXP_LINSPACE_R)*Te_prof(end)];
Ti_prof  = [Ti_prof exp(-EXP_LINSPACE_R)*Ti_prof(end)];
vtor_prof  = [vtor_prof exp(-EXP_LINSPACE_R)*vtor_prof(end)];

WALL_TEMP=0.06

Te_prof=max(Te_prof,WALL_TEMP);
Ti_prof=max(Ti_prof,WALL_TEMP);


calculate_collision_profiles

psi_norm


nuNBie=(1/(3*(2*pi)^1.5)).*(sqrt(me)*eV^4/epsilon0^2)./(mD*eV^1.5).*(ne_prof)./((Te_prof).^1.5).*log_lambda_ie;

Ecrit=Te_prof.*(mD./me).^(1/3).*(3*sqrt(pi)/4).^(2/3);

INJECTION_ENERGY=80*1e3
taus_80keV=(1./(1.5.*nuNBie)).*log(1+(0.5*INJECTION_ENERGY./Ecrit).^1.5);

figure(3)
hold on
plot(psi_norm,taus_80keV);
legend('taus 80keV')


%%
rho_pol=sqrt(psi_norm);
roh_pol_mp=sqrt(scale_psiN_mp);
ne_mp_prof=interp1(rho_pol,ne_prof,roh_pol_mp);
Te_mp_prof=interp1(rho_pol,Te_prof,roh_pol_mp);


wall_pos=1.2

figure;
hold on
grid on
plot(scale_R_mp,ne_mp_prof*1e-20)
plot(scale_R_mp,Te_mp_prof*1e-3)

plot([wall_pos wall_pos],[0 5],'k--')
plot([Rsep Rsep],[0 5],'b--')

xlabel('R [m]')

legend('ne','Te')


Rvalue=1.15:0.00005:1.17;
% Ln=LT
plot(Rvalue,1e-20*ne_sep*exp((Rsep-Rvalue)/((1/FAC_LN_LT)*lambda_T)),'m--')

%%
if SAVEFILE==1
        save('./data_tokamak/pressure_profile.mat','-append','psi_norm','n0_prof','ne_prof','ni_prof','Te_prof','Ti_prof','vtor_prof','nuNBie','nuNBii_v3','vthi_prof','vthe_prof','taus_80keV')
end
