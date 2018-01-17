%27.1.2012: Creates input files for finesse from METIS data
clc
load('./data_tokamak/physics_constants.mat')

%INPUT PARAMETERS
ind=200; %number of moments in boundary decomposition
PLOT_BOUNDARY=1
PLOT_PROFILES=1

vtorfac=0; %to zero out the toroidal rotation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
    load('METIS_profiles.mat')
catch
    disp('No METIS data file avilable!')
end

clear index
% using the vertically centered data to ease Fourier decomposition
boundR=METISgeo.Rsepc;
boundZ=METISgeo.Zsepc;

r0=METISgeo.R0;
z0=METISgeo.Zaxis;

rmin=METISgeo.amr;

x=boundR-r0;
y=boundZ;

% parameters for alpha and epsilon inputs of FINESSE
a0 = 0.5*(max(boundR)-min(boundR));
B0 = METISgeo.B0;
psi1 = max(abs(METISdata.psi_metis))-min(abs(METISdata.psi_metis));
R0 = METISgeo.R0;

% fac=2*a0; 
fac=2/a0;

[theta,rho]=cart2pol(x,y);

if PLOT_BOUNDARY
    figure; 
    plot(x,y,'.',0,0,'b'); axis equal; hold on
    polar(theta,rho)
end

%Fourier transform

rho0=interp1([theta(end),theta(1)],[rho(end),rho(1)],0);
% theta(1)=0;
% rho(1)=rho0;

theta(theta<0)=theta(theta<0)+2*pi;

% theta=[0 theta ];
% rho=[rho0 rho ];
% 
Dtheta=0.00001;
t=0:Dtheta:2*pi;
r=interp1(theta,rho,t);

clear r_test a b c cmin

c0=1/pi*trapz(t,r,2)/2*fac;
r_test=c0;
for m=1:ind
	a(m)=1/pi*trapz(t,r.*cos(m.*t),2);
	b(m)=1/pi*trapz(t,r.*sin(m.*t),2);
	c(m)=1/2*(a(m)-i*b(m))*fac;
	cmin(m)=1/2*(a(m)+i*b(m));
	%r_test=r_test+c(m).*exp(i.*t.*m)+cmin(m).*exp(-i.*t.*m);
	r_test=r_test+real(2*c(m).*exp(i.*t.*m));
end

if PLOT_BOUNDARY
    polar(t,r_test/fac,'r');
%     legend('METIS data','polar','polar test')
end



%%%WRITE THE DATA FILES OUT%%%

dlmwrite('.\boundary.dat',ind+1,' ');
dlmwrite('.\boundary.dat',[c0 0],'-append','delimiter',' ');
dlmwrite('.\boundary.dat',[real(c)' imag(c)'],'-append','delimiter',' ');

%% Now dealing with other FINESSE input parameters

alpha = a0^2 * B0 / psi1
epsilon = a0 / R0

Omega_tor_profile=METISdata.Omega_profile;
Fdia_profile=METISdata.Fdia_profile;
jphi_profile=METISdata.jphi_profile;
% Ptot_profile=METISdata.Ptot_profile;
% the pressure will be calculated from the polynomials
% to avoid inconsistencies
Ptot_profile=METISdata.Ptot_profile;
dPtot_profile=METISdata.dPtot_profile;
psi_norm_profile=METISdata.psi_norm_metis;
psi_profile=METISdata.psi_metis;
iRavg_profile=METISdata.iRavg_profile;

F2finesse=Fdia_profile.^2-Fdia_profile(1)^2; 
% Pfinesse=Ptot_profile/max(Ptot_profile);
Ofinesse=Omega_tor_profile/max(Omega_tor_profile);

F2poly_coefs=polyfit(psi_norm_profile,F2finesse,4);
% Ppoly_coefs=polyfit(psi_norm_profile,Pfinesse,5);
Opoly_coefs=polyfit(psi_norm_profile,Omega_tor_profile,4);

F2recalc=Fdia_profile(1).^2+polyval(F2poly_coefs,psi_norm_profile);
dF2recalc=gradient(F2recalc,psi_profile);
% Frecalc=sqrt(F2recalc);
% dFrecalc=gradient(Frecalc,psi_profile);


dPfinesse=(1./iRavg_profile).*(jphi_profile-0.5*(iRavg_profile/mu0).*dF2recalc);

Ptot_integ=dPtot_profile*0;
for index=1:length(Ptot_integ)-1
    Ptot_integ(end-index)=Ptot_integ(end-index+1)+(psi_profile(end-index)-psi_profile(end-index+1))*dPtot_profile(end-index+1);
end

Pfinesse=dPfinesse*0;
for index=1:length(Pfinesse)-1
    Pfinesse(end-index)=Pfinesse(end-index+1)+(psi_profile(end-index)-psi_profile(end-index+1))*dPfinesse(end-index+1);
end
Ppoly_coefs=polyfit(psi_norm_profile,Pfinesse/max(Pfinesse),6);



%%
% F2poly_coefs=[0.0    -0.14   -0.0064   0.195  -0.1137 ];F2poly_coefs=F2poly_coefs(end:-1:1)
% Ppoly_coefs=[  1.0   -0.91   1.515  -3.9528   2.7498 -0.5 0.1];Ppoly_coefs=Ppoly_coefs(end:-1:1)
% 
% pol_F2=F2poly_coefs;
% pol_P =Ppoly_coefs;

if PLOT_PROFILES
    figure
    subplot(2,1,1)
    hold on ; grid on;
    plot(psi_norm_profile,Fdia_profile.^2);
    plot(psi_norm_profile,Fdia_profile(1).^2+polyval(F2poly_coefs,psi_norm_profile),'r');
    title('F2dia')
    subplot(2,1,2)
    hold on ; grid on;
    plot(psi_norm_profile,Ptot_profile);
    plot(psi_norm_profile,Ptot_profile(1)*polyval(Ppoly_coefs,psi_norm_profile),'r');
    title('Ptot')
    legend('exp','FINESSE')
end


disp('you can copy/paste this in FINESSE input file')
Fpoly_coefs_FINESSE=F2poly_coefs(end:-1:1)
Ppoly_coefs_FINESSE=Ppoly_coefs(end:-1:1)

save('METIS_finesse_data.mat','alpha','epsilon','F2poly_coefs','Ppoly_coefs','Opoly_coefs')


