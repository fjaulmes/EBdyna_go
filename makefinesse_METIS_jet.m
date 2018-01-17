%27.1.2012: Creates input files for finesse from METIS data
clc

%INPUT PARAMETERS
ind=256; %number of moments in boundary decomposition
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
psi1 = max(abs(METISdata.psi_metis));
R0 = METISgeo.R0;

% fac=2*a0; 
fac=2/a0;

[theta,rho]=cart2pol(x,y);

if PLOT_BOUNDARY
    figure; 
    plot(x,y,'.',0,0,'r'); axis equal; hold on
    polar(theta,rho)
end

%Fourier transform

rho0=interp1([theta(end),theta(1)],[rho(end),rho(1)],0);
theta(1)=[];
rho(1)=[];

theta(theta<0)=theta(theta<0)+2*pi;

theta=[0 theta ];
rho=[rho0 rho ];

t=0:0.001:2*pi;
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
    polar(t,r); hold on; polar(t,r_test/fac,'r');
    legend('METIS data','polar','polar test')
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
Ptot_profile=METISdata.Ptot_profile;
psi_norm_profile=METISdata.psi_norm_metis;

F2finesse=Fdia_profile.^2-Fdia_profile(1)^2; 
Pfinesse=Ptot_profile/max(Ptot_profile);
Ofinesse=Omega_tor_profile/max(Omega_tor_profile);

F2poly_coefs=polyfit(psi_norm_profile,F2finesse,5);
Ppoly_coefs=polyfit(psi_norm_profile,Pfinesse,5);
Opoly_coefs=polyfit(psi_norm_profile,Omega_tor_profile,5);

if PLOT_PROFILES
    figure
    subplot(2,1,1)
    hold on ; grid on;
    plot(psi_norm_profile,F2finesse);
    plot(psi_norm_profile,polyval(F2poly_coefs,psi_norm_profile));
    title('F2dia')
    subplot(2,1,2)
    hold on ; grid on;
    plot(psi_norm_profile,Pfinesse);
    plot(psi_norm_profile,polyval(Ppoly_coefs,psi_norm_profile));
    title('Ptot')
end

disp('you can copy/paste this in FINESSE input file')
Fpoly_coefs_FINESSE=F2poly_coefs(end:-1:1)
Ppoly_coefs_FINESSE=Ppoly_coefs(end:-1:1)

save('METIS_finesse_data.mat','alpha','epsilon','F2poly_coefs','Ppoly_coefs','Opoly_coefs')


