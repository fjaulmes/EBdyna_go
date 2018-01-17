%27.1.2012: Creates input files for finesse from CRONOS data and HELENA
%equilibrium

%t=50s is at index 161 for JET(6) 2000 and 2001

%SO CALLED INPUT PARAMETERS
ind=12; %number of moments in boundary decomposition

gammaE=0.15; 
%desired gammaE at x=0.33. This sets a multiplier on the entire rotation function which is then
%normalized to the mach number

%vtorfac=0.3/0.0493; %first number is desired gamE at x=0.33;

vtorfac=0; %to zero out the toroidal rotation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear i
tind=161;
boundR=data.geo.R(tind,:);
boundZ=data.geo.Z(tind,:);



r0=data.geo.r0(tind);
z0=data.geo.z0(tind);

%r0=(max(boundR)+min(boundR))/2;
%z0=(max(boundZ)+min(boundZ))/2;

b0=data.geo.b0(tind)
rmin=data.geo.a(tind);
csou=(data.prof.te(tind,:)*1.6e-19/2/1.67e-27).^0.5;


x=data.geo.R(tind,:)-r0;% x=x/max(x);
y=data.geo.Z(tind,:)-z0;% y=y/max(x);

%because of the error of normalization for coefficients in FINESSE
a0 = 0.5*(max(max(data.geo.R))-min(min(data.geo.R)))
fac=2*a0; 
%fac=1;

[theta,rho]=cart2pol(x,y);

%figure; plot(x,y,'.',0,0,'*tind'); axis equal; hold on
%polar(theta,rho)

%messy fourier transform
theta(1)=[];
rho(1)=[];
rho0=interp1([theta(end),theta(1)],[rho(end),rho(1)],0);

theta(theta<0)=theta(theta<0)+2*pi;

theta=[0 theta 2*pi];
rho=[rho0 rho rho0];

t=0:0.01:2*pi;
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

figure;
polar(t,r); hold on; polar(t,r_test/fac,'r'); 

psi1=abs(data.prof.psi(tind,end)-data.prof.psi(tind,1))
psinorm=data.prof.psi(tind,:)-data.prof.psi(tind,1);
psinorm=psinorm./psinorm(end);

p=data.prof.ptot(tind,:);
ti=data.prof.ti(tind,:);

F=data.equi.F(tind,:);
vtor=data.prof.vtor_exp(tind,:);
%Fprime=

x=linspace(0,1,101);
psi=data.prof.psi(tind,:);
psid1=data.prof.psid1(tind,:);
test=diff(F); test1=[test(1) test]; test2=[test test(end)];
dF=test1/2+test2/2; 

test=diff(psi); test1=[test(1) test]; test2=[test test(end)];
dpsi=test1/2+test2/2; 
dFdpsi=dF./dpsi; 



%figure; plot(psinorm,dFdpsi); hold on
%dFdpsi=smooth(dFdpsi',10);
%plot(psinorm,dFdpsi,'r');

beta=p./(b0^2/2/(4*pi*1e-7));
 
grho=data.equi.grho(tind,:);
rhomax=data.equi.rhomax(tind);
xedge=0.5;
Tiedge=ti(101);  Tiped=ti(xedge*100+1);

invgradTi=r0./data.prof.lti(tind,:)*fac;
 
tiprof=ti;
for k=1:xedge*100+1
tiprof(k)=Tiped*exp(rhomax./r0*trapz(k/100-0.01:0.01:xedge,1./grho(k:xedge*100+1).*invgradTi(k:xedge*100+1),2));
end

figure; plot(x,ti,x,tiprof,'--');

q=data.prof.q(tind,:);

vtornew=vtor*vtorfac;

xrho=data.equi.rhoRZ(tind,:)./data.equi.rhoRZ(tind,101);
Rout=interp1(xrho,data.equi.R(tind,:,1),x);

Rinner=(Rout-Rout(1)); Rinner=Rinner/Rinner(end);

test=diff(vtornew); test1=[test(1) test]; test2=[test test(end)];
dvtor=-(test1/2+test2/2); 

test=diff(Rout); test1=[test(1) test]; test2=[test test(end)];

dr=(test1/2+test2/2); 
dvtordr=dvtor./dr;

% gamE=dvtordr.*Rinner./q2004(tind,:)./ (csou./r0);

M=vtornew*r0./(sqrt(5/3)*csou);

% figure; plot(x,gamE); title('gamE');
figure; plot(x,M); title('M')

ne=data.prof.ne(tind,:);
Bphi=data.equi.BPHI(tind,:,1);
Bphix=interp1(xrho,Bphi,x);
Ftest=Bphix.*Rout;

Ffinesse=F.^2-F(1)^2; %Ffinesse=Ffinesse/Ffinesse(1);
%%%WRITE THE DATA FILES OUT%%%
test=diff(Ffinesse); test1=[test(1) test]; test2=[test test(end)];
test1/2+test2/2; 
test=diff(psinorm); test1=[test(1) test]; test2=[test test(end)];
dpsi=test1/2+test2/2; 
dFdpsi=dF./dpsi;


alpha = a0^2 * b0 / psi1
epsilon = a0 / r0

dlmwrite('.\boundary.dat',ind+1,' ');
dlmwrite('.\boundary.dat',[c0 0],'-append','delimiter',' ');
dlmwrite('.\boundary.dat',[real(c)' imag(c)'],'-append','delimiter',' ');

dlmwrite('.\delta1.dat',101,' ');
dlmwrite('.\delta1.dat',[psinorm' Ffinesse'],'-append','delimiter',' ');

dlmwrite('.\delta2.dat',101,' ');
dlmwrite('.\delta2.dat',abs([psinorm' p'/p(1)]),'-append','delimiter',' ');

dlmwrite('.\delta3.dat',101,' ');
dlmwrite('.\delta3.dat',abs([psinorm' ne'/ne(1)]),'-append','delimiter',' ');

dlmwrite('.\delta4.dat',101,' ');
dlmwrite('.\delta4.dat',abs([psinorm' vtor'/vtor(1)]),'-append','delimiter',' ');

dlmwrite('.\basicdata.dat',[b0;r0;psi1;z0]);
dlmwrite('.\profiles.dat',[psinorm',(F.^2)',(F.*dFdpsi)',p',vtor',tiprof']);

