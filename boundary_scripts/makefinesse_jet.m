%27.1.2012: Creates input files for finesse from CRONOS data and HELENA
%equilibrium

%t=50s is at index 161 for JET(6) 2000 and 2001

%SO CALLED INPUT PARAMETERS
clear all 
close all
load('jetdata_fabien_85383.mat')


ind=32; %number of moments in boundary decomposition

B0=3.3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear i x y theta rho
tind=565; 
tDpre=565; 
tepre=142; 
tmpre=594; 

boundR=double(squeeze(Rbnd(:,tind)));
boundZ=double(squeeze(Zbnd(:,tind)));

%r0=data.geo.r0(tind);
%z0=data.geo.z0(tind);
a=(max(boundR)-min(boundR))/2


r0=(max(boundR)+min(boundR))/2
z0=(max(boundZ)+min(boundZ))/2

B0=Bphi(tmpre)

% csou=(data.prof.te(tind,:)*1.6e-19/2/1.67e-27).^0.5;

x=boundR-r0;% x=x/max(x);
y=boundZ-z0;% y=y/max(x);

% a = 0.5*(max(max(data.geo.R))-min(min(data.geo.R)))
fac=1/a; 
% fac=1/a0; %.0583; 
% fac=1;

[theta,rho]=cart2pol(x,y);

%figure; plot(x,y,'.',0,0,'*tind'); axis equal; hold on
%polar(theta,rho)

%messy fourier transform
%theta(end)=[];
%rho(end)=[];
%rho1=interp1([theta(end),theta(1)],[rho(end),rho(1)],0);

% theta(1)=0;
theta(theta<0)=theta(theta<0)+2*pi;
[theta,Itheta] = sort(theta);
rho=rho(Itheta);

[theta_extra,Itheta1,Itheta2] = unique(theta);
rho_extra=rho(Itheta1);

theta_extra(1)=0;
theta_extra(end)=2*pi;
%theta=[theta' 2*pi];
%rho=[rho' rho1];

t=0:0.001:2*pi;
r=interp1(theta_extra,rho_extra,t); r=double(r);

clear r_test a0 b c cmin

c0=2*(1/pi*trapz(t,r,2)/2)*fac;
r_test=0.5*c0;
for m=1:ind
	acoef(m)=1/pi*trapz(t,r.*cos(m.*t),2);
	b(m)=1/pi*trapz(t,r.*sin(m.*t),2);
	c(m)=2*1/2*(acoef(m)-i*b(m))*fac;
	cmin(m)=1/2*(acoef(m)+i*b(m));
	%r_test=r_test+c(m).*exp(i.*t.*m)+cmin(m).*exp(-i.*t.*m);
	r_test=r_test+real(c(m).*exp(i.*t.*m));
end

figure;
polar(t,r); hold on; polar(t,r_test/fac,'r'); 

psi1=1

% psinorm=(data.prof.psi(tind,:)-data.prof.psi(tind,1))/2/pi; 
% psi1=-psinorm(end);
% psinorm=psinorm./psinorm(end);

nD_prof=squeeze(nD(:,tDpre));

%p=(data.prof.pe(tind,:)+data.prof.pion(tind,:))*data.gene.wdia(tind)./data.gene.wth(tind);

te_prof=squeeze(Te(14:end-2,tepre));

% F=data.equi.F(tind,:);
% vtor=data.prof.vtor_exp(tind,:);
%Fprime=

% x=linspace(0,1,101);
% psi=data.prof.psi(tind,:)/2/pi;
% 
% dFdpsi=gradient(F,psi);

beta=mean(te_prof)*mean(nD_prof)/(B0^2/2/(4*pi*1e-7));
 
% grho=data.equi.grho(tind,:);
% rhomax=data.equi.rhomax(tind);

% q=data.prof.q(tind,:);
% 
% vtorfac=0;
% vtornew=vtor*vtorfac;
% vtor=vtor*vtorfac;
% 
% xrho=data.equi.rhoRZ(tind,:)./data.equi.rhoRZ(tind,101);
% Rout=interp1(xrho,data.equi.R(tind,:,1),x);
% 
% ne=data.prof.ne(tind,:);
% Bphi=data.equi.BPHI(tind,:,1);
% Bphix=interp1(xrho,Bphi,x);
% Ftest=Bphix.*Rout;

% Ffinesse=F.^2-F(1)^2; %Ffinesse=Ffinesse/Ffinesse(1);
% dFdpsi=gradient(Ffinesse,psinorm);


alpha = a^2 * B0 / psi1

% disp('Ffinesse coef = ');
% disp((eps*alpha^2)/(a^2 * B0^2));
% Ffinesse=(eps*alpha^2)/(a^2 * B0^2)*Ffinesse;

%%%WRITE THE DATA FILES OUT%%%


dlmwrite('.\boundary.dat',ind+1,' ');
dlmwrite('.\boundary.dat',[c0 0],'-append','delimiter',' ');
dlmwrite('.\boundary.dat',[real(c)' imag(c)'],'-append','delimiter',' ');


xnD=linspace(0,1,length(nD_prof));
xnD=xnD.^2;
figure(2)
plot(xnD,nD_prof/nD_prof(1))


x=linspace(0,1,length(te_prof));
x=x.^2;
figure(3)
plot(x,te_prof/te_prof(1))

nD_prof_recalc=interp1(xnD,nD_prof,x)';

% x=linspace(0,1,length(te_prof));
% x=x.^2;
figure(6);
grid on; hold on
% plot(x,nD_prof_recalc/(nD_prof_recalc(1)))

plot(x,te_prof.*nD_prof_recalc/(te_prof(1)*nD_prof_recalc(1)))

% dlmwrite('.\delta1.dat',101,' ');
% dlmwrite('.\delta1.dat',[psinorm' Ffinesse'],'-append','delimiter',' ');
% 
% dlmwrite('.\delta2.dat',101,' ');
% dlmwrite('.\delta2.dat',abs([psinorm' p'/p(1)]),'-append','delimiter',' ');
% 
% dlmwrite('.\delta3.dat',101,' ');
% dlmwrite('.\delta3.dat',abs([psinorm' ne'/ne(1)]),'-append','delimiter',' ');
% 
% dlmwrite('.\delta4.dat',101,' ');
% if vtor(1) > 0 
% 	dlmwrite('.\delta4.dat',abs([psinorm' vtor'/vtor(1)]),'-append','delimiter',' ');
% else
% 	dlmwrite('.\delta4.dat',abs([psinorm' vtor']),'-append','delimiter',' ');
% end
% 	
% dlmwrite('.\basicdata.dat',[b0;r0;z0;psi1],'delimiter',' ');
% dlmwrite('.\profiles.dat',[psinorm',(F.^2)',(F.*dFdpsi)',p',vtor',ti',q']);
% dlmwrite('.\profiles.dat',[psinorm',(F.^2)',(F.*dFdpsi)',p',vtor',ti',q']);

