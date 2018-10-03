load('../data_tokamak/B_fields.mat', 'Bphi0')
ksi_dot=gradient(ksi0_evol_lin,time_scale_lin*(4*1e-4));
rx_dot=gradient(rx_evol_lin,time_scale_lin*(4*1e-4));
u_sep=ksi_dot-rx_dot;
r_core=rx_evol_lin-ksi0_evol_lin;
f1=u_sep./ksi_dot;
f2=rx_dot./ksi_dot;

[minval TRANSITION_FRAME]=min(r_core);
TRANSITION_FRAME=length(Bstar1_evol_lin)



%%
C0=sqrt(1/(epsilon0*mu0));
Ne1=interp1(1:Nradial,Ne_profile,psi_rank_q1)
Te1=interp1(1:Nradial,Te_profile,psi_rank_q1)
vthe1=sqrt(2*Te1/me)
omegape=sqrt((Ne1*eV^2)/(epsilon0*me))

dvolume1_evol_lin=gradient(volume1_evol_lin,1:length(volume1_evol_lin));
dvolume2_evol_lin=gradient(volume2_evol_lin,1:length(volume1_evol_lin));
dvolume3_evol_lin=gradient(volume3_evol_lin,1:length(volume1_evol_lin));
dvolume1_evol_lin=abs(dvolume1_evol_lin);

E1=0.5*Bstar1_evol_lin.^2/mu0;
E2=0.5*(dvolume2_evol_lin./dvolume1_evol_lin).*Bstar2_evol_lin.^2/mu0;
Efinal=0.5*Bstar3_evol_lin.^2/mu0;

% deltaE=(r_core(1:length(E1))/r_value_q1_mean).*(E1);
deltaE=(f1(1:length(E1)).*E1+f2(1:length(E1)).*E2)-Efinal;
% deltaE=(E1+E2)-Efinal;
vA_out=sqrt(2*mu0*deltaE)/sqrt(mu0*Ne1*mD);
% l_out=(Bphi0./sqrt(2*mu0*deltaE))*r_value_q1_mean;
l_out=R0/(1-q_initial_profile(1));

tau_star=r_value_q1_mean./vA_out;

delta_evol=((C0/omegape)*(1./tau_star+(vthe1./l_out)).^(1/2)).*(tau_star.^(1/2));
% figure(2)
% hold on
% plot(time_scale_lin(1:length(vA_out)),delta_evol./tau_star)
% plot(time_scale_lin(1:length(vA_out)),u_sep(1:length(vA_out)),'r')

delta_avg=mean(delta_evol(8:end-2));

Delta_SP=2*delta_evol.*vA_out./u_sep(1:length(vA_out));
% Delta_SP=2*delta_avg*vA_out./u_sep_avg;
% plot(time_scale_lin(3:length(vA_out)),Delta_SP(3:length(vA_out)))

%%
% representation of the circles
DX_lin=0.0001;
Xlinscale=(0:DX_lin:1);
acoef=zeros(length(Xlinscale),1);

deltaarray=zeros(length(Xlinscale),1);
delta_posx=zeros(length(r_core),1);
delta_acoef=zeros(length(r_core),1);
delta_theta0=zeros(length(r_core),1);

linearray=zeros(1,length(Xlinscale));
Ycore=zeros(length(Xlinscale),length(r_core));
Yrx=zeros(length(Xlinscale),length(r_core));
for t=1:length(r_core)
    Ycore(:,t)=sqrt(r_core(t)^2-Xlinscale.^2)-r_core(t);
    Yrx(:,t)=sqrt(rx_evol_lin(t)^2-Xlinscale.^2)-rx_evol_lin(t);
end
Yrx(imag(Yrx(:,:))~=0)=0;
Ycore(imag(Ycore(:,:))~=0)=0;
% Yrx(~isreal(Yrx))=0;
% Ycore(~isreal(Ycore))=0;
for t=2:TRANSITION_FRAME
    acoef(2:end-1)=(Yrx(2:end-1,t)+rx_evol_lin(t))./Xlinscale(2:end-1)';
    XI1=zeros(length(Xlinscale),1);
    YI1=zeros(length(Xlinscale),1);
    XI2=zeros(length(Xlinscale),1);
    YI2=zeros(length(Xlinscale),1);
    length_core=floor(r_core(t)/DX_lin);
    for x=2:length_core
        if acoef(x)~=0
            linearray=acoef(x)*Xlinscale-rx_evol_lin(t);
            linearray=linearray';
            XI1(x)=interp1((Ycore(1:length_core,t)-linearray(1:length_core)),Xlinscale(1:length_core),0,'cubic');
            YI1(x)=interp1(Xlinscale,Ycore(:,t),XI1(x));
            XI2(x)=interp1((Yrx(1:length_core,t)-linearray(1:length_core)),Xlinscale(1:length_core),0,'cubic');
            YI2(x)=interp1(Xlinscale,Yrx(:,t),XI2(x));
        end
    end
    deltaarray=sqrt((YI2-YI1).^2+(XI2-XI1).^2);
    [valmin deltaindex]=min(abs(delta_avg-deltaarray(:)));
    delta_posx(t)=deltaindex;
    delta_acoef(t)=acoef(deltaindex);
    delta_theta0(t)=pi/2-atan(delta_acoef(t));
end


delta_theta0=delta_theta0';

%%
vA_out_normal=cos(0.5*delta_theta0(1:length(vA_out))).*vA_out;

tau_star=(delta_theta0(1:length(vA_out)).*rx_evol_lin(1:length(vA_out)))./vA_out_normal;

delta_evol=0.65*((C0/omegape)*(1./tau_star+(vthe1./l_out)).^(1/2)).*(tau_star.^(1/2));

delta_avg=mean(delta_evol(8:end-2));
% u_sep_avg=mean(u_sep);

Delta_SP=2*delta_evol.*vA_out_normal./ksi_dot(1:length(vA_out));



%%
figure(2)
hold on
% plot(time_scale_lin(1:length(vA_out)),2*delta_evol.*vA_out./(2*delta_theta0(1:length(vA_out)).*rx_evol_lin(1:length(vA_out))))
% plot(time_scale_lin(1:length(vA_out)),ksi_dot(1:length(vA_out)),'r')

DELTA_TIME=(time_scale_lin(2)-time_scale_lin(1))*4*1e-4;
ksi_dot_recalc=2*delta_evol.*vA_out_normal./(2*delta_theta0(1:length(vA_out)).*rx_evol_lin(1:length(vA_out)));
ksi_dot_recalc(1)=0;
ksi_dot_recalc(TRANSITION_FRAME+1)=0.5*(ksi_dot_recalc(TRANSITION_FRAME));
ksi_dot_recalc=0.98*ksi_dot_recalc;

ksi_dot_recalc(TRANSITION_FRAME:length(time_scale_lin))=0;
ksi_recalc_evol=zeros(length(time_scale_lin),1);
for t=2:length(time_scale_lin)
    ksi_recalc_evol(t)=ksi_recalc_evol(t-1)+ksi_dot_recalc(t-1)*DELTA_TIME;
end

plot(time_scale_lin,ksi_recalc_evol)
plot(time_scale_lin,ksi0_evol_lin,'r')
