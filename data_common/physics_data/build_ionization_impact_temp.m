%CX against D at maxwellian temperature

SAVEFILE=0

isotope_mass=2;  % amu

Eb=[1 ...
5*10^2 ...
1*10^03  ...
2*10^03  ...
5*10^03 ...
1*10^04 ...
2*10^04 ...
5*10^04 ...
1*10^05 ...
2*10^05 ...
5*10^05 ...
1*10^06 ...
2*10^06 ...
5*10^06 ];


sigma_ionization_val=[0 ...
1.46*10^-20 ...
1.46*10^-19 ...
1.02*10^-18 ...
6.24*10^-18 ...
1.94*10^-17 ...
6.73*10^-17 ...
1.43*10^-16 ...
1.10*10^-16 ...
6.99*10^-17 ...
3.48*10^-17 ...
1.94*10^-17 ...
1.05*10^-17 ...
4.62*10^-18]*1e-4;  % (m^{2})


Eb_values_eV_simple=Eb*isotope_mass;


Eb_values_eV=linspace(0,2.1e5,1e3+1);

Eb_values_keV=1e-3*Eb_values_eV;

b1_io=2.016*1e-3;
b2_io=3.7154;
b3_io=3.989*1e-2;
b4_io=3.1413*1e-1;
b5_io=2.1254;
b6_io=6.399*1e3;
b7_io=6.1897*1e1;
b8_io=9.2731*1e3;

sigma_ionization_val2=1e-20.*b1_io.*((Eb_values_keV.^b2_io.*exp(-b3_io.*Eb_values_keV)./(1+b4_io.*Eb_values_keV.^b5_io))+b6_io.*exp(-b7_io./Eb_values_keV).*log(1+b8_io.*Eb_values_keV)./Eb_values_keV)
sigma_ionization_val2(1)=0;

Eb_values_eV=Eb_values_eV*isotope_mass;

%%
% figure(6)
loglog(Eb_values_eV_simple,sigma_ionization_val,'linewidth',2);

grid on
hold on

loglog(Eb_values_eV,sigma_ionization_val2,'linestyle','--','linewidth',2);
set(gca,'fontsize',16)

title('Ionization impact (D+) with non excitated states')
ylabel('\sigma [m^{2}]')
xlabel('Eb [eV]')
xlim([1 1e6])
ylim([10^-25 1.2*10^-19])

%%
figure(9)
vb=Eb*0;
if isotope_mass==2
    vb=sqrt(2*Eb_values_eV*eV/mD);
elseif isotope_mass==1
    vb=sqrt(2*Eb_values_eV*eV/mH);
end
rate_ionization=sigma_ionization_val2.*vb;

loglog(Eb_values_eV,rate_ionization,'k','linewidth',3);
grid on
hold on
set(gca,'fontsize',16)
ylabel('rate [m^3 s^{-1}]')
xlabel('Eb [eV]')

title('Ionization rate (D+) with non excitated states')
% ylim([10^-20 1.2*10^-15])




%%
% maxwellian distribution
% figure(9)

% T0_values_eV=linspace(1e3,10e3,11); % (background ions)
% T0_values_eV=[100 1e3 10e3]*2;
% T0_values_eV=[100 1e3 10*1e3];
T0_values_eV=[100 500 1*1e3 5e3 10e3];
T0_values=(eV)*T0_values_eV; % in eV (background neutrals)

sigma_ionization_T0=zeros(length(sigma_ionization_val2),length(T0_values)+1);
sigma_ionization_T0(:,1)=sigma_ionization_val2;

for ind_T0=1:length(T0_values)
    T0=T0_values(ind_T0);
    T0_ev=T0/eV
    
     v0=linspace(-3.8e6,3.8e6,4e2);
     vr=linspace(0,3.8e6,4e2);
%    v0=linspace(0,4.8e6,4e3);
    if isotope_mass==2
        fv0=(mD/(2*pi*T0))^0.5*exp(-mD.*v0.^2./(2*T0));
        fvr=(mD/(2*pi*T0))*exp(-mD.*vr.^2./(2*T0)).*(2*pi*vr);
%         fv0_2=(mD/(2*pi*T0))^0.5*exp(-mD.*v0.^2./(2*T0));
    elseif isotope_mass==1
        fv0=(mH/(2*pi*T0))^0.5*exp(-mH.*v0.^2./(2*T0));
        fvr=(mH/(2*pi*T0))*exp(-mH.*vr.^2./(2*T0)).*(2*pi*vr);
    end
    dv0=v0(2)-v0(1);
    
    fv0_2D=repmat(fv0,length(vr),1)';
%     fvr_2D=repmat(fvr,length(v0),1);
    v0_2D=repmat(v0,length(vr),1)';
    vr_2D=repmat(vr,length(v0),1);
    vbeam=v0_2D*0;
    
    
    mb_cum=cumtrapz(v0,fv0);
    mb_cum2=cumtrapz(vr,fvr);
%     plot(vr,mb_cum2)
    
    % plot(v0,fv0)
    
    %%
    sigma_ion_mb=sigma_ionization_val2*0;
    
    for ind_vb=1:length(vb)
        vbeam=sqrt((vb(ind_vb)-v0_2D).^2+vr_2D.^2);
        if isotope_mass==2
            Ebeam=0.5*mD*vbeam.^2/eV;
        else
            Ebeam=0.5*mH*vbeam.^2/eV;
        end
        Ebeam=min(Ebeam,max(Eb_values_eV));
        sigma_ion_vb=interp1(Eb_values_eV,sigma_ionization_val2,Ebeam,'*linear');
        sigma_ion_1D=trapz(v0,fv0_2D.*sigma_ion_vb.*vbeam);
        sigma_ion_mb(ind_vb)=trapz(vr,fvr.*sigma_ion_1D);
%         fv0=(mD/(2*pi*T0))^0.5*exp(-mD.*(v0+vb(ind_vb)).^2./(2*T0));
%         fv0_2=(mD/(2*pi*T0))^0.5*exp(-mD.*(v0-vb(ind_vb)).^2./(2*T0));
%         sigma_cx_mb(ind_vb)=trapz(v0,(fv0-fv0_2).*sigma_cx_vb.*vbeam.^2)./vb(ind_vb);
    end
    sigma_ion_mb(1)
    sigma_ionization_T0(:,ind_T0+1)=sigma_ion_mb;

    
    %%
    % loglog(vbeam,sigma_cx_vb,'linewidth',2);
    % grid on
    % set(gca,'fontsize',16)
    % xlim([1 1e5])
    %
    %%
    loglog(Eb_values_eV,sigma_ion_mb,'linewidth',2);
    grid on
    set(gca,'fontsize',16)
    hold on
    
    
end


T0_values_eV=[0 ; T0_values_eV'];

% legend(num2str(T0_values_eV))
% legend('T_0 = 0 eV','T_0 = 100 eV','T_0 = 1 keV','T_0 = 10 keV')
legend('T_i = 0 eV','T_0 = 100 eV','T_0 = 500 eV','T_0 = 1 keV','T_0 = 5 keV')

xlim([100 1e5])
ylabel('$\left< \sigma v \right> [\mathrm{m}^{3}\mathrm{s}^{-1}]$','interpreter','latex')
xlabel('Eb [eV]')
ylim([2*10^-19 1.2*10^-13]);


%%

Alog=100;
Emin_log=20;
Eb_values_log=Emin_log*log(Eb_values_eV/Alog);

Eb_values_log_lin=linspace(0,150,151);
sigma_ionD_Ti_loglog=zeros(length(Eb_values_log_lin),length(T0_values)+1);

for ind_T0=1:length(T0_values)
    sigma_ionD_Ti_loglog(:,ind_T0)=interp1(Eb_values_eV,sigma_ionization_T0(:,ind_T0),Alog*exp(Eb_values_log_lin/Emin_log));
end

% figure;
% plot(Alog*exp(Eb_values_log_lin/Emin_log),sigma_cxD_T0_log)


%%

% this is using the same scale as the CX values
% T0 (here Ti) between 10 eV and 10 keV

Emin_eV=50;
Elog_scaling=20;
Eb_values_log=Elog_scaling*log(Eb_values_eV/Emin_eV);

T0min_eV=10;
T0log_scaling=3;
T0_values_log=T0log_scaling*log10(T0_values_eV/T0min_eV+T0min_eV);


Eb_values_log_lin=linspace(0,180,181);
T0_values_log_lin=linspace(0,9,10);

sigma_ionD_Ti_log=zeros(length(Eb_values_log_lin),length(T0_values_eV));
sigma_ionD_Ti_loglog=zeros(length(Eb_values_log_lin),length(T0_values_log_lin));

for ind_T0=1:length(T0_values_eV)
    sigma_ionD_Ti_log(:,ind_T0)=interp1(Eb_values_eV,sigma_ionization_T0(:,ind_T0),Emin_eV*exp(Eb_values_log_lin/Elog_scaling));
end

for ind_E=1:length(Eb_values_log_lin)
    sigma_ionD_Ti_loglog(ind_E,:)=interp1(T0_values_eV,sigma_ionD_Ti_log(ind_E,:)',T0min_eV*10.^(T0_values_log_lin/T0log_scaling));
end

if SAVEFILE
save('./sigma_cxD_T0_log.mat','-append','sigma_ionD_Ti_loglog')
end

%%
figure
plot(Emin_eV*exp(Eb_values_log_lin/Elog_scaling),sigma_ionD_Ti_loglog)