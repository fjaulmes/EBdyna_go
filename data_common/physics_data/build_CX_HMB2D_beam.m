%CX against D at maxwellian temperature

SAVEFILE=0


isotope_mass=2;  % amu

vb_H=linspace(0.01,0.01+8e6,1e3);

if isotope_mass==2
    Eb=0.5*(mD/eV).*vb_H.^2;
elseif  isotope_mass==1
    Eb=0.5*(mH/eV).*vb_H.^2;
end

Eb_keV=Eb*1e-3;

A1_CX=3.2345
A2_CX=235.88
A3_CX=0.038371
A4_CX=3.8068e-6;
A5_CX=1.1832e-10;
A6_CX=2.3713;

sigma_cx=10^-20 * (A1_CX .*log (A2_CX./Eb_keV + A6_CX))./(1+A3_CX.*Eb_keV+A4_CX.*Eb_keV.^3.5 + A5_CX.*Eb_keV.^ (5.4));


Eb_values_eV=Eb*isotope_mass;


%%
figure

loglog(Eb*isotope_mass,sigma_cx,'linewidth',2);
grid on
hold on
set(gca,'fontsize',16)

title('CX (D) with non excitated states')
ylabel('\sigma [m^{2}]')
xlabel('Eb [eV]')
xlim([100 1e5])
ylim([10^-20 1.2*10^-18])

%%
figure
vb=vb_H
if isotope_mass==2
    vb=sqrt(2*Eb_values_eV*eV/mD);
elseif isotope_mass==1
    vb=sqrt(2*Eb_values_eV*eV/mH);
end
sigma_cx_vb_T0=sigma_cx.*vb;

loglog(Eb_values_eV,sigma_cx_vb_T0,'linewidth',2);
grid on
hold on
set(gca,'fontsize',16)

title('CX D+ + D for various Maxwellian D temp.')

% maxwellian distribution

%%
T0_values_eV=linspace(1e3,6e3,6); % (background neutrals)
% T0_values_eV=[100 1e3 10e3]*2;
T0_values_eV=[100 1e3 10*1e3];
T0_values_eV=[100 500 1*1e3 3e3 5e3 10e3];
T0_values=(eV)*T0_values_eV; % in eV (background neutrals)

sigma_cxD_T0=zeros(length(sigma_cx),length(T0_values)+1);
sigma_cxD_T0(:,1)=sigma_cx_vb_T0;

for ind_T0=1:length(T0_values)
    T0=T0_values(ind_T0);
    T0_ev=T0/eV
    
     v0=linspace(-4.0e6,4.0e6,5e2);
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
    sigma_cx_mb=sigma_cx*0;
    
    for ind_vb=1:length(vb)
        vbeam=sqrt((vb(ind_vb)-v0_2D).^2+vr_2D.^2);
        if isotope_mass==2
            Ebeam=0.5*mD*vbeam.^2/eV;
        else
            Ebeam=0.5*mH*vbeam.^2/eV;
        end
        Ebeam=min(Ebeam,max(Eb_values_eV));
        sigma_cx_vb=interp1(Eb_values_eV,sigma_cx,Ebeam);
        sigma_cx_mb_1D=trapz(v0,fv0_2D.*sigma_cx_vb.*vbeam);
        sigma_cx_mb(ind_vb)=trapz(vr,fvr.*sigma_cx_mb_1D);
%         fv0=(mD/(2*pi*T0))^0.5*exp(-mD.*(v0+vb(ind_vb)).^2./(2*T0));
%         fv0_2=(mD/(2*pi*T0))^0.5*exp(-mD.*(v0-vb(ind_vb)).^2./(2*T0));
%         sigma_cx_mb(ind_vb)=trapz(v0,(fv0-fv0_2).*sigma_cx_vb.*vbeam.^2)./vb(ind_vb);
    end
    sigma_cxD_T0(:,ind_T0+1)=sigma_cx_mb;

    
    %%
    % loglog(vbeam,sigma_cx_vb,'linewidth',2);
    % grid on
    % set(gca,'fontsize',16)
    % xlim([1 1e5])
    %
    %%
    loglog(Eb_values_eV,sigma_cx_mb,'linewidth',2);
    grid on
    set(gca,'fontsize',16)
    
    
end
%%
T0_values_eV=[0 ; T0_values_eV'];

% legend(num2str(T0_values_eV))
% legend('T_0 = 0 eV','T_0 = 100 eV','T_0 = 1 keV','T_0 = 10 keV')
legend('T_0 = 0 eV','T_0 = 100 eV','T_0 = 500 eV','T_0 = 1 keV','T_0 = 3 keV','T_0 = 5 keV','T_0 = 10 keV')

xlim([100 1e5])
ylabel('$\left< \sigma v \right> [\mathrm{m}^{3}\mathrm{s}^{-1}]$','interpreter','latex')
xlabel('Eb [eV]')
ylim([2*10^-14 1.2*10^-13])


%%

Emin_eV=50;
Elog_scaling=20;
Eb_values_log=Elog_scaling*log(Eb_values_eV/Emin_eV);

T0min_eV=10;
T0log_scaling=3;
T0_values_log=T0log_scaling*log10(T0_values_eV/T0min_eV+T0min_eV);


Eb_values_log_lin=linspace(0,180,181);
T0_values_log_lin=linspace(0,9,10);

sigma_cxD_T0_log=zeros(length(Eb_values_log_lin),length(T0_values_eV));
sigma_cxD_T0_loglog=zeros(length(Eb_values_log_lin),length(T0_values_log_lin));

for ind_T0=1:length(T0_values_eV)
    sigma_cxD_T0_log(:,ind_T0)=interp1(Eb_values_eV,sigma_cxD_T0(:,ind_T0),Emin_eV*exp(Eb_values_log_lin/Elog_scaling));
end

for ind_E=1:length(Eb_values_log_lin)
    sigma_cxD_T0_loglog(ind_E,:)=interp1(T0_values_eV,sigma_cxD_T0_log(ind_E,:)',T0min_eV*10.^(T0_values_log_lin/T0log_scaling));
end

% figure;
% plot(Emin_eV*exp(Eb_values_log_lin/Elog_scaling),sigma_cxD_T0_log)

% save('./sigma_cxD_T0_log.mat','sigma_cxD_T0_log','Eb_values_log','T0_values_eV','Elog_scaling','Emin_eV')
% save('./sigma_cxD_T0_log.mat','sigma_cxD_T0_log','T0_values_eV','Eb_values_log_lin','Elog_scaling','Emin_eV')
if SAVEFILE
save('./sigma_cxD_T0_log.mat','sigma_cxD_T0_loglog','T0_values_log_lin','T0log_scaling','T0min_eV','Eb_values_log_lin','Elog_scaling','Emin_eV')
end

