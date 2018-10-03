close all;





figure(1)
set(gca,'FontSize',16);
title('Counter-Passing / Co-passing');



PRECESS_LAMBDA_BIN_SIZE=0.04;
precession_lambda_bins=(0:PRECESS_LAMBDA_BIN_SIZE:31*PRECESS_LAMBDA_BIN_SIZE);
lambda_values=(0.5*PRECESS_LAMBDA_BIN_SIZE:PRECESS_LAMBDA_BIN_SIZE:30.5*PRECESS_LAMBDA_BIN_SIZE);

precession_Ekin_bins=[0 400 1000 2000 3600]*1e3;
Ekin_values=[200 700 1500 2800]*1e3;

lambda_precess_value=zeros(length(precession_lambda_bins)-1,length(precession_Ekin_bins)-1);

for Ebin=1:(length(precession_Ekin_bins)-1)
    PRECESS_EKIN_POP=((alphas_Ekin>precession_Ekin_bins(Ebin)).*(alphas_Ekin<=precession_Ekin_bins(Ebin+1)));
    for bin=1:(length(precession_lambda_bins)-1)
        PRECESS_LAMBDA_POP_CO=find((alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*PRECESS_EKIN_POP.*CO_PASSING_POP);
        PRECESS_LAMBDA_POP_COUNTER=find((alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*PRECESS_EKIN_POP.*COUNTER_PASSING_POP);
        lambda_precess_value(bin,Ebin)=(mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_COUNTER))/mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP_CO)))));
        lambda_precess_value(bin,Ebin)=lambda_precess_value(bin,Ebin)*(mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_CO))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP_COUNTER)))));
    end
end


hold on
grid on
plot(lambda_values,lambda_precess_value(:,1),'g','LineWidth',1);
plot(lambda_values,lambda_precess_value(:,2),'b','LineWidth',2);
plot(lambda_values,lambda_precess_value(:,3),'k','LineWidth',3);
plot(lambda_values,lambda_precess_value(:,4),'r','LineWidth',4);
plot(lambda_values,lambda_values*0+1,'y--','LineWidth',2);
plot([1 1],[-1 max(max(lambda_precess_value(1:end,:)))],'y--','LineWidth',2);
xlim([0.1 0.9])

Ekin_values=Ekin_values*1e-3;
legend(strcat('Ekin=',num2str(Ekin_values(1)),'keV'),strcat('Ekin=',num2str(Ekin_values(2)),'keV'),strcat('Ekin=',num2str(Ekin_values(3)),'keV'),strcat('Ekin=',num2str(Ekin_values(4)),'keV'));

xlabel('\lambda_{0}')
% ylabel('(\omega_{precess} (COUNTER/CO)')
ylabel('(\omega_{ratios} (COUNTER/CO)')









figure(3)
set(gca,'FontSize',16);
title('Trapped');


precession_Ekin_bins=[0 400 1000 2000 3600]*1e3;
Ekin_values=[200 700 1500 2800]*1e3;

lambda_precess_value=zeros(length(precession_lambda_bins)-1,length(precession_Ekin_bins)-1);

for Ebin=1:(length(precession_Ekin_bins)-1)
    PRECESS_EKIN_POP=((alphas_Ekin>precession_Ekin_bins(Ebin)).*(alphas_Ekin<=precession_Ekin_bins(Ebin+1)));
    for bin=1:(length(precession_lambda_bins)-1)
        PRECESS_LAMBDA_POP=find((alphas_lambda>precession_lambda_bins(bin)).*(alphas_lambda<=precession_lambda_bins(bin+1)).*PRECESS_EKIN_POP.*ALL_TRAPPED_POP);
        lambda_precess_value(bin,Ebin)=(mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP))./mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP)))));
    end
end



set(gca,'FontSize',16);

hold on
grid on
hold on
grid on
plot(lambda_values,lambda_precess_value(:,1),'g','LineWidth',1);
plot(lambda_values,lambda_precess_value(:,2),'b','LineWidth',2);
plot(lambda_values,lambda_precess_value(:,3),'k','LineWidth',3);
plot(lambda_values,lambda_precess_value(:,4),'r','LineWidth',4);
plot(lambda_values,lambda_values*0+1,'y--','LineWidth',2);
plot([1 1],[-1 max(max(lambda_precess_value(1:end,:)))],'y--','LineWidth',2);
xlim([lambda_values(1) lambda_values(end)])

Ekin_values=Ekin_values*1e-3;
legend(strcat('Ekin=',num2str(Ekin_values(1)),'keV'),strcat('Ekin=',num2str(Ekin_values(2)),'keV'),strcat('Ekin=',num2str(Ekin_values(3)),'keV'),strcat('Ekin=',num2str(Ekin_values(4)),'keV'));

xlabel('\lambda_{0}')
ylabel('(\omega_{precess}/\omega_{\psi})')


