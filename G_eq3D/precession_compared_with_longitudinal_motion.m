close all;


PRECESS_RADIAL_BIN_SIZE=0.1;
precession_radial_bins=(0:PRECESS_RADIAL_BIN_SIZE:9*PRECESS_RADIAL_BIN_SIZE);
radial_values=(0.5*PRECESS_RADIAL_BIN_SIZE:PRECESS_RADIAL_BIN_SIZE:8.5*PRECESS_RADIAL_BIN_SIZE);

precession_Ekin_bins=[0 400 1000 2000 3600]*1e3;
Ekin_values=[200 700 1500 2800]*1e3;

radial_precess_value=zeros(length(precession_radial_bins)-1,length(precession_Ekin_bins)-1);

for Ebin=1:(length(precession_Ekin_bins)-1)
    PRECESS_EKIN_POP=((alphas_Ekin>precession_Ekin_bins(Ebin)).*(alphas_Ekin<=precession_Ekin_bins(Ebin+1)));
    for bin=1:(length(precession_radial_bins)-1)
        PRECESS_LAMBDA_POP=find((r_avg'>precession_radial_bins(bin)).*(r_avg'<=precession_radial_bins(bin+1)).*PRECESS_EKIN_POP.*CO_PASSING_POP);
        radial_precess_value(bin,Ebin)=mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP)));
    end
end



figure(1)
set(gca,'FontSize',16);
title('Co-Passing');

hold on
grid on
hold on
grid on
plot(radial_values,radial_precess_value(:,1),'g--','LineWidth',1);
plot(radial_values,radial_precess_value(:,2),'b--','LineWidth',2);
plot(radial_values,radial_precess_value(:,3),'k--','LineWidth',3);
plot(radial_values,radial_precess_value(:,4),'r--','LineWidth',4);
plot([r_value_q1_mean r_value_q1_mean],[0 max(max(radial_precess_value(1:end-1,:)))],'y--','LineWidth',2);
xlim([radial_values(1) radial_values(end-1)])
ylim([0 max(max(radial_precess_value(1:end-2,:)))])

Ekin_values=Ekin_values*1e-3;
legend(strcat('Ekin=',num2str(Ekin_values(1)),'keV'),strcat('Ekin=',num2str(Ekin_values(2)),'keV'),strcat('Ekin=',num2str(Ekin_values(3)),'keV'),strcat('Ekin=',num2str(Ekin_values(4)),'keV'));

xlabel('r_{avg}')
ylabel('\omega_{\psi} and \omega_{pre}');





PRECESS_RADIAL_BIN_SIZE=0.1;
precession_radial_bins=(0:PRECESS_RADIAL_BIN_SIZE:9*PRECESS_RADIAL_BIN_SIZE);
radial_values=(0.5*PRECESS_RADIAL_BIN_SIZE:PRECESS_RADIAL_BIN_SIZE:8.5*PRECESS_RADIAL_BIN_SIZE);

precession_Ekin_bins=[0 400 1000 2000 3600]*1e3;
Ekin_values=[200 700 1500 2800]*1e3;

radial_precess_value=zeros(length(precession_radial_bins)-1,length(precession_Ekin_bins)-1);

for Ebin=1:(length(precession_Ekin_bins)-1)
    PRECESS_EKIN_POP=((alphas_Ekin>precession_Ekin_bins(Ebin)).*(alphas_Ekin<=precession_Ekin_bins(Ebin+1)));
    for bin=1:(length(precession_radial_bins)-1)
        PRECESS_LAMBDA_POP=find((r_avg'>precession_radial_bins(bin)).*(r_avg'<=precession_radial_bins(bin+1)).*PRECESS_EKIN_POP.*CO_PASSING_POP);
        radial_precess_value(bin,Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP)));
    end
end


hold on
grid on
plot(radial_values,radial_precess_value(:,1),'g','LineWidth',1);
plot(radial_values,radial_precess_value(:,2),'b','LineWidth',2);
plot(radial_values,radial_precess_value(:,3),'k','LineWidth',3);
plot(radial_values,radial_precess_value(:,4),'r','LineWidth',4);
plot([r_value_q1_mean r_value_q1_mean],[0 max(max(radial_precess_value(1:end-1,:)))],'y--','LineWidth',2);
xlim([radial_values(1) radial_values(end-1)])
ylim([0 max(max(radial_precess_value(1:end-1,:)))])

Ekin_values=Ekin_values*1e-3;
legend(strcat('Ekin=',num2str(Ekin_values(1)),'keV'),strcat('Ekin=',num2str(Ekin_values(2)),'keV'),strcat('Ekin=',num2str(Ekin_values(3)),'keV'),strcat('Ekin=',num2str(Ekin_values(4)),'keV'));

xlabel('r_{avg}')







PRECESS_RADIAL_BIN_SIZE=0.1;
precession_radial_bins=(0:PRECESS_RADIAL_BIN_SIZE:9*PRECESS_RADIAL_BIN_SIZE);
radial_values=(0.5*PRECESS_RADIAL_BIN_SIZE:PRECESS_RADIAL_BIN_SIZE:8.5*PRECESS_RADIAL_BIN_SIZE);

precession_Ekin_bins=[0 400 1000 2000 3600]*1e3;
Ekin_values=[200 700 1500 2800]*1e3;

radial_precess_value=zeros(length(precession_radial_bins)-1,length(precession_Ekin_bins)-1);

for Ebin=1:(length(precession_Ekin_bins)-1)
    PRECESS_EKIN_POP=((alphas_Ekin>precession_Ekin_bins(Ebin)).*(alphas_Ekin<=precession_Ekin_bins(Ebin+1)));
    for bin=1:(length(precession_radial_bins)-1)
        PRECESS_LAMBDA_POP=find((r_avg'>precession_radial_bins(bin)).*(r_avg'<=precession_radial_bins(bin+1)).*PRECESS_EKIN_POP.*COUNTER_PASSING_POP);
        radial_precess_value(bin,Ebin)=mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP)));
    end
end



figure(2)
set(gca,'FontSize',16);
title('Counter-Passing');

hold on
grid on
hold on
grid on
plot(radial_values,radial_precess_value(:,1),'g--','LineWidth',1);
plot(radial_values,radial_precess_value(:,2),'b--','LineWidth',2);
plot(radial_values,radial_precess_value(:,3),'k--','LineWidth',3);
plot(radial_values,radial_precess_value(:,4),'r--','LineWidth',4);
plot([r_value_q1_mean r_value_q1_mean],[0 max(max(radial_precess_value(1:end-1,:)))],'y--','LineWidth',2);
xlim([radial_values(1) radial_values(end-1)])
ylim([0 max(max(radial_precess_value(1:end-2,:)))])

Ekin_values=Ekin_values*1e-3;
legend(strcat('Ekin=',num2str(Ekin_values(1)),'keV'),strcat('Ekin=',num2str(Ekin_values(2)),'keV'),strcat('Ekin=',num2str(Ekin_values(3)),'keV'),strcat('Ekin=',num2str(Ekin_values(4)),'keV'));

xlabel('r_{avg}')
ylabel('\omega_{\psi} and \omega_{pre}');





PRECESS_RADIAL_BIN_SIZE=0.1;
precession_radial_bins=(0:PRECESS_RADIAL_BIN_SIZE:9*PRECESS_RADIAL_BIN_SIZE);
radial_values=(0.5*PRECESS_RADIAL_BIN_SIZE:PRECESS_RADIAL_BIN_SIZE:8.5*PRECESS_RADIAL_BIN_SIZE);

precession_Ekin_bins=[0 400 1000 2000 3600]*1e3;
Ekin_values=[200 700 1500 2800]*1e3;

radial_precess_value=zeros(length(precession_radial_bins)-1,length(precession_Ekin_bins)-1);

for Ebin=1:(length(precession_Ekin_bins)-1)
    PRECESS_EKIN_POP=((alphas_Ekin>precession_Ekin_bins(Ebin)).*(alphas_Ekin<=precession_Ekin_bins(Ebin+1)));
    for bin=1:(length(precession_radial_bins)-1)
        PRECESS_LAMBDA_POP=find((r_avg'>precession_radial_bins(bin)).*(r_avg'<=precession_radial_bins(bin+1)).*PRECESS_EKIN_POP.*COUNTER_PASSING_POP);
        radial_precess_value(bin,Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP)));
    end
end


hold on
grid on
plot(radial_values,radial_precess_value(:,1),'g','LineWidth',1);
plot(radial_values,radial_precess_value(:,2),'b','LineWidth',2);
plot(radial_values,radial_precess_value(:,3),'k','LineWidth',3);
plot(radial_values,radial_precess_value(:,4),'r','LineWidth',4);
plot([r_value_q1_mean r_value_q1_mean],[0 max(max(radial_precess_value(1:end-1,:)))],'y--','LineWidth',2);
xlim([radial_values(1) radial_values(end-1)])
ylim([0 max(max(radial_precess_value(1:end-1,:)))])

Ekin_values=Ekin_values*1e-3;
legend(strcat('Ekin=',num2str(Ekin_values(1)),'keV'),strcat('Ekin=',num2str(Ekin_values(2)),'keV'),strcat('Ekin=',num2str(Ekin_values(3)),'keV'),strcat('Ekin=',num2str(Ekin_values(4)),'keV'));

xlabel('r_{avg}')






PRECESS_RADIAL_BIN_SIZE=0.1;
precession_radial_bins=(0:PRECESS_RADIAL_BIN_SIZE:9*PRECESS_RADIAL_BIN_SIZE);
radial_values=(0.5*PRECESS_RADIAL_BIN_SIZE:PRECESS_RADIAL_BIN_SIZE:8.5*PRECESS_RADIAL_BIN_SIZE);

precession_Ekin_bins=[0 400 1000 2000 3600]*1e3;
Ekin_values=[200 700 1500 2800]*1e3;

radial_precess_value=zeros(length(precession_radial_bins)-1,length(precession_Ekin_bins)-1);

for Ebin=1:(length(precession_Ekin_bins)-1)
    PRECESS_EKIN_POP=((alphas_Ekin>precession_Ekin_bins(Ebin)).*(alphas_Ekin<=precession_Ekin_bins(Ebin+1)));
    for bin=1:(length(precession_radial_bins)-1)
        PRECESS_LAMBDA_POP=find((r_avg'>precession_radial_bins(bin)).*(r_avg'<=precession_radial_bins(bin+1)).*PRECESS_EKIN_POP.*ALL_TRAPPED_POP);
        radial_precess_value(bin,Ebin)=mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP)));
    end
end



figure(3)
set(gca,'FontSize',16);
title('Trapped');

hold on
grid on
hold on
grid on
plot(radial_values,radial_precess_value(:,1),'g--','LineWidth',1);
plot(radial_values,radial_precess_value(:,2),'b--','LineWidth',2);
plot(radial_values,radial_precess_value(:,3),'k--','LineWidth',3);
plot(radial_values,radial_precess_value(:,4),'r--','LineWidth',4);
plot([r_value_q1_mean r_value_q1_mean],[0 max(max(radial_precess_value(1:end-1,:)))],'y--','LineWidth',2);
xlim([radial_values(1) radial_values(end-1)])
ylim([0 max(max(radial_precess_value(1:end-2,:)))])

Ekin_values=Ekin_values*1e-3;
legend(strcat('Ekin=',num2str(Ekin_values(1)),'keV'),strcat('Ekin=',num2str(Ekin_values(2)),'keV'),strcat('Ekin=',num2str(Ekin_values(3)),'keV'),strcat('Ekin=',num2str(Ekin_values(4)),'keV'));

xlabel('r_{avg}')
ylabel('\omega_{\psi} and \omega_{pre}');





PRECESS_RADIAL_BIN_SIZE=0.1;
precession_radial_bins=(0:PRECESS_RADIAL_BIN_SIZE:9*PRECESS_RADIAL_BIN_SIZE);
radial_values=(0.5*PRECESS_RADIAL_BIN_SIZE:PRECESS_RADIAL_BIN_SIZE:8.5*PRECESS_RADIAL_BIN_SIZE);

precession_Ekin_bins=[0 400 1000 2000 3600]*1e3;
Ekin_values=[200 700 1500 2800]*1e3;

radial_precess_value=zeros(length(precession_radial_bins)-1,length(precession_Ekin_bins)-1);

for Ebin=1:(length(precession_Ekin_bins)-1)
    PRECESS_EKIN_POP=((alphas_Ekin>precession_Ekin_bins(Ebin)).*(alphas_Ekin<=precession_Ekin_bins(Ebin+1)));
    for bin=1:(length(precession_radial_bins)-1)
        PRECESS_LAMBDA_POP=find((r_avg'>precession_radial_bins(bin)).*(r_avg'<=precession_radial_bins(bin+1)).*PRECESS_EKIN_POP.*ALL_TRAPPED_POP);
        radial_precess_value(bin,Ebin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP)));
    end
end



set(gca,'FontSize',16);

hold on
grid on
hold on
grid on
plot(radial_values,radial_precess_value(:,1),'g','LineWidth',1);
plot(radial_values,radial_precess_value(:,2),'b','LineWidth',2);
plot(radial_values,radial_precess_value(:,3),'k','LineWidth',3);
plot(radial_values,radial_precess_value(:,4),'r','LineWidth',4);
plot([r_value_q1_mean r_value_q1_mean],[0 max(max(radial_precess_value(1:end-1,:)))],'y--','LineWidth',2);
plot(radial_values,radial_values.*0+omega_crash,'y--','LineWidth',2);
xlim([radial_values(1) radial_values(end-1)])
ylim([0 max(max(radial_precess_value(1:end-2,:)))])

Ekin_values=Ekin_values*1e-3;
legend(strcat('Ekin=',num2str(Ekin_values(1)),'keV'),strcat('Ekin=',num2str(Ekin_values(2)),'keV'),strcat('Ekin=',num2str(Ekin_values(3)),'keV'),strcat('Ekin=',num2str(Ekin_values(4)),'keV'));

xlabel('r_{avg}')



