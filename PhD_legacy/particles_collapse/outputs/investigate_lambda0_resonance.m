close all


PRECESS_LAMBDA_BIN_SIZE=0.02;
precession_lambda_bins=(40*PRECESS_LAMBDA_BIN_SIZE:PRECESS_LAMBDA_BIN_SIZE:60*PRECESS_LAMBDA_BIN_SIZE);
lambda_values=(40.5*PRECESS_LAMBDA_BIN_SIZE:PRECESS_LAMBDA_BIN_SIZE:59.5*PRECESS_LAMBDA_BIN_SIZE);

lambda_precess_value=zeros(length(precession_lambda_bins)-1,1);
lambda_longitudinal_value=zeros(length(precession_lambda_bins)-1,1);
Delta_pphi_value=zeros(length(precession_lambda_bins)-1,1);

for bin=1:(length(precession_lambda_bins)-1)
    PRECESS_LAMBDA_POP=find((psi_pos0<=psi_rank_q1).*(alphas_lambda0>precession_lambda_bins(bin)).*(alphas_lambda0<=precession_lambda_bins(bin+1)).*ALL_TRAPPED_POP);
    % lambda_precess_value(bin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP)))/mean(abs(omega_psi_avg(PRECESS_LAMBDA_POP)));
    lambda_precess_value(bin)=mean(abs(omega_precess_avg(PRECESS_LAMBDA_POP)));
    lambda_longitudinal_value(bin)=mean((omega_psi_avg(PRECESS_LAMBDA_POP)));
    Delta_pphi_value(bin)=mean(abs(Delta_pphi(PRECESS_LAMBDA_POP)));
end



figure(1)
subplot(2,1,1)
set(gca,'FontSize',22);
hold on
grid on
plot(lambda_values,Delta_pphi_value,'b','LineWidth',2)
xlim([0.82 1.05])
xlabel('\lambda_0')
ylabel('\Delta (p_\phi)')


subplot(2,1,2)
set(gca,'FontSize',22);
hold on
grid on
plot(lambda_values,lambda_values.*0+omega_crash,'y--','LineWidth',2)
% plot(lambda_values,lambda_values.*0-omega_crash,'y--','LineWidth',2)
plot(lambda_values,lambda_precess_value,'r+','LineWidth',3)
plot(lambda_values,lambda_precess_value,'r','LineWidth',3)
% plot(lambda_values,abs(lambda_precess_value)./abs(lambda_longitudinal_value),'k--','LineWidth',3)
xlim([0.82 1.05])
xlabel('\lambda_0')
ylabel('\omega_{v_D}')

