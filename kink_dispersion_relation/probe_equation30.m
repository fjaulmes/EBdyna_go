L_values_neg=-(0:0.001:4);
L_values=(0:0.001:1);
plot(L_values_neg,real((sqrt(1+i*L_values_neg.^2)+sqrt(i)*L_values_neg).^4));hold on
plot(L_values,real((sqrt(1+i*L_values.^2)+sqrt(i)*L_values).^4));hold on