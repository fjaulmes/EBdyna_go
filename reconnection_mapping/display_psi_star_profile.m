%close all;

figure(1);
set(gca,'FontSize',16);
hold on;
grid on;
%if (f-2*round(0.5*f))==0
plot([-scale_r(NR:-1:1)  scale_r(2:NR)]/rmix,[psi_star_2D(NR:-1:1,Nomega) ; psi_star_2D(2:NR,1)],'Linewidth',1);
%end
xlabel('r/r_{mix}')
ylabel('\Psi _* (T.m^{-2})')
axis([-1 1 -0.0005 0.005])


der_psi_star2=zeros(1,NR);
for(r=2:NR-1)
    Dr_value=scale_r(r+1)-scale_r(r-1);
    der_psi_star2(r)=(psi_star_2D(r+1,Nomega)-psi_star_2D(r-1,Nomega))/Dr_value;
end

der_f_rank=round(0.5*f);

time_scale_der(der_f_rank)=time_value;
xi_scale_der(der_f_rank)=abs(ksi0);
a1_coef_evol(der_f_rank)=a1_coef;
a2_coef_evol(der_f_rank)=a2_coef;

figure(2);
subplot(2,1,2)
set(gca,'FontSize',16);
grid on;
hold on;
plot(2:NR-1,-der_psi_star2(2:NR-1),'r');
grid on;

der_psi_star=zeros(1,NR);
for(r=2:NR-1)
    Dr_value=scale_r(r+1)-scale_r(r-1);
    der_psi_star(r)=(psi_star_2D(r+1,1)-psi_star_2D(r-1,1))/Dr_value;
end
ylim([-0.05 0.05])

pos_r_nu_cont13=round(x_nu_cont13);
pos_r_nu_cont13_initial=round(x_nu_cont13_initial);
%ksi0<0
r13_value=rx+2*ksi0

if r13_value>=0
    der_cont13_1(der_f_rank)=interp1((1:NR-1),der_psi_star(1:NR-1),0.998*x_nu_cont13);
    der_cont13_3(der_f_rank)=interp1((1:NR-1),der_psi_star(1:NR-1),1.002*x_nu_cont13);
else
    der_cont13_1(der_f_rank)=interp1((1:NR-1),-der_psi_star2(1:NR-1),0.998*x_nu_cont13);
    der_cont13_3(der_f_rank)=interp1((1:NR-1),-der_psi_star2(1:NR-1),1.002*x_nu_cont13);
end
der_cont23_3(der_f_rank)=interp1(scale_r(1:NR-1),der_psi_star(1:NR-1),0.998*rx);
der_cont23_2(der_f_rank)=interp1(scale_r(1:NR-1),dPsih_final(1:NR-1),1.002*rx);

%figure(2);
subplot(2,1,1)
set(gca,'FontSize',16);
ylim([-0.05 0.05])

hold on;
plot(2:NR-1,der_psi_star(2:NR-1));
grid on;
hold off;