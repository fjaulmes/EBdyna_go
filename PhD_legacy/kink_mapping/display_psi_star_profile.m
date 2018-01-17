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
ylim([-0.03 0.03])

%figure(2);
subplot(2,1,1)
set(gca,'FontSize',16);
ylim([-0.03 0.03])

hold on;
plot(2:NR-1,der_psi_star(2:NR-1));
grid on;
hold off;