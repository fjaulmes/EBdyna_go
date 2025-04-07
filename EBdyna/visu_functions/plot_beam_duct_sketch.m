figure
set(gca,'fontsize',16)
RSHIFT=-0.65

grid on 
grid minor
hold on
axis equal
plot_circle(RSHIFT,0,0.325);
plot_circle(RSHIFT,0,0.79)

% main part of BD
plot([0.65 + RSHIFT; 3],[-0.45 ; -0.45],'b','linewidth',3)
plot([0.45 + RSHIFT; 3],[-0.65 ; -0.65],'b','linewidth',3)

% Rtan
plot([RSHIFT ; 3],[-0.55 ; -0.55],'r--','linewidth',2)

% TC measurements inside BD
plot([0.5 ; 0.5],[-0.65 ; -0.45],'g-.','linewidth',3)


ylim([-1.2 1.2])

xlabel('R_{BD} [m]')
ylabel('[m]')