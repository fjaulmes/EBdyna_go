
clear all;
close all;
warning off MATLAB:griddata:DuplicateDataPoints;
format long;

initialize_folder_names;
calc_psi_star_evolution_PR;

figure(8);
subplot(2,1,1);
set(gca,'FontSize',22);
grid on;
hold on;
plot(radial_r_value,q_initial,'b--','LineWidth',1.5);
xlim([0 r_value_q1_mean*1.6]);
xlabel('r (m)');
ylabel('q');
%plot(pos_Bstar_final,q_final_interp,'r');
plot(radial_r_value,q_final,'r','LineWidth',1.5);
ylim([0.65 1.6]);

legend('before collapse','after collapse');

figure(8);
subplot(2,1,2);
set(gca,'FontSize',22);
plot(radial_r_value,Psih,'b--','LineWidth',1.5);
grid on;
hold on;
plot(radial_r_value,Psih_final,'r','LineWidth',1.5);
xlim([0 r_value_q1_mean*1.6]);
xlabel('r (m)');
ylabel('\Psi_* (T.m^{-2})');
legend('\psi_{*-}','\psi_{*+}');
ylim([-0.006 0.2])


r_q1=interp1(1:Nradial,radial_r_value_flux,xPsih_max)
initialize_map_r_omega;

xPsih_max=interp1(scale_r,1:NR,r_q1,'cubic');
xPsih_zero=interp1(Psih(round(xPsih_max):NR),round(xPsih_max):NR,0,'cubic');
xMixing_radius=ceil(xPsih_zero);

% Derivative of initial Psih curve

dPsih=zeros(1,size_r);

for (x=2:NR-1)
    Dr_value=Dr_avg(x+1)+Dr_avg(x);
    dPsih(x)=(Psih(x+1)-Psih(x-1))/Dr_value;
end
dPsih(1)=0.5*dPsih(2);

% Derivative of final Psih curve

dPsih_final=zeros(1,size_r);

for x=2:NR-1
    Dr_value=Dr_avg(x+1)+Dr_avg(x);
    dPsih_final(x)=(Psih(x+1)-Psih(x-1))/Dr_value;
end
dPsih_final(1)=0.5*dPsih_final(2);


psi_star_2D=zeros(NR,Nomega);



%alpha=atan(2);
alpha=0.00*pi;
%filename='reconnection_movie\t0p00.bmp';


ksi0=0;

pos_rmix=xPsih_zero;
rmix=1;



% for (omega=1:Nomega)
%     psi_star_2D(:,omega)=interp1(scale_X(1:xPsih_zero),Psih(1:xPsih_zero),scale_r);
% end
for (omega=1:Nomega)
    psi_star_2D(:,omega)=Psih(1:NR);
end

map_region1=zeros(NR,Nomega);

for (r=1:NR)
    for (omega=1:Nomega)
        if (r > xPsih_max)
            map_region1(r,omega)=0;
        else
            map_region1(r,omega)=1;
        end
    end
end

% Create grids and convert polar coordinates to rectangular
[THETA,RR] = meshgrid(scale_omega,scale_r);
[XX,YY] = pol2cart(THETA,RR);

pos_psi_rx=xPsih_max;

% figure(4)
% surf(XX,YY,psi_star_2D,'edgecolor','none');
% view(0,90);
%imagesc(scale_r,scale_omega,psi)

% F=getframe;
% [im,map] = frame2im(F);    %Return associated image data 
% 
% imwrite(im,filename,'bmp');
