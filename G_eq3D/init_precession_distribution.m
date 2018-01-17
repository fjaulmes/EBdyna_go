close all;

posX_output=Xpos_gc_output;
posZ_output=Zpos_gc_output;
PART_POP=(~alphas_ejected);

X_output_pop=posX_output(:,PART_POP);
Z_output_pop=posZ_output(:,PART_POP);
clear posX_output posZ_output 

vpll_output_pop=vpll_output(:,PART_POP);

sign_v=zeros(Nalphas_simulated,1);
max_vpll_avg=zeros(Nalphas_simulated,1);
for number=1:Nalphas_simulated
    [max_vpll max_index]=max(abs(vpll_output_pop(:,number)));
    max_vpll_avg(number)=vpll_output_pop(max_index,number);
%     sign_v(number)=sign(vpll_output_corr(max_index,number));
end
sign_v=sign(max_vpll_avg);


calculate_vD_XZ_map;
calculate_gradZ_theta_map;

bX_output_pop=interp2(scale_X,scale_Z,bX_XZ_map',X_output_pop,Z_output_pop,'*linear');
bZ_output_pop=interp2(scale_X,scale_Z,bZ_XZ_map',X_output_pop,Z_output_pop,'*linear');
bphi_output_pop=sqrt(bX_output_pop.^2+bZ_output_pop.^2);

theta_output_pop=interp2(scale_X,scale_Z,theta_XZsmall_map',X_output_pop,Z_output_pop,'*nearest');
GZ_theta_output_pop=interp2(scale_X,scale_Z,gradZ_theta_map',X_output_pop,Z_output_pop,'*linear');

q_output_pop=interp2(scale_X,scale_Z,q_initial_XZsmall_map',X_output_pop,Z_output_pop,'*linear');
r_output_pop=sqrt((X_output_pop).^2+Z_output_pop.^2);




% r_output=interp1(1:257,radial_r_value_flux,r_output,'*linear');
vDZ_output_corr=interp2(scale_X,scale_Z,vD_Z_XZ_map',X_output_pop,Z_output_pop,'*linear');
vDphi_output_corr=interp2(scale_X,scale_Z,vD_phi_XZ_map',X_output_pop,Z_output_pop,'*linear');




% load('initialG_800keV_flat_pre_collapse.mat');
disp('Total number of particles simulates');
disp(Nalphas_simulated);
mid_X_axis=mid_X+round(X_axis/DX)
Bavg=Btot_XZ_map(mid_X_axis, mid_Z)
% alphas_Bfield=interp2(scale_X,scale_Z,Btot_XZ_map',alphas_pos_x,alphas_pos_z,'*linear');
alphas_v=sqrt(2*(eV/mHe)*alphas_Ekin);
% alphas_vD_Z=interp2(scale_X,scale_Z,vD_Z_XZ_map',alphas_pos_x,alphas_pos_z,'*linear');
r_avg=mean(r_output_pop(:,:),1)';
q_avg=mean(q_output_pop(:,:),1)';
alphas_kappa=sqrt((alphas_Ekin*R0+Bavg*alphas_mm.*(r_avg-R0))./(2*alphas_mm.*r_avg*Bavg));
% alphas_kappa=sqrt(((alphas_Ekin+Bavg*alphas_mm)*(R0-r_avg))./(2*alphas_mm.*r_avg*Bavg));
delta_orbit=0.5*(mHe/eV)*(alphas_v.*q_avg)./(ZHe*Bavg*sqrt(r_avg/R0));
delta_potato=R0*(2*mHe*q_avg.*alphas_v./(R0*ZHe*eV*Bavg)).^(2/3);

avg_dist_X_inf=round((r_avg-delta_orbit)/DX);
avg_dist_X_sup=round((r_avg+delta_orbit)/DX);
Bmax_co_passing=Btot_XZ_map(max(mid_X-avg_dist_X_inf,1), mid_Z);
Bmax_counter_passing=Btot_XZ_map(max(mid_X-avg_dist_X_sup,1), mid_Z);

lambda_co_passing=Bmax_co_passing.*alphas_mm./(alphas_Ekin);
lambda_counter_passing=Bmax_counter_passing.*alphas_mm./(alphas_Ekin);
% lambda_co_passing=Bavg*R0*alphas_mm./(alphas_Ekin.*(R0-r_avg+delta_orbit));
% lambda_counter_passing=Bavg*R0*alphas_mm./(alphas_Ekin.*(R0-r_avg-delta_orbit));
PLUS=find(sign_v==1);
MINUS=find(sign_v==-1);
HFS_crossing=zeros(Nalphas_simulated,1);
minX_value_pop=zeros(Nalphas_simulated,1);
% theta_output_pop=theta_output(:,PART_POP);
vpll_output_pop=vpll_output(:,PART_POP);

clear theta_output

for number=1:Nalphas_simulated
    [HFS_theta_value HFS_theta_pos]=min(abs(theta_output_pop(:,number)-pi));
    [minX_value minX_pos]=min(X_output_pop(:,number));
    minX_value_pop(number)=minX_value;
    % to be adjusted according to guiding center position

    if (HFS_theta_value<0.1)&&(X_output_pop(HFS_theta_pos,number)<=X_axis)
        HFS_crossing(number)=1;
%     if mod(number,1000)==0
%         theta_output_pop(HFS_theta_pos,number)
%         X_output_pop(HFS_theta_pos,number)
%     end
    end
end

alphas_lambda=Bavg*alphas_mm./alphas_Ekin;



ENERGY_POPULATION=find(alphas_Ekin>1*1e3);




disp('EVALUATED_POPULATION')
disp(size(ENERGY_POPULATION,1));


alphas_pos_psi=interp2(scale_X,scale_Z,psi_norm_XZsmall_map',alphas_pos_x,alphas_pos_z,'*linear');
alphas_value_psi=interp2(scale_X,scale_Z,psi_XZsmall_map',alphas_pos_x,alphas_pos_z,'*linear');
theta_value=interp2(scale_X,scale_Z,theta_XZsmall_map',alphas_pos_x,alphas_pos_z,'*nearest');


time_scale_corr=time_scale';
end_ts=size(vpll_output_pop,1);

% phipos_output_wrap=wrap2pi(phipos_output');

% r_output=sqrt((posX_output-X_axis).^2+posZ_output.^2);

% q_output_pop=q_output(:,PART_POP);
% r_output_pop=r_output(:,PART_POP);
% vDZ_output_corr=vDZ_output(:,PART_POP);
% vDphi_output_corr=vDphi_output(:,PART_POP);

% mu_output_pop=atan(r_output_pop./(q_output_pop.*(R0+X_output_pop)));

% omega_output=theta_output-phipos_output;
% omega_output_corr=omega_output(:,PART_POP);
% phipos_output_corr=phipos_output(:,PART_POP);

for ts=1:end_ts
    vDZ_output_corr(ts,:)=(vDZ_output_corr(ts,:)').*(2*alphas_Ekin-Bavg*alphas_mm)/(ZHe);
    vDphi_output_corr(ts,:)=(vDphi_output_corr(ts,:)').*(2*alphas_Ekin-Bavg*alphas_mm)/(ZHe);
end

clear   vDZ_output vDphi_output omega_output

B_output_pop=interp2(scale_X,scale_Z,theta_XZsmall_map',X_output_pop,Z_output_pop,'*nearest');


figure(1)
set(gca,'FontSize',22);

hold on
grid on;

NXpos=histc(alphas_pos_x(ENERGY_POPULATION),(-18:18)*0.06);xlim([-0.9 0.9]);

plot((-18:18)*0.06,NXpos,'r');

title('axial position (X)');


figure(2)
set(gca,'FontSize',22);

hold on;
grid on;

vpll_bin_size=0.05*1e7;
vpll_range=(-1.6:0.05:1.6)*1e7;
vpll_range=vpll_range';

Nvpll=histc(alphas_vpll(ENERGY_POPULATION),vpll_range-0.5*vpll_bin_size);

plot(vpll_range,Nvpll,'r');
title('v_{||}');


figure(3)
set(gca,'FontSize',22);
hold on;
grid on;

radial_bin_size=9;
psi_bin_pos=(1:radial_bin_size:254);
Npsi=histc(alphas_psi(ENERGY_POPULATION),psi_bin_pos-0.5*radial_bin_size);

plot(psi_bin_pos,Npsi,'r');
title('flux surface position (\psi)');



figure(4)
set(gca,'FontSize',22);
hold on;
grid on;

lambda_bin_size=0.05;
lambda_bins=(0:lambda_bin_size:1.2);

Nlambda=histc(alphas_lambda(ENERGY_POPULATION),lambda_bins-0.5*lambda_bin_size);

plot(lambda_bins,Nlambda,'r');
title('pitch angle evolution');


% alphas_vpll0_map=griddata(psi_pos0,X0,abs(vpll0),PSI_xx,XX_psi);


figure(5)
set(gca,'FontSize',22);
hold on;
grid on;

p_phi_scale=(-2:0.25:4);
p_phi_bin_size=0.25;

p_phi=mHe*(alphas_pos_x(ENERGY_POPULATION)+R0).*(alphas_vpll(ENERGY_POPULATION))-2*eV*(alphas_value_psi(ENERGY_POPULATION));

Npphi=histc(p_phi/eV,p_phi_scale-0.5*p_phi_bin_size);


plot(p_phi_scale,Npphi,'r','LineWidth',2);
title('canonical toroidal angular momentum');




figure(6)
set(gca,'FontSize',22);
hold on;
grid on;

kappa_bin_size=0.1;
kappa_bin_pos=(kappa_bin_size:kappa_bin_size:101*kappa_bin_size)-0.5*kappa_bin_size;

Nkappa=histc(alphas_kappa(ENERGY_POPULATION),kappa_bin_pos-0.5*kappa_bin_size);
plot(kappa_bin_pos(1:end-1),Nkappa(1:end-1),'r');

title('\kappa');




figure(7)
set(gca,'FontSize',22);
hold on;
grid on;

delta_bin_size=0.02;
delta_bin_pos=(0:delta_bin_size:51*delta_bin_size);

Ndelta=histc(delta_orbit(ENERGY_POPULATION),delta_bin_pos-0.5*delta_bin_size);
Nravg=histc(r_avg(ENERGY_POPULATION),delta_bin_pos-0.5*delta_bin_size);
Ndiff=histc(r_avg(ENERGY_POPULATION)-2*delta_orbit(ENERGY_POPULATION),delta_bin_pos-10.5*delta_bin_size);

plot(delta_bin_pos(1:end-1),Ndelta(1:end-1),'r');
plot(delta_bin_pos(1:end-1),Nravg(1:end-1),'b--');
plot(delta_bin_pos(1:end-1)-10*delta_bin_size,Ndiff(1:end-1),'k--');
title('\delta and <r>');

sign_change_pop=zeros(Nalphas_simulated,1);

for number=1:Nalphas_simulated
    sign_vpll0=sign(vpll_output_pop(1,number));
    SIGN_CHANGE=find(0.5*abs(sign(vpll_output_pop(:,number))-sign_vpll0));
    if (~isempty(SIGN_CHANGE))
        sign_change_pop(number)=1;
    end
end

ALL_TRAPPED_POP=(sign_change_pop.*(~HFS_crossing));
ALL_PASSING_POP=(~sign_change_pop.*(HFS_crossing));
TRAPPED_CO_PASSING_POP=(sign_change_pop.*(~HFS_crossing).*(sign_v==1));
TRAPPED_COUNTER_PASSING_POP=(sign_change_pop.*(~HFS_crossing).*(sign_v==-1));
POTATOES_POP=(sign_change_pop.*(HFS_crossing));
STAGNATION_POP=((~sign_change_pop).*(~HFS_crossing).*(minX_value_pop>=X_axis));
ALL_PASSING_POP=ALL_PASSING_POP+((~sign_change_pop).*(~HFS_crossing).*(minX_value_pop<X_axis));
ALL_PASSING_POP(ALL_PASSING_POP>1)=1;
CO_PASSING_POP=((~sign_change_pop).*(HFS_crossing).*(sign_v==1));
CO_PASSING_POP=CO_PASSING_POP+((~sign_change_pop).*(~HFS_crossing).*(minX_value_pop<X_axis));
CO_PASSING_POP(CO_PASSING_POP>1)=1;
COUNTER_PASSING_POP=((~sign_change_pop).*(HFS_crossing).*(sign_v==-1));

ALL_TRAPPED=find(ALL_TRAPPED_POP);
ALL_PASSING=find(ALL_PASSING_POP);
TRAPPED_CO_PASSING=find(TRAPPED_CO_PASSING_POP);
TRAPPED_COUNTER_PASSING=find(TRAPPED_COUNTER_PASSING_POP);
POTATOES=find(POTATOES_POP);
STAGNATION=find(STAGNATION_POP);
CO_PASSING=find(CO_PASSING_POP);
COUNTER_PASSING=find(COUNTER_PASSING_POP);



disp('---------------------------------------------');
disp('-- GROUPS OF PARTICLE ORBITS (PERCENTAGES) --');
disp('---------------------------------------------');
% disp('EXOTIC_LIST=');
% disp(100*length(EXOTIC_LIST)/Nalphas_simulated);
disp('POTATOES_LIST=');
disp(100*length(POTATOES)/Nalphas_simulated);
disp('STAGNATION_LIST=');
disp(100*length(STAGNATION)/Nalphas_simulated);

disp('TRAPPED_COUNTER_PASSING=');
disp(100*length(TRAPPED_COUNTER_PASSING)/Nalphas_simulated);
disp('TRAPPED_CO_PASSING=');
disp(100*length(TRAPPED_CO_PASSING)/Nalphas_simulated);

disp('CO_PASSING=');
disp(100*length(CO_PASSING)/Nalphas_simulated);
disp('COUNTER_PASSING=');
disp(100*length(COUNTER_PASSING)/Nalphas_simulated);

disp('TOTAL=');
disp(100*(length(COUNTER_PASSING)+length(CO_PASSING)+length(TRAPPED_CO_PASSING)+length(TRAPPED_COUNTER_PASSING)+...
    length(STAGNATION)+length(POTATOES))/Nalphas_simulated);

% Now we need to scan the data
% to get an integer number of bounce period for the trapped
% so we need to mark the new end of data for precession 
% calculations


deriv_vpll_pop=zeros(end_ts,Nalphas_simulated);

for number=1:Nalphas_simulated
    deriv_vpll_pop(3:end-2,number)=0.25*(vpll_output_pop(5:end,number)-vpll_output_pop(1:end-4,number));
end


last_pos_pop=zeros(Nalphas_simulated,1);

for n=1:length(ALL_TRAPPED)
    number=ALL_TRAPPED(n);
%     x0_pos=X_output_pop(1,number);
%     z0_pos=Z_output_pop(1,number);
    v0_pos=vpll_output_pop(1,number);
%     dist_ini_pos=sqrt((x0_pos-X_output_pop(3:end,number)).^2+(z0_pos-Z_output_pop(3:end,number)).^2);
    dist_vpll_pos=abs(v0_pos-vpll_output_pop(3:end,number));
    min_dist=find(dist_vpll_pos(3:end)<0.05*v0_pos);
    if isempty(min_dist)
        min_dist=find(dist_vpll_pos<abs(0.1*v0_pos));
    end
    if isempty(min_dist)
        min_dist=end_ts-2;
    end
    last_period_pos=min_dist(end)+2;
%     [min_dist last_period_pos]=min(dist_ini_pos(end:-1:3));
    last_pos_pop(number)=last_period_pos;
end

% figure(1)
% hold on
% 
% for n=1:5000:length(ALL_TRAPPED)
%     number=ALL_TRAPPED(n);
%     plot(1:last_pos_pop(number),vpll_output_pop(1:last_pos_pop(number),number))
% end
