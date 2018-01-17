

% non ejected particles only

Nalphas_simulated=length(vpll_output_pop)

% clear posX_output posZ_output 


sign_v=zeros(Nalphas_simulated,1);
max_vpll_avg=zeros(Nalphas_simulated,1);
max_vpll_traj=zeros(Nalphas_simulated,1);
min_vpll_traj=zeros(Nalphas_simulated,1);
max_vpll_pos=zeros(Nalphas_simulated,1);
min_vpll_neg=zeros(Nalphas_simulated,1);
for number=1:Nalphas_simulated
    [max_vpll max_index]=max(abs(vpll_output_pop(:,number)));
    max_vpll_avg(number)=vpll_output_pop(max_index,number);
    max_vpll_traj(number)=max(vpll_output_pop(:,number));
	min_vpll_traj(number)=min(vpll_output_pop(:,number));
%     sign_v(number)=sign(vpll_output_corr(max_index,number));
end
sign_v=sign(max_vpll_avg);
% min_vpll_traj=min_vpll_traj';
max_vpll_pos=round(0.5*((sign(max_vpll_traj))+1));
min_vpll_neg=round(0.5*((-sign(min_vpll_traj))+1));
sigma_vpll=((max_vpll_pos+min_vpll_neg)==2);

calculate_vD_XZ_map_post;
calculate_gradZ_theta_map;
calculate_gradX_theta_map;


theta_output_pop=interp2(scale_X,scale_Z,theta_XZsmall_map',X_output_pop,Z_output_pop,'*nearest');
GZ_theta_output_pop=interp2(scale_X,scale_Z,gradZ_theta_map',X_output_pop,Z_output_pop,'*linear');
GX_theta_output_pop=interp2(scale_X,scale_Z,gradX_theta_map',X_output_pop,Z_output_pop,'*linear');

q_output_pop=interp2(scale_X,scale_Z,q_final_XZsmall_map',X_output_pop,Z_output_pop,'*linear');
r_output_pop=interp2(scale_X,scale_Z,psi_norm_XZsmall_map',X_output_pop,Z_output_pop,'*linear');
r_output_pop=interp1(1:257,radial_r_value_flux,r_output_pop,'*linear');

% r_output_pop=sqrt((X_output_pop-X_axis).^2+Z_output_pop.^2);



% r_output=interp1(1:257,radial_r_value_flux,r_output,'*linear');
vDX_output_corr=interp2(scale_X,scale_Z,vD_X_XZ_map',X_output_pop,Z_output_pop,'*linear');
vDZ_output_corr=interp2(scale_X,scale_Z,vD_Z_XZ_map',X_output_pop,Z_output_pop,'*linear');
vDphi_output_corr=interp2(scale_X,scale_Z,vD_phi_XZ_map',X_output_pop,Z_output_pop,'*linear');


clear vD_phi_XZ_map vD_Z_XZ_map psi_norm_XZsmall_map q_final_XZsmall_map gradZ_theta_map theta_XZsmall_map

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

% lambda_co_passing=Bmax_co_passing.*alphas_mm./(alphas_Ekin);
% lambda_counter_passing=Bmax_counter_passing.*alphas_mm./(alphas_Ekin);
% lambda_co_passing=Bavg*R0*alphas_mm./(alphas_Ekin.*(R0-r_avg+delta_orbit));
% lambda_counter_passing=Bavg*R0*alphas_mm./(alphas_Ekin.*(R0-r_avg-delta_orbit));
PLUS=find(sign_v==1);
MINUS=find(sign_v==-1);
HFS_crossing=zeros(Nalphas_simulated,1);
sigma_LFS=zeros(Nalphas_simulated,1);
delta_theta_pi_pop=zeros(Nalphas_simulated,1);
minX_value_pop=zeros(Nalphas_simulated,1);
maxX_value_pop=zeros(Nalphas_simulated,1);
% theta_output_pop=theta_output(:,PART_POP);

clear theta_output

for number=1:Nalphas_simulated
    [HFS_theta_value HFS_theta_pos]=min(abs(theta_output_pop(:,number)-pi));
    [minX_value minX_pos]=min(X_output_pop(:,number));
    [maxX_value maxX_pos]=max(X_output_pop(:,number));
    minX_value_pop(number)=minX_value;
    maxX_value_pop(number)=maxX_value;
    delta_theta_pi_pop(number)=HFS_theta_value;
end
HFS_crossing=((minX_value_pop<=X_axis).*(delta_theta_pi_pop<=(0.045-0.045*r_avg/a)));
sigma_LFS=(maxX_value_pop>=X_axis);
not_sigma_LFS=(maxX_value_pop<=X_axis);

clear delta_theta_pi_pop 

alphas_lambda=Bavg*alphas_mm./alphas_Ekin;






disp('EVALUATED_POPULATION')
disp(length(alphas_Ekin));


%alphas_pos_psi=interp2(scale_X,scale_Z,psi_norm_XZsmall_map',alphas_pos_x,alphas_pos_z,'*linear');
%alphas_value_psi=interp2(scale_X,scale_Z,psi_XZsmall_map',alphas_pos_x,alphas_pos_z,'*linear');
%theta_value=interp2(scale_X,scale_Z,theta_XZsmall_map',alphas_pos_x,alphas_pos_z,'*nearest');


time_scale_corr=time_scale';
end_ts=size(vpll_output_pop,1);

B_output_pop=interp2(scale_X,scale_Z,Btot_XZ_map',X_output_pop,Z_output_pop,'*linear');


for ts=1:end_ts
    vDX_output_corr(ts,:)=(vDX_output_corr(ts,:)).*(2*alphas_Ekin'-B_output_pop(ts,:).*alphas_mm')/(ZHe);
    vDZ_output_corr(ts,:)=(vDZ_output_corr(ts,:)).*(2*alphas_Ekin'-B_output_pop(ts,:).*alphas_mm')/(ZHe);
    vDphi_output_corr(ts,:)=(vDphi_output_corr(ts,:)).*(2*alphas_Ekin'-B_output_pop(ts,:).*alphas_mm')/(ZHe);
end

%clear   vDZ_output vDphi_output omega_output

B_avg=mean(B_output_pop(:,:),1)';
clear B_output_pop;



%sign_change_pop=zeros(Nalphas_simulated,1);

%for number=1:Nalphas_simulated
%    sign_vpll0=sign(vpll_output_pop(1,number));
%    SIGN_CHANGE=find(0.5*abs(sign(vpll_output_pop(:,number))-sign_vpll0));
%    if (~isempty(SIGN_CHANGE))
%        sign_change_pop(number)=1;
%   end
%end
sign_change_pop=sigma_vpll;

nb_STAGNATION_HFS=length(find((~sign_change_pop).*(~sigma_LFS).*(HFS_crossing)))

ALL_TRAPPED_POP=(sign_change_pop.*(~HFS_crossing));
ALL_PASSING_POP=(~sign_change_pop.*(HFS_crossing).*(sigma_LFS));
TRAPPED_MINUS_POP=(sign_change_pop.*(~HFS_crossing).*(sign_v==-1));
TRAPPED_PLUS_POP=(sign_change_pop.*(~HFS_crossing).*(sign_v==1));
POTATOES_POP=(sign_change_pop.*(HFS_crossing).*(sigma_LFS));
STAGNATION_POP=((~sign_change_pop).*(minX_value_pop>=X_axis).*sigma_LFS)+((~sign_change_pop).*(not_sigma_LFS));
% STAGNATION_POP(STAGNATION_POP>1)=1;
ALL_PASSING_POP=ALL_PASSING_POP+((~sign_change_pop).*(~HFS_crossing).*(minX_value_pop<X_axis));
ALL_PASSING_POP(ALL_PASSING_POP>1)=1;
CO_PASSING_POP=((~sign_change_pop).*(HFS_crossing).*(sign_v==1).*(sigma_LFS));
CO_PASSING_POP=CO_PASSING_POP+((sigma_LFS).*(~sign_change_pop).*(~HFS_crossing).*(minX_value_pop<X_axis).*(sign_v==1));
CO_PASSING_POP(CO_PASSING_POP>1)=1;
COUNTER_PASSING_POP=((~sign_change_pop).*(HFS_crossing).*(sigma_LFS).*(sign_v==-1));
COUNTER_PASSING_POP=COUNTER_PASSING_POP+((sigma_LFS).*(~sign_change_pop).*(~HFS_crossing).*(minX_value_pop<X_axis).*(sign_v==-1));
COUNTER_PASSING_POP(COUNTER_PASSING_POP>1)=1;
ALL_PASSING_POP=CO_PASSING_POP+COUNTER_PASSING_POP;
ALL_PASSING_POP(ALL_PASSING_POP>1)=1;

ALL_TRAPPED=find(ALL_TRAPPED_POP);
ALL_PASSING=find(ALL_PASSING_POP);
TRAPPED_PLUS=find(TRAPPED_PLUS_POP);
TRAPPED_MINUS=find(TRAPPED_MINUS_POP);
POTATOES=find(POTATOES_POP);
STAGNATION=find(STAGNATION_POP);
CO_PASSING=find(CO_PASSING_POP);
COUNTER_PASSING=find(COUNTER_PASSING_POP);

% clear minX_value_pop



disp('---------------------------------------------');
disp('-- GROUPS OF PARTICLE ORBITS (PERCENTAGES) --');
disp('---------------------------------------------');
% disp('EXOTIC_LIST=');
% disp(100*length(EXOTIC_LIST)/Nalphas_simulated);
disp('POTATOES_LIST=');
disp(100*length(POTATOES)/Nalphas_simulated);
disp('STAGNATION_LIST=');
disp(100*length(STAGNATION)/Nalphas_simulated);

disp('TRAPPED_MINUS=');
disp(100*length(TRAPPED_MINUS)/Nalphas_simulated);
disp('TRAPPED_PLUS=');
disp(100*length(TRAPPED_PLUS)/Nalphas_simulated);

disp('CO_PASSING=');
disp(100*length(CO_PASSING)/Nalphas_simulated);
disp('COUNTER_PASSING=');
disp(100*length(COUNTER_PASSING)/Nalphas_simulated);

disp('TOTAL=');
disp(100*(length(COUNTER_PASSING)+length(CO_PASSING)+length(TRAPPED_PLUS)+length(TRAPPED_MINUS)+...
    length(STAGNATION)+length(POTATOES))/Nalphas_simulated);
	


calculate_precession_distribution_post

save_distribution_precession_data_post

