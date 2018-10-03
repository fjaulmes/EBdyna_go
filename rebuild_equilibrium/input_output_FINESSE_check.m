load('../data_tokamak/physics_constants.mat')
load('../data_tokamak/pressure_profile.mat')
load('../data_tokamak/B_fields.mat')
load('../data_tokamak/flux_geometry.mat')
load('../data_tokamak/tokamak_map_dimensions')
load('../data_tokamak/tokamak_PR_map')
load('../finesse_tokamak_parameters')


%%
psi_norm_scale=-psi_scale/psi_scale(1)+1;

F_PR_map=(Btor_PR_map.*Rpos_PR_map);
F_profile=mean(F_PR_map(1:NP-1,:),1);
F2_profile=F_profile.^2;
F2_profile_shift=F2_profile-F2_profile(1);


% pol_F2=[0. -3.3  6.9  -8.4 3.8  ]
% pol_P=[1. -1 -4.7 11 -6.5 0.25  ]
% pol_N=[ 1. -1.2   8.9  -32.8  45.7 -21.5 ]


pol_F2(1:end)=pol_F2(end:-1:1);
pol_P(1:end)=pol_P(end:-1:1);
pol_N(1:end)=pol_N(end:-1:1);

psi_range=((1:Nradial)-1)/256;

F2_input_profile=polyval(pol_F2,psi_range);
P_input_profile=polyval(pol_P,psi_range)*P0;
Ne_input_profile=polyval(pol_N,psi_range)*Ne0;

F2_input_profile=F2_input_profile-F2_input_profile(1);

figure(1);
grid on;
hold on;
set(gca,'fontsize',20);
plot(psi_range,F2_input_profile/max(abs(F2_input_profile)),'linewidth',3);
plot(psi_norm_scale,F2_profile_shift/max(abs(F2_profile_shift)),'r--','linewidth',2);

figure(2);
grid on;
hold on;
set(gca,'fontsize',20);
plot(psi_range,P_input_profile/max(abs(P_input_profile)),'linewidth',3);
plot(psi_norm_scale,P_initial_profile/max(abs(P_initial_profile)),'r--','linewidth',2);
