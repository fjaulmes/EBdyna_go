
% considering only co-passing particles for the moment

load('initialG_alphas_vA1_pre_collapse.mat')
load('initialG_alphas_vA1_precession_stats')

alphas_Ekin=real(alphas_Ekin);
alphas_mm=max(alphas_mm,0);

PART_POP_TAE=find((alphas_psi>=0.9*pTAE_inf).*(alphas_psi<=1.1*pTAE_sup).*CO_PASSING_POP.*(abs(alphas_vpll)<1e7));
alphas_pos_x1=alphas_pos_x(PART_POP_TAE);
alphas_pos_z1=alphas_pos_z(PART_POP_TAE);
alphas_pos_phi1=alphas_pos_phi(PART_POP_TAE);
alphas_Ekin1=alphas_Ekin(PART_POP_TAE);
alphas_mm1=alphas_mm(PART_POP_TAE);
alphas_vpll1=alphas_vpll(PART_POP_TAE);
alphas_psi1=alphas_psi(PART_POP_TAE);
v_X1=v_X(PART_POP_TAE);
v_Z1=v_Z(PART_POP_TAE);
v_phi1=v_phi(PART_POP_TAE);
alphas_pphi01=alphas_pphi0(PART_POP_TAE);

Nalphas_simulated1=length(alphas_vpll1)



load('initialG_alphas_vA2_pre_collapse.mat')
load('initialG_alphas_vA2_precession_stats')

alphas_Ekin=real(alphas_Ekin);
alphas_mm=max(alphas_mm,0);

PART_POP_TAE=find((alphas_psi>=0.9*pTAE_inf).*(alphas_psi<=1.1*pTAE_sup).*CO_PASSING_POP.*(abs(alphas_vpll)<1e7));
alphas_pos_x2=alphas_pos_x(PART_POP_TAE);
alphas_pos_z2=alphas_pos_z(PART_POP_TAE);
alphas_pos_phi2=alphas_pos_phi(PART_POP_TAE);
alphas_Ekin2=alphas_Ekin(PART_POP_TAE);
alphas_mm2=alphas_mm(PART_POP_TAE);
alphas_vpll2=alphas_vpll(PART_POP_TAE);
alphas_psi2=alphas_psi(PART_POP_TAE);
v_X2=v_X(PART_POP_TAE);
v_Z2=v_Z(PART_POP_TAE);
v_phi2=v_phi(PART_POP_TAE);
alphas_pphi02=alphas_pphi0(PART_POP_TAE);

Nalphas_simulated2=length(alphas_vpll2)



load('initialG_alphas_vA3_pre_collapse.mat')
load('initialG_alphas_vA3_precession_stats')

alphas_Ekin=real(alphas_Ekin);
alphas_mm=max(alphas_mm,0);

PART_POP_TAE=find((alphas_psi>=0.9*pTAE_inf).*(alphas_psi<=1.1*pTAE_sup).*CO_PASSING_POP.*(abs(alphas_vpll)<1e7));
alphas_pos_x3=alphas_pos_x(PART_POP_TAE);
alphas_pos_z3=alphas_pos_z(PART_POP_TAE);
alphas_pos_phi3=alphas_pos_phi(PART_POP_TAE);
alphas_Ekin3=alphas_Ekin(PART_POP_TAE);
alphas_mm3=alphas_mm(PART_POP_TAE);
alphas_vpll3=alphas_vpll(PART_POP_TAE);
alphas_psi3=alphas_psi(PART_POP_TAE);
v_X3=v_X(PART_POP_TAE);
v_Z3=v_Z(PART_POP_TAE);
v_phi3=v_phi(PART_POP_TAE);
alphas_pphi03=alphas_pphi0(PART_POP_TAE);

Nalphas_simulated3=length(alphas_vpll3)





load('initialG_alphas_vA4_pre_collapse.mat')
load('initialG_alphas_vA4_precession_stats')

alphas_Ekin=real(alphas_Ekin);
alphas_mm=max(alphas_mm,0);

PART_POP_TAE=find((alphas_psi>=0.9*pTAE_inf).*(alphas_psi<=1.1*pTAE_sup).*CO_PASSING_POP.*(abs(alphas_vpll)<1e7));
alphas_pos_x4=alphas_pos_x(PART_POP_TAE);
alphas_pos_z4=alphas_pos_z(PART_POP_TAE);
alphas_pos_phi4=alphas_pos_phi(PART_POP_TAE);
alphas_Ekin4=alphas_Ekin(PART_POP_TAE);
alphas_mm4=alphas_mm(PART_POP_TAE);
alphas_vpll4=alphas_vpll(PART_POP_TAE);
alphas_psi4=alphas_psi(PART_POP_TAE);
v_X4=v_X(PART_POP_TAE);
v_Z4=v_Z(PART_POP_TAE);
v_phi4=v_phi(PART_POP_TAE);
alphas_pphi04=alphas_pphi0(PART_POP_TAE);

Nalphas_simulated4=length(alphas_vpll4)




load('initialG_alphas_vA5_pre_collapse.mat')
load('initialG_alphas_vA5_precession_stats')

alphas_Ekin=real(alphas_Ekin);
alphas_mm=max(alphas_mm,0);

PART_POP_TAE=find((alphas_psi>=0.9*pTAE_inf).*(alphas_psi<=1.1*pTAE_sup).*CO_PASSING_POP.*(abs(alphas_vpll)<1e7));
alphas_pos_x5=alphas_pos_x(PART_POP_TAE);
alphas_pos_z5=alphas_pos_z(PART_POP_TAE);
alphas_pos_phi5=alphas_pos_phi(PART_POP_TAE);
alphas_Ekin5=alphas_Ekin(PART_POP_TAE);
alphas_mm5=alphas_mm(PART_POP_TAE);
alphas_vpll5=alphas_vpll(PART_POP_TAE);
alphas_psi5=alphas_psi(PART_POP_TAE);
v_X5=v_X(PART_POP_TAE);
v_Z5=v_Z(PART_POP_TAE);
v_phi5=v_phi(PART_POP_TAE);
alphas_pphi05=alphas_pphi0(PART_POP_TAE);

Nalphas_simulated5=length(alphas_vpll5)




load('initialG_alphas_vA6_pre_collapse.mat')
load('initialG_alphas_vA6_precession_stats')

alphas_Ekin=real(alphas_Ekin);
alphas_mm=max(alphas_mm,0);

PART_POP_TAE=find((alphas_psi>=0.9*pTAE_inf).*(alphas_psi<=1.1*pTAE_sup).*CO_PASSING_POP.*(abs(alphas_vpll)<1e7));
alphas_pos_x6=alphas_pos_x(PART_POP_TAE);
alphas_pos_z6=alphas_pos_z(PART_POP_TAE);
alphas_pos_phi6=alphas_pos_phi(PART_POP_TAE);
alphas_Ekin6=alphas_Ekin(PART_POP_TAE);
alphas_mm6=alphas_mm(PART_POP_TAE);
alphas_vpll6=alphas_vpll(PART_POP_TAE);
alphas_psi6=alphas_psi(PART_POP_TAE);
v_X6=v_X(PART_POP_TAE);
v_Z6=v_Z(PART_POP_TAE);
v_phi6=v_phi(PART_POP_TAE);
alphas_pphi06=alphas_pphi0(PART_POP_TAE);

Nalphas_simulated6=length(alphas_vpll6)




load('initialG_alphas_vA7_pre_collapse.mat')
load('initialG_alphas_vA7_precession_stats')

alphas_Ekin=real(alphas_Ekin);
alphas_mm=max(alphas_mm,0);

PART_POP_TAE=find((alphas_psi>=0.9*pTAE_inf).*(alphas_psi<=1.1*pTAE_sup).*CO_PASSING_POP.*(abs(alphas_vpll)<1e7));
alphas_pos_x7=alphas_pos_x(PART_POP_TAE);
alphas_pos_z7=alphas_pos_z(PART_POP_TAE);
alphas_pos_phi7=alphas_pos_phi(PART_POP_TAE);
alphas_Ekin7=alphas_Ekin(PART_POP_TAE);
alphas_mm7=alphas_mm(PART_POP_TAE);
alphas_vpll7=alphas_vpll(PART_POP_TAE);
alphas_psi7=alphas_psi(PART_POP_TAE);
v_X7=v_X(PART_POP_TAE);
v_Z7=v_Z(PART_POP_TAE);
v_phi7=v_phi(PART_POP_TAE);
alphas_pphi07=alphas_pphi0(PART_POP_TAE);

Nalphas_simulated7=length(alphas_vpll7)




load('initialG_alphas_vA8_pre_collapse.mat')
load('initialG_alphas_vA8_precession_stats')

alphas_Ekin=real(alphas_Ekin);
alphas_mm=max(alphas_mm,0);

PART_POP_TAE=find((alphas_psi>=0.9*pTAE_inf).*(alphas_psi<=1.1*pTAE_sup).*CO_PASSING_POP.*(abs(alphas_vpll)<1e7));
alphas_pos_x8=alphas_pos_x(PART_POP_TAE);
alphas_pos_z8=alphas_pos_z(PART_POP_TAE);
alphas_pos_phi8=alphas_pos_phi(PART_POP_TAE);
alphas_Ekin8=alphas_Ekin(PART_POP_TAE);
alphas_mm8=alphas_mm(PART_POP_TAE);
alphas_vpll8=alphas_vpll(PART_POP_TAE);
alphas_psi8=alphas_psi(PART_POP_TAE);
v_X8=v_X(PART_POP_TAE);
v_Z8=v_Z(PART_POP_TAE);
v_phi8=v_phi(PART_POP_TAE);
alphas_pphi08=alphas_pphi0(PART_POP_TAE);

Nalphas_simulated8=length(alphas_vpll8)





alphas_pos_x=[alphas_pos_x1 ; alphas_pos_x2 ; alphas_pos_x3 ; alphas_pos_x4 ; alphas_pos_x5 ; alphas_pos_x6 ; alphas_pos_x7 ; alphas_pos_x8];
alphas_pos_z=[alphas_pos_z1 ; alphas_pos_z2 ; alphas_pos_z3 ; alphas_pos_z4 ; alphas_pos_z5 ; alphas_pos_z6 ; alphas_pos_z7 ; alphas_pos_z8];
alphas_pos_phi=[alphas_pos_phi1 ; alphas_pos_phi2 ; alphas_pos_phi3 ; alphas_pos_phi4 ; alphas_pos_phi5 ; alphas_pos_phi6 ; alphas_pos_phi7 ; alphas_pos_phi8];

alphas_Ekin=[alphas_Ekin1 ; alphas_Ekin2 ; alphas_Ekin3 ; alphas_Ekin4 ; alphas_Ekin5 ; alphas_Ekin6 ; alphas_Ekin7 ; alphas_Ekin8];
alphas_mm=[alphas_mm1 ; alphas_mm2 ; alphas_mm3 ; alphas_mm4 ; alphas_mm5 ; alphas_mm6 ; alphas_mm7 ; alphas_mm8];
alphas_psi=[alphas_psi1 ; alphas_psi2 ; alphas_psi3 ; alphas_psi4 ; alphas_psi5 ; alphas_psi6 ; alphas_psi7 ; alphas_psi8];
alphas_vpll=[alphas_vpll1 ; alphas_vpll2 ; alphas_vpll3 ; alphas_vpll4 ; alphas_vpll5 ; alphas_vpll6 ; alphas_vpll7 ; alphas_vpll8];

v_X=[v_X1 ; v_X2 ; v_X3 ; v_X4 ; v_X5 ; v_X6 ; v_X7 ; v_X8];
v_Z=[v_Z1 ; v_Z2 ; v_Z3 ; v_Z4 ; v_Z5 ; v_Z6 ; v_Z7 ; v_Z8];
v_phi=[v_phi1 ; v_phi2 ; v_phi3 ; v_phi4 ; v_phi5 ; v_phi6 ; v_phi7 ; v_phi8];

alphas_pphi0=[alphas_pphi01 ; alphas_pphi02 ; alphas_pphi03 ; alphas_pphi04 ; alphas_pphi05 ; alphas_pphi06 ; alphas_pphi07 ; alphas_pphi08];


Nalphas_simulated=Nalphas_simulated1+Nalphas_simulated2+Nalphas_simulated3+Nalphas_simulated4+Nalphas_simulated5+Nalphas_simulated6+Nalphas_simulated7+Nalphas_simulated8

save('initialG_alphas_TAE_large.mat', 'mHe', 'ZHe', 'density_part_ratio', 'Nalphas_simulated', ...
'alphas_pos_x', 'alphas_pos_z', 'alphas_pos_phi', 'alphas_Ekin', ...
'alphas_mm', 'alphas_vpll', 'alphas_psi', 'v_X', 'v_Z', 'v_phi', 'alphas_pphi0');
