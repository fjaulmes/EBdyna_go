load('initial_NBI60keV_R_pre_collapse_all.mat')

Epll=0.5*(mHe/eV)*alphas_vpll.^2;
Eperp=max(alphas_Ekin-Epll,0);
alphas_vperp=sqrt(2*(eV/mHe)*Eperp);
alphas_vpll=alphas_vpll+R0*alphas_momentum;

alphas_vtot=sqrt(alphas_vperp.^2+alphas_vpll.^2);
alphas_Ekin=0.5*(mHe/eV)*alphas_vtot.^2;


Epll=0.5*(mHe/eV)*alphas_vpll.^2;
Eperp=max(alphas_Ekin-Epll,0);
alphas_vperp=sqrt(2*(eV/mHe)*Eperp);
alphas_vtot=sqrt(2*(eV/mHe)*alphas_Ekin);

% SAVENAME='initial_NBI60keV_Rlab_pre_collapse_all.mat'
% save(SAVENAME,'mHe','ZHe','alphas_momentum','pos_X_gc','pos_Z_gc','alphas_pos_x','alphas_pos_z','alphas_pos_phi','v_X','v_Z','v_phi','alphas_vpll','alphas_mm','alphas_psi','alphas_Ekin','alphas_pphi0','Nalphas_simulated');


%%
load('NBI60keV_R_fc1p6h1p6_all.mat')
% fixing rotation direction for this case
% rotation_profile_end=-interp1((0:256)/256,OMEGA_profile_interp_end,rho_tor_scale);
% alphas_pos_psi=interp2(scale_X,scale_Z,psi_norm_XZsmall_map',pos_X_gc,pos_Z_gc);
% alphas_momentum=interp1(1:Nradial,rotation_profile_end,alphas_pos_psi);

% to reflect the lowering of global momentum during the sawtooth.
ADJUST_MOMENTUM_FINAL=0.9
alphas_momentum=ADJUST_MOMENTUM_FINAL*alphas_momentum;

PART_POP=find(~alphas_ejected);

Epll=0.5*(mHe/eV)*alphas_vpll.^2;
Eperp=max(alphas_Ekin-Epll,0);
alphas_vperp=sqrt(2*(eV/mHe)*Eperp);
alphas_vpll=alphas_vpll+R0*alphas_momentum(1:length(alphas_vpll));

alphas_vtot=sqrt(alphas_vperp.^2+alphas_vpll.^2);
alphas_Ekin=0.5*(mHe/eV)*alphas_vtot.^2;


Epll=0.5*(mHe/eV)*alphas_vpll.^2;
Eperp=max(alphas_Ekin-Epll,0);
alphas_vperp=sqrt(2*(eV/mHe)*Eperp);
alphas_vtot=sqrt(2*(eV/mHe)*alphas_Ekin);
SAVENAME='NBI60keV_Rlab_fc1p6h1p6_all.mat'

save(SAVENAME,'mHe','ZHe','alphas_momentum','alphas_Ekin', 'alphas_mm', 'alphas_pphi0', 'alphas_pos_x', 'alphas_pos_z', 'alphas_pos_phi', 'alphas_vpll',  'alphas_ejected', 'alphas_psi_value_corr',  'time', 'frame_rank_precise', 'alphas_omega', 'alphas_Etot', 'pos_X_gc', 'pos_Z_gc');
