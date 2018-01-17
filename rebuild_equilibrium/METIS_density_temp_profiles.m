

%refine the pressure from the value of the betas

P_finesse_profile=finesse_data(1:Nradial,end-31);
P_finesse_profile=P_finesse_profile';

Pnorm_XZ_map=griddata(finesse_data(:,1),finesse_data(:,2),finesse_data(:,end-31),XX,ZZ,'cubic');
Pnorm_XZ_map(isnan(Pnorm_XZ_map)) = 0; 
Pnorm_XZ_map=Pnorm_XZ_map';

Pnorm_XZ_map=Pnorm_XZ_map/max(max(Pnorm_XZ_map));

% 
BpolX_initial_map=BR_XZ_map;
BpolZ_initial_map=BZ_XZ_map;

Bpol_XZ_map=sqrt(BpolX_initial_map.^2+BpolZ_initial_map.^2);
if psi_profile(1)>0
Bpol_XZ_map=Bpol_XZ_map.*(psi_XZ_map>=0);
Btor_XZ_map=Bphi_XZ_map.*(psi_XZ_map>=0);
Pnorm_XZ_map=Pnorm_XZ_map.*(psi_XZ_map>=0);
else
Bpol_XZ_map=Bpol_XZ_map.*(psi_XZ_map<=0);
Btor_XZ_map=Bphi_XZ_map.*(psi_XZ_map<=0);
Pnorm_XZ_map=Pnorm_XZ_map.*(psi_XZ_map<=0);
end

B2_profile=radial_r_value_flux*0;
Bpol2_profile=radial_r_value_flux*0;
Ppoly_profile=radial_r_value_flux*0;
volume_profile=radial_r_value_flux*0;

plasma_beta_profile=radial_r_value_flux*0;
plasma_beta_pol_profile=radial_r_value_flux*0;

for (n=2:Nradial)
    for(x=X1_Nradial(n):X2_Nradial(n))
        Z_max_up=Z_psi_fit_up(n,x);
        Z_max_down=Z_psi_fit_down(n,x);
        % surface integral calculation
        for z=floor(Z_max_down/DX+mid_Z):ceil(Z_max_up/DX+mid_Z)
            B2_profile(n)=B2_profile(n)+(Btor_XZ_map(x,z)^2)*(DX*DX*2*pi*((x-mid_X)*DX+R0));
            Bpol2_profile(n)=Bpol2_profile(n)+(Bpol_XZ_map(x,z)^2)*(DX*DX*2*pi*((x-mid_X)*DX+R0));
            Ppoly_profile(n)=Ppoly_profile(n)+Pnorm_XZ_map(x,z)*(DX*DX*2*pi*((x-mid_X)*DX+R0));
            volume_profile(n)=volume_profile(n)+(DX*DX*2*pi*((x-mid_X)*DX+R0));
        end
    end
    plasma_beta_profile(n)=2*mu0*Ppoly_profile(n)/B2_profile(n);
    plasma_beta_pol_profile(n)=2*mu0*Ppoly_profile(n)/Bpol2_profile(n);
end

B2_profile=B2_profile/volume_profile(end);
Bpol2_profile=Bpol2_profile/volume_profile(end);
Ppoly_profile=Ppoly_profile/volume_profile(end);

plasma_beta_recalc=plasma_beta_profile(end)
plasma_beta_pol_recalc=plasma_beta_pol_profile(end)

P0_recalc=0.5*(plasma_beta_tor/plasma_beta_recalc+plasma_beta_pol/plasma_beta_pol_recalc)


disp('P0_recalc is used to scale the final pressure profile....')

% refresh the data with the new pressure value
P_initial_profile=P_finesse_profile*P0_recalc/P_finesse_profile(1);
P0=P0_recalc;


%%

if psi_profile(1)<0
    rho_pol_profile=sqrt(1+psi_profile/abs(psi_profile(1)));
else
    rho_pol_profile=sqrt(psi_profile/abs(psi_profile(1)));
end
load('../METIS_profiles.mat');
Ne_profile=interp1(METISdata.rho_pol_metis, METISdata.ne_profile,rho_pol_profile);
Te_profile=interp1(METISdata.rho_pol_metis, METISdata.Te_profile,rho_pol_profile);
Ni_profile=interp1(METISdata.rho_pol_metis, METISdata.ni_profile,rho_pol_profile);
Ti_profile=interp1(METISdata.rho_pol_metis, METISdata.Ti_profile,rho_pol_profile);
Ptot_profile=interp1(METISdata.rho_pol_metis, METISdata.Ptot_profile,rho_pol_profile);
Pbulk_profile=Ne_profile.*Te_profile*eV+Ni_profile.*Ti_profile*eV;

P0_METIS=Ptot_profile(1)

% This new Pfast_profile is based on the equilibrium solution for P
% and NOT the METIS data : we use this difference to get an exact match
% for total pressure from GS solution
Pfast_profile=Ptot_profile-Pbulk_profile;
disp('ratio of fast particle presure over electron pressure:')
BETA_ALPHAS_profile=Pfast_profile./(Te_profile.*Ne_profile*eV);
BETA_ALPHAS=mean(BETA_ALPHAS_profile)

Ne0=Ne_profile(1)
Te0=Te_profile(1)*eV;



disp('core temperature = ')
disp(Te0/eV)
disp('in eV')
