

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
if SIGN_TOROIDAL_FIELD<0
Bpol_XZ_map=Bpol_XZ_map.*(psi_XZ_map>0);
Btor_XZ_map=Bphi_XZ_map.*(psi_XZ_map>0);
Pnorm_XZ_map=Pnorm_XZ_map.*(psi_XZ_map>0);
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

% BETA_ALPHAS=0.2
TE_FRAC=1/(2+BETA_ALPHAS)

pol_P_reverse=pol_P(end:-1:1);
pol_N_reverse=pol_N(end:-1:1);
pol_T_reverse=pol_T(end:-1:1);

% psi_range=((1:Nradial)-1)/256;
if SIGN_TOROIDAL_FIELD==1
	psi_range=1-psi_scale/min(psi_scale);
else
	psi_range=1-psi_scale/max(psi_scale);
end

Te_profile=polyval(pol_T_reverse,psi_range)*Te0;
Pinput_profile=polyval(pol_P_reverse,psi_range)*P0;
Ne_profile=polyval(pol_N_reverse,psi_range)*Ne0;

%
%
% Te_profile=0.5*Pinput_profile./Ne_profile;
% Te_profile=(polyval(pol_N,psi_range)*Te0)/eV;
% Te_profile=0.5*Pinput_profile./Ne_profile;
Te_profile=0.5*Te_profile+0.5*TE_FRAC*P_initial_profile./Ne_profile;

% Ne_profile=0.5*((0.5)*P_initial_profile./Te_profile+Ne_profile);
Ne_profile=TE_FRAC*P_initial_profile./Te_profile;
[Ne_max Ne_max_pos]=max(Ne_profile);
Ne_profile(1:Ne_max_pos)=max(Ne_profile(1:Ne_max_pos),Ne_max);
Te_profile=TE_FRAC*P_initial_profile./Ne_profile;
% Te_profile=0.5*(0.5*Pinput_profile./Ne_profile+Te_profile);

[Te_min Te_min_pos]=min(Te_profile(1:end-2));
Te_profile(Te_min_pos:end)=min(Te_profile(Te_min_pos:end),Te_min);
[Te_max Te_max_pos]=max(Te_profile(1:end-2));
Te_profile(1:Te_max_pos)=max(Te_profile(1:Te_max_pos),Te_max);
Ne_profile=TE_FRAC*P_initial_profile./Te_profile;

Ne_profile=max(Ne_profile,0);

RESCALING_TEMPERATURE_FACTOR=(1.0)
Te_min_pos=round(1.1*Te_min_pos)
RESCALING_TEMPERATURE_FACTOR_END=(Nradial:-1:Te_min_pos)/Nradial

Te_profile(Te_min_pos:end)=(RESCALING_TEMPERATURE_FACTOR_END).*Te_profile(Te_min_pos:end);
Ne_profile(Te_min_pos:end)=(1./RESCALING_TEMPERATURE_FACTOR_END).*Ne_profile(Te_min_pos:end);
Te_profile=RESCALING_TEMPERATURE_FACTOR*Te_profile;
Ne_profile=(1./RESCALING_TEMPERATURE_FACTOR)*Ne_profile;

P_initial_profile=(2+BETA_ALPHAS)*Te_profile.*Ne_profile;

Ne0=Ne_profile(1)
Te0=Te_profile(1);



disp('core temperature = ')
disp(Te0/eV)
disp('in eV')
