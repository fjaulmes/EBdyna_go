n=1;

for(p=1:NP)
    for(r=1:Nradial)
        %pressure_norm_PR_map(p,r)=finesse_data(n,end-31);
        pressure_norm_PR_map(p,r)=P_initial_profile(r);
        n=n+1;
    end
end
% Pnorm_XZ_map=griddata(finesse_data(:,1),finesse_data(:,2),finesse_data(:,end-31),XX,ZZ,'cubic');
% Pnorm_XZ_map(isnan(Pnorm_XZ_map)) = 0; 
% Pnorm_XZ_map=Pnorm_XZ_map';
% 
% Pnorm_XZ_map=Pnorm_XZ_map/max(max(Pnorm_XZ_map));
% 

gd_data=reshape((pressure_norm_PR_map(:,:)'),NP*Nradial,1);
P_XZ_map=griddata(finesse_data(:,1),finesse_data(:,2),gd_data,XX,ZZ,'cubic');
P_XZ_map(isnan(P_XZ_map))=0;
P_XZ_map=P_XZ_map';

Pnorm_XZ_map=P_XZ_map/max(max(P_XZ_map));

% Bpol_XZ_map=sqrt(BR_XZ_map.^2+BZ_XZ_map.^2);
Bpol_XZ_map=sqrt(BpolX_initial_map.^2+BpolZ_initial_map.^2);

%plasma_beta_recalc=mean(2*mu0*volume_flux_diff.*P_initial_profile)/mean(volume_flux_diff.*Btot_radial_profile.^2);
%plasma_beta_pol_recalc=mean(2*mu0*volume_flux_diff.*P_initial_profile)/mean(volume_flux_diff.*Bpol_radial_profile.^2);


B2_profile=radial_r_value_flux*0;
Bpol2_profile=radial_r_value_flux*0;
Ppoly_profile=radial_r_value_flux*0;
volume_profile=radial_r_value_flux*0;


plasma_beta_profile=radial_r_value_flux*0;
plasma_beta_pol_profile=radial_r_value_flux*0;

Pfinesse_profile=finesse_data(1:257,9);
%
% now more precisely

for (n=2:Nradial)
%     B2_profile(n)=B2_profile(n-1);
%     Bpol2_profile(n)=Bpol2_profile(n-1);
%     Ppoly_profile(n)=Ppoly_profile(n-1);
    %I_flux(n)=0;
    for(x=X1_Nradial(n):X2_Nradial(n))
        %Z_max=round(0.5*(Z_psi(n,x)+Z_psi_fit(n,x)))+mid_Z+1;
        Z_max_up=Z_psi_fit_up(n,x);
        Z_max_down=Z_psi_fit_down(n,x);
        % toroidal flux calculation
        for z=floor(Z_max_down/DX+mid_Z):ceil(Z_max_up/DX+mid_Z)
            B2_profile(n)=B2_profile(n)+(Bphi_XZ_map(x,z)^2)*(DX*DX*2*pi*((x-mid_X)*DX+R0));
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

plasma_beta_recalc=P0*plasma_beta_profile(end)
plasma_beta_pol_recalc=P0*plasma_beta_pol_profile(end)
