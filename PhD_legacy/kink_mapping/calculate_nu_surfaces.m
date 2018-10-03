clear S3 S3_filt;


%%
S3=zeros(1,N_NU);

for (n=1:N_NU)
    %%
    nu_test=nu_values(n);
    clear nu_surf_radius;
    r_count=1;
    nu_surf_radius=zeros(Nomega*N_NU,1);
    
    for (omega=1:Nomega)    
        for (r=2:pos_psi_rx)
            r_value=scale_r(r);
            omega_value=scale_omega(omega);
            if ((nu(r,omega)>(nu_test-Dnu)) && (nu(r,omega)<=(nu_test)));
                Dr1=Dr_avg(r);
                if (omega==1) || (omega==Nomega)
                    %nu_surf_radius(r_count)=sqrtg_S_avg(r)*Domega*Dr1;
                    nu_surf_radius(r_count)=dl_avg(r)*Dr1;
                else
                    Dr2=Dr1;
                    %nu_surf_radius(r_count)=sqrtg_S_avg(r)*Domega*Dr1+sqrtg_S_avg(r)*Domega*Dr2;
                    nu_surf_radius(r_count)=dl_avg(r)*Dr1+dl_avg(r)*Dr2;
                end
                r_count=r_count+1;
            end
        end
    end
    S3(n)=sum(nu_surf_radius);
end
%%
S3_filt=S3;

% filter profile of surfaces to avoid numerical noise
S3_filt(1)=0;
%  for n=2:N_NU-1
%      S3_filt(n)=S3_filt(n-1)+S3(n);
%  end
for n=2:N_NU-1
  S3_filt(n)=S3_filt(n-1)+S3(n);
end
S3_filt(N_NU)=S3_filt(N_NU-1)+S3(N_NU);

%%
S3=S3_filt;


S3_pol=polyfit(nu_values(1:end),S3(1:end),1);

S3_filt=polyval(S3_pol,nu_values);


S3_min=S3_filt(1);
S3=S3_filt-S3_min;

S3_max=S3(end);


S3=(S3/S3_max)*area_region3;


psi3_nu=interp1(surf_flux_precise(1:NR),Psih(1:NR),S3,'cubic');


