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

% S3_pol=polyfit(nu_values(1:end),S3(1:end),1);
% 
% S3_filt=polyval(S3_pol,nu_values);

     
if (f>8*PRECISE_MESH)
     S3_pol=polyfit(nu_values(1:end),S3(1:end),3);
    
     S3_filt=polyval(S3_pol,nu_values);
% elseif (f>=2)
%      S3_pol=polyfit(nu_values(1:end),S3(1:end),2);
%     
%      S3_filt=polyval(S3_pol,nu_values);
%      %S3_filt=interp1(nu_values(1:2:end),S3(1:2:end),nu_values,'cubic');
else
     S3_pol=polyfit(nu_values(1:end),S3(1:end),1);
    
     S3_filt=polyval(S3_pol,nu_values);
end
% S3_filt=interp1(nu_values(1:5:end),S3(1:5:end),nu_values,'cubic');
S3_min=S3_filt(1);
% S3=S3_filt;
S3=S3_filt-S3_min;

S3_max=S3(end);
% if f<=f_transition
% 	if f<round(0.55*f_final)
% 		S3_max=S3_max*(1.0001-0.00001*f)
% 	elseif f<=f_final
% 	%     S3_max=S3_max*max((1.0009+0.0001*f),1.0011)
% 		S3_max=S3_max*max((1.0001-0.00001*round(0.55*f_final))+0.000015*(f-round(0.6*f_final)),1.0000)
% 	else 
% 		S3_max=S3_max*1.000
% 	end
% end
%r_psi_final_rx=interp1((1:pos_psi_rx),scale_r(1:pos_psi_rx),x_psi_final_rx);

%S3=(S3/S3_max)*interp1((1:size_r),surf_flux(1:size_r),x_psi_final_rx,'cubic');



S3=(S3/S3_max)*area_region3;
%S3=S3+S3_min;
%S3=area_region3/S3_max;

%L3=sqrt(S3/pi);

%%

psi3_nu=interp1(surf_flux_precise(1:NR),Psih_final(1:NR),S3,'pchip');
Dr_avg_pos_psi_rx=interp1(1:NR,Dr_avg(1:NR),pos_psi_rx,'pchip');

if (f<FRAME_LIMIT_SUP_FOR_CONTINUITY)
    der_Psih_rx=interp1(1:NR-1,dPsih(1:NR-1),x_nu_cont13,'pchip');
    
    psi3_nu_max=max(psi3_nu);
    psi3_nu_min=min(psi3_nu);
    psi3_nu_norm=psi3_nu-psi3_nu_min;
    psi3_nu_norm=psi3_nu_norm/max(psi3_nu_norm);
    
    psi3_nu_max=psi3_nu_max-(psi_rx+0.001*der_Psih_rx*Dr_avg_pos_psi_rx);
    psi3_nu=psi3_nu_max*psi3_nu_norm+psi_rx+0.001*der_Psih_rx*Dr_avg_pos_psi_rx;
end
%psi3_nu_max=max(psi3_nu);
%psi3_nu=((psi3_nu-psi3_nu_min)/(psi3_nu_max-psi3_nu_min))*(max_psi3-min_psi3)+min_psi3;


