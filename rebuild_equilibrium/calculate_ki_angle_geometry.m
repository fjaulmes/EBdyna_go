% ***************************************************************
% Define polar coordinates centered on the shifted magnetic axis
% ki gives an approximate value of the poloidal angle
% ***************************************************************
%
clear theta_scale


% % also toroidal flux calculation
% 
% tor_flux_profile_der=q_profile.*gradient(psi_scale,radial_r_value_flux);
% tor_flux_profile=tor_flux_profile_der*0;
% for (r=2:Nradial)
%     tor_flux_profile(r)=tor_flux_profile(r-1)+tor_flux_profile_der(r)*(surf_flux(r)-surf_flux(r-1));
% end


DP=(2*pi)/(NP-1);
for(p=1:NP)
    theta_scale(p)=(p-1)*DP;    
end


% now ki angle mapping

ki_XZ_map=zeros(NZ,NZ);

for (x=inf_X:sup_X)
    for (z=inf_Z:sup_Z)
        if(r_radial_mask_map(x,z)==1)
            xi_value=xi_map(x,z);
            r_value=r_value_map(x,z);
            if (z>=mid_Z)
                ki_XZ_map(x,z)=acos((X_scale(x)-xi_value)/r_value);
            else
                ki_XZ_map(x,z)=2*pi-acos((X_scale(x)-xi_value)/r_value);
            end
            
        end              
    end
end

ki_XZ_map(Raxis_pos,mid_Z)=0;

NKI=NP;
Dki=2*pi/(NKI-1);

ki_values=((1:NKI)-1)*Dki;




ki_data=reshape(ki_XZ_map',NZ*NZ,1);
X_data=reshape(XX,NZ*NZ,1);
Z_data=reshape(ZZ,NZ*NZ,1);

ki_PR_map=griddata(X_data,Z_data,ki_data,X_PR_map,Z_PR_map,'cubic');
ki_PR_map=max(ki_PR_map,0);
ki_PR_map=min(ki_PR_map,2*pi);
ki_PR_map(:,1)=zeros(NP,1);
ki_PR_map(:,Nradial)=zeros(NP,1)+pi;



for (r=1:Nradial)
    for (k=1:NP)
        theta_ki_psi_map(k,r)=interp1(theta_scale,ki_PR_map(:,r),ki_values(k));
    end
end
  

% theta_data=reshape(theta_XZ_map',NZ*NZ,1);


    