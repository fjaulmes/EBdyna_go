% ***************************************************************
% Define polar coordinates centered on the shifted magnetic axis
% ki gives an approximate value of the poloidal angle
% ***************************************************************
%
clear theta_scale




DP=(2*pi)/(NP-1);
for(p=1:NP)
    theta_scale(p)=(p-1)*DP;    
end


% now ki angle mapping

ki_XZ_map=zeros(length(X_scale),length(Z_scale));

for (x=1:length(X_scale))
    for (z=1:length(Z_scale))

            xi_value=X_axis;
            r_value=radial_XZ_map(x,z);
            ki_XZ_map(x,z)=atan(abs(Z_scale(z)-Z_axis)/abs(X_scale(x)-xi_value));
            if (X_scale(x)<X_axis)
                ki_XZ_map(x,z)=pi-ki_XZ_map(x,z);
            end
            if (Z_scale(z)<Z_axis)
                ki_XZ_map(x,z)=2*pi-ki_XZ_map(x,z);
            end            
           
    end
end


