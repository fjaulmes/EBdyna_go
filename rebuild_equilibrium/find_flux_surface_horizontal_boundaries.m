% *************************************************************
% Discretisation to a finite number of flux values
% *************************************************************


X1_Nradial=zeros(1,Nradial);
X2_Nradial=zeros(1,Nradial);


psi2D_mask=psi2D.*mask_XZ_map;

% finding what the inner (left) and outer (right)
% boundaries for the external flux surface
    if SIGN_CO_CURRENT_FIELD>0
		[eps x_inner ] = min(mask_XZ_map(mid_X:-1:1,mid_Z));
	else
		[eps x_inner ] = min(mask_XZ_map(mid_X:-1:1,mid_Z));
	end
x_inner=mid_X-x_inner+1;
    if SIGN_CO_CURRENT_FIELD>0
		[eps x_outer ] = min(mask_XZ_map(mid_X:NR,Z_axis_pos));
	else
		[eps x_outer ] = min(mask_XZ_map(mid_X:NR,Z_axis_pos));
	end
x_outer=x_outer+mid_X-1;


% now building mapping in psi
% of the inner (X1) and outer (X2) boundaries

Nradial_min=12;
NP_half=round(NP/2);

X1_Nradial(1)=Raxis_pos;
X2_Nradial(1)=Raxis_pos;
for (n=2:Nradial-1)

    [eps X1_Nradial(n)]=min(abs(psi2D_mask(x_inner-1:Raxis_pos-1,Z_axis_pos)-psi_profile(n))); 
    [eps X2_Nradial(n)]=min(abs(psi2D_mask(Raxis_pos+1:x_outer+1,Z_axis_pos)-psi_profile(n))); 
    X1_Nradial(n)=x_inner+X1_Nradial(n)-2;
    X2_Nradial(n)=X2_Nradial(n)+Raxis_pos;
end
X1_Nradial(Nradial)=x_inner;
X2_Nradial(Nradial)=x_outer;


% *************************************************************
% Adjust curves for profiles with Nradial values
% *************************************************************

X1_Nradial(1)=Raxis_pos;
X2_Nradial(1)=Raxis_pos;


X1_Nradial=round(X1_Nradial);
X2_Nradial=round(X2_Nradial);


