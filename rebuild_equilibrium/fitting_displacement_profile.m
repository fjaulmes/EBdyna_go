% *************************************************************
% Discretisation to a finite number of flux values
% *************************************************************
%

% D_flat_axis_Nradial=zeros(1,Nradial)+2;
% D_flat_axis_Nradial(2)=1;

psi2D_mask=psi2D.*mask_XZ_map;



%the first two psi values are not relevant
% xi_psi_Nradial(2)=xi_psi_Nradial(1);

for (n=2:Nradial)
    [max_val pos]=max(X_PR_map(:,n));
    Xmax_pos=X_PR_map(pos,n);
    [min_val pos]=min(X_PR_map(:,n));
    Xmin_pos=X_PR_map(pos,n);
    xi_psi_Nradial(n)=0.5*(Xmin_pos+Xmax_pos);
end


xi_pol_psi=polyfit(psi_scale,xi_psi_Nradial,4);
xi_psi_fit=polyval(xi_pol_psi,psi_scale);

xi_psi_Nradial=xi_psi_fit;
axis_pos=round(xi_psi_Nradial/DX)+mid_X;

