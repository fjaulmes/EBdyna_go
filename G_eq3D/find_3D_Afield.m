function [field] = find_3D_Afield(x,request )
global maps dim par
%find_RMP_field Summary of this function goes here
%   Detailed explanation goes here

X_ind=((x(:,1,:)-dim.R0)*dim.DX_inv)+dim.mid_Xzero;
Z_ind=( x(:,2,:)        *dim.DZ_inv)+dim.mid_Z;
phi_ind=mod(x(:,3,:),2*pi/dim.n3D.symm)  *dim.Dphi_inv+1;

for i=1:length(request)
switch par.coord_syst
    case 'flux'
        psi_norm=ba_interp2(maps(1).psi_norm_XZ,Z_ind,X_ind,'linear');
        psi_ind=dim.n3D.ind_2D_to_3D(psi_norm);
        if par.interp_scheme==1
            [X_ind,Z_ind] = interp_index_list_2D ([dim.size_X,dim.size_Z],X_ind,Z_ind);     % Return the indexes / slopes as reps. X_ind and Z_ind
        end
        theta = interpolate_theta_XZ(X_ind,Z_ind);
        theta_ind=theta*dim.Dtheta_inv+1;

        field.(request{i})=ba_interp3(maps(1).n3D.(request{i}),psi_ind,theta_ind,phi_ind,'linear');
    case 'toroidal'
        field.(request{i})=ba_interp3(maps(1).n3D.(request{i}),Z_ind,X_ind,phi_ind,'linear');
end
end


end