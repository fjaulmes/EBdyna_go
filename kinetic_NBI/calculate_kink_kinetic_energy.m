phi_rank_kinetic_energy=0;
EX_PR_map=squeeze(Efield_X_map_phi(phi_rank,:,:));
EZ_PR_map=squeeze(Efield_Z_map_phi(phi_rank,:,:));

    psi_star_dot_PR_map_phi=zeros(NB_THETA,NB_THETA,size_r);
    for index=1:NB_THETA
        psi_star_dot_PR_map_phi(index,:,:)=rotate_map_phi(psi_star_dot_omega_map,index)';
    end
    
    grad_Phi_tor_PR_map=squeeze(grad_Phi_tor_map_phi(phi_rank,:,:));
    psi_star_dot_PR_map=squeeze(psi_star_dot_PR_map_phi(phi_rank,:,:));
    
    Ephi_PR_map=psi_star_dot_PR_map./Rpos_PR_map(:,1:size_r)-grad_Phi_tor_PR_map;
    
% Ephi_PR_map=squeeze(Efield_phi_map_phi(phi_rank,:,:));
Etot_PR_map=sqrt(EX_PR_map.^2+EZ_PR_map.^2+Ephi_PR_map.^2);

Efield_phi_map_phi(phi_rank,:,:)=Etot_PR_map;

ksi_dot_PR_map=Etot_PR_map./B_PR_map;
ksi_dot_PR_map(isnan(ksi_dot_PR_map))=0;


ksi_dot_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,ksi_dot_PR_map(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
ksi_dot_XZ_map(isnan(ksi_dot_XZ_map))=0;
ksi_dot_XZ_map=ksi_dot_XZ_map';


for (x=3:sizeX-2)
%     for (z=3:sizeZ-2)
        phi_rank_kinetic_energy=phi_rank_kinetic_energy+...
            0.5*DPHI*DX*DX*mD*sum(ion_density_XZ_map(x,3:sizeZ-2).*Rpos_XZsmall_map(x,3:sizeZ-2).*ksi_dot_XZ_map(x,3:sizeZ-2).^2);
%     end
end

phi_rank_kinetic_energy