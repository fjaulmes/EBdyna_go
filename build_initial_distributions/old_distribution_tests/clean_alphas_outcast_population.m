        alphas_psi=interp2(scale_X,scale_Z,radial_XZsmall_map',alphas_pos_x,alphas_pos_z,'nearest');
        % correct the particles that are out of simulation domain
        % by giving them a position randomly in the initial distribution
        %outcast=[find(alphas_psi>=scale_X(end-2));find(alphas_pos_z>=scale_Z(end-2));find(alphas_pos_x<=scale_X(3));find(alphas_pos_z<=scale_Z(3))];
        outcast=find(alphas_psi>simulation_size_r+1);
        recast=ceil(rand(size(outcast))*Nalphas_simulated-1)+1;
        % we only give them a new position and 
        % keep their initial energy
        alphas_pos_x(outcast)=X0(recast);
        alphas_pos_z(outcast)=Z0(recast);
        alphas_pos_phi(outcast)=rand(size(recast))*2*pi;
        alphas_mm(outcast)=(alphas_Ekin(outcast)-0.5*(mHe*alphas_vpll(outcast).^2)/eV)./interp2(scale_X,scale_Z,Btot_XZ_map',alphas_pos_x(outcast),alphas_pos_z(outcast),'*linear');

        alphas_psi=interp2(scale_X,scale_Z,radial_XZsmall_map',alphas_pos_x,alphas_pos_z,'nearest');
