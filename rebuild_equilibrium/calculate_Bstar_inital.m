

BstarX_XZ_map_ini=zeros(NZ,NZ);
BstarZ_XZ_map_ini=zeros(NZ,NZ);
BstarX_PR_map=zeros(NP,Nradial);
BstarZ_PR_map=zeros(NP,Nradial);

if NOQ1SURF==0
	psi2D=psi_star_XZ_map;


	gpsi_R=zeros(NZ,NZ);
	gpsi_Z=zeros(NZ,NZ);

	for (x=3:NZ-2)
		for (z=3:NZ-2)
			gpsi_R(x,z)=(1/12)*(-psi2D(x+2,z)+psi2D(x-2,z))+(2/3)*(psi2D(x+1,z)-psi2D(x-1,z));
			gpsi_Z(x,z)=(1/12)*(-psi2D(x,z+2)+psi2D(x,z-2))+(2/3)*(psi2D(x,z+1)-psi2D(x,z-1));
	%         gpsi_Z(x,z)=0.5*(psi2D(x,z+1)-psi2D(x,z-1));
		end
	end
	gpsi_R=gpsi_R/DX;
	gpsi_Z=gpsi_Z/DX;

	BstarX_XZ_map_ini=-gpsi_Z./Rpos_map;
	BstarZ_XZ_map_ini=gpsi_R./Rpos_map;

	BstarX_XZ_map_ini=BstarX_XZ_map_ini.*mask_XZ_map;
	BstarZ_XZ_map_ini=BstarZ_XZ_map_ini.*mask_XZ_map;


	%%
	BH_data=reshape(BstarX_XZ_map_ini(:,:)',NZ*NZ,1);
	BstarX_PR_map=griddata(X_scale_data,Z_scale_data,BH_data,Rpos_PR_map-R0,Z_PR_map,'cubic');

	BH_data=reshape(BstarZ_XZ_map_ini(:,:)',NZ*NZ,1);
	BstarZ_PR_map=griddata(X_scale_data,Z_scale_data,BH_data,Rpos_PR_map-R0,Z_PR_map,'cubic');

end

Bstar_PR_map=sqrt(BstarX_PR_map.^2+BstarZ_PR_map.^2);

        


