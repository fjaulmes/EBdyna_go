
gradP_X=zeros(sizeX,sizeZ);
gradP_Z=zeros(sizeX,sizeZ);

% smooth_B_map=smooth_small_map(Btot_XZ_map);



for (x=3:sizeX-2)
    for (z=3:sizeZ-2)
        gradP_X(x,z)=(1/12)*(-NBI_Phot_XZ_map_small(x+2,z)+NBI_Phot_XZ_map_small(x-2,z))+(2/3)*(NBI_Phot_XZ_map_small(x+1,z)-NBI_Phot_XZ_map_small(x-1,z));
        gradP_Z(x,z)=(1/12)*(-NBI_Phot_XZ_map_small(x,z+2)+NBI_Phot_XZ_map_small(x,z-2))+(2/3)*(NBI_Phot_XZ_map_small(x,z+1)-NBI_Phot_XZ_map_small(x,z-1));
    end
end

gradP_X=gradP_X/DX;
gradP_Z=gradP_Z/DX;



