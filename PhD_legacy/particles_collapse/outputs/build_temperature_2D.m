function [ posX_scale posZ_scale n_2D_map temp_2D_map ] = build_temperature_2D(posX,posZ,Ekin_dist,XINF,XSUP,ZINF,ZSUP,DELTA_2D_RES)
% creates a temperature map according to positions and energy of particles

nb_X=round((XSUP-XINF)/DELTA_2D_RES);
nb_Z=round((ZSUP-ZINF)/DELTA_2D_RES);

n_2D_map=zeros(nb_X,nb_Z);
temp_2D_map=zeros(nb_X,nb_Z);
posX_scale=linspace(XINF,XSUP,nb_X)+0.5*DELTA_2D_RES;
posZ_scale=linspace(ZINF,ZSUP,nb_Z)+0.5*DELTA_2D_RES;



for (x=1:nb_X)
    for (z=1:nb_Z)
        POP_2D_CELL=find((posX>=(posX_scale(x)-0.5*DELTA_2D_RES)).*(posZ>=(posZ_scale(z)-0.5*DELTA_2D_RES)).* ...
            (posX<(posX_scale(x)+0.5*DELTA_2D_RES)).*(posZ<(posZ_scale(z)+0.5*DELTA_2D_RES)));
        n_cell=length(POP_2D_CELL);
        Ekin_cell=sum(Ekin_dist(POP_2D_CELL));
        n_2D_map(x,z)=n_cell;
        temp_2D_map(x,z)=Ekin_cell;
    end
end

end

