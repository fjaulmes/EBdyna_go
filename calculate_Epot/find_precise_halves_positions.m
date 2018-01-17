% determines values of position of reference axis intersections on curved contour

function [ Nlength_half_first Nlength_half Nlength_end ] = find_precise_halves_positions( Nlength, X_coords, Z_coords, Line_ref_X, Line_ref_Z, NOFFSET, MIN_LENGTH_HALF )

    Nlength_half_first=-1;
    Nlength_half=-1;
    Nlength_end=-1;
    Nmiddle=round(0.5*Nlength);
    Nbegin=1;


    Nlength_half_first=eval_intersection_ref_axis_pot(Nbegin,Nbegin+Nmiddle-NOFFSET-1,X_coords,Z_coords,Line_ref_X,Line_ref_Z);

    if (Nlength_half_first==1)||(Nlength_half_first==Nmiddle-NOFFSET)
        Nlength_half_first=eval_intersection_ref_axis_pot(Nbegin+Nmiddle-NOFFSET-1,Nbegin+Nmiddle-MIN_LENGTH_HALF-1,X_coords,Z_coords,Line_ref_X,Line_ref_Z);
        Nlength_half_first=Nlength_half_first+Nmiddle-NOFFSET-1;
    end

    if (Nlength_half_first==Nmiddle-NOFFSET-1)||(Nlength_half_first==Nmiddle-MIN_LENGTH_HALF)
        Nlength_half=eval_intersection_ref_axis_pot(Nbegin+NOFFSET,Nbegin+Nlength-NOFFSET,X_coords,Z_coords,Line_ref_X,Line_ref_Z);
        Nlength_half=Nlength_half+NOFFSET;
        Nlength_half_first=-1;
    end

    if (Nlength_half_first==-1)
        Nlength_end=eval_intersection_ref_axis_pot(Nbegin+Nmiddle+NOFFSET,Nbegin+Nlength-1,X_coords,Z_coords,Line_ref_X,Line_ref_Z);
        if (Nlength_end>1)
            Nlength_end=Nlength_end+Nmiddle+NOFFSET;
        else
            Nlength_end=-1;
        end
    end



    if (Nlength_half==-1)||(Nlength_half+Nlength_end>Nlength)
        if (Nlength_end>Nlength_half_first)
            if (Nlength_half_first>3)
                Nlength_half=eval_intersection_ref_axis_pot(Nbegin+Nlength_half_first+2,Nbegin+Nlength_end-MIN_LENGTH_HALF,X_coords,Z_coords,Line_ref_X,Line_ref_Z);
                Nlength_half=Nlength_half+Nlength_half_first+2;
            else
                Nlength_half=eval_intersection_ref_axis_pot(Nbegin+NOFFSET,Nbegin+Nlength_end-MIN_LENGTH_HALF,X_coords,Z_coords,Line_ref_X,Line_ref_Z);
                Nlength_half=Nlength_half+NOFFSET;
            end
        else

            if (Nlength-MIN_LENGTH_HALF)>Nlength_half_first+NOFFSET
                Nlength_half=eval_intersection_ref_axis_pot(Nbegin+Nlength_half_first+NOFFSET,Nbegin+Nlength-MIN_LENGTH_HALF,X_coords,Z_coords,Line_ref_X,Line_ref_Z);
                Nlength_half=Nlength_half+Nlength_half_first+NOFFSET;
            else
                Nlength_half=eval_intersection_ref_axis_pot(Nbegin+Nlength_half_first+2,Nbegin+Nlength-MIN_LENGTH_HALF,X_coords,Z_coords,Line_ref_X,Line_ref_Z);
                Nlength_half=Nlength_half+Nlength_half_first+2;
            end

        end
    end


    if (Nlength_end==Nlength_half)
        Nlength_end=eval_intersection_ref_axis_pot(Nbegin+Nlength_half+1,Nbegin+Nlength-1,X_coords,Z_coords,Line_ref_X,Line_ref_Z);
        if (Nlength_end>1)
            Nlength_end=Nlength_end+Nlength_half+1;
        else
            Nlength_end=-1;
        end
    end


    if  ((Nlength_half_first+Nlength_half)>Nlength)
        Nlength_end=Nlength_half;
        Nlength_half=Nlength_half_first;
        Nlength_half_first=-1;
    end



    if (Nlength_end==-1)&&(Nlength_half_first==-1)
        Nlength_half_first=eval_intersection_ref_axis_pot(Nbegin,Nbegin+NOFFSET,X_coords,Z_coords,Line_ref_X,Line_ref_Z);
    end

    if (Nlength_end~=-1)
        Nlength_half_first=-1;
    end

end

