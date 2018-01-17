
E_potential_RZ_map=zeros(NZ,NZ);
E_potential_PR_map=zeros(NP,Nradial);

clear E_potential_X_data E_potential_Z_data E_potential_XZ_data
clear x_value_map z_value_map integ_value_map dl_value_map index_value_map
clear Nlength_data Nhalf_data Nend_data Nstart_data
clear Cont_psi_star hc 
clear final_psi_values Psih_final_recalc


% parameters for the contour extraction

SMALL_PSI_SEP_RANK=8;
MIN_LENGTH=280;
MIN_LENGTH_HALF=round((0.5*MIN_LENGTH)/2);


% extract contour from psi star XZ map


Npsi_values=round(xPsih_zero+1);

% Psih_final_recalc=max(Psih_final,-0.00001);
Psih_final_recalc=max(psi_star_final./(1+0.01*(0:Nradial-1)/(Nradial-1)),0);

defining_psi_value_final_values;

if (psi_sep_rank>2)
    psi_sep_rank=psi_sep_rank+1;
end

if (psi_sep_rank>1)
    final_psi_values(psi_sep_rank)=0.01*final_psi_values(psi_sep_rank-1)+0.99*psi_limit13;
end

% DELTA_PSI_SEP=1e-6;
% We need to get closer to the separatrix psi value
if EXPULSION_FRAME_NUMBER>1
    if (f>=EXPULSION_FRAME_NUMBER-1) && (f<=EXPULSION_FRAME_NUMBER+1)
        SEP_PSI_COEF=0.7;
    else
        SEP_PSI_COEF=0.85;
    end
else
    SEP_PSI_COEF=0.8;
end
if (psi_sep_rank>1)
    if (psi_sep_rank>=3)
        if (psi_sep_rank>=4)
            final_psi_values(psi_sep_rank-2)=0.5*final_psi_values(psi_sep_rank-3)+0.5*final_psi_values(psi_sep_rank-2);
        end
        final_psi_values(psi_sep_rank-1)=0.5*final_psi_values(psi_sep_rank-2)+0.5*final_psi_values(psi_sep_rank-1);
    end
    final_psi_values(psi_sep_rank)=(1-SEP_PSI_COEF)*final_psi_values(psi_sep_rank-1)+SEP_PSI_COEF*final_psi_values(psi_sep_rank);
end

% adjust psi_star_XZ_zoom_map so psi star values outside
% are not taken into account for the contours description
psi_star_XZ_zoom_map=psi_star_XZ_zoom_map_copy;

if (psi_sep_rank>=2)
    for (x=1:NZ_zoom)
        for (z=1:NZ_zoom)
            if (psi_star_XZ_zoom_map(x,z)<psi_limit13)
                psi_star_XZ_zoom_map(x,z) = -2*psi_star_final(1);
            end
        end
    end
end

[val row_ind] =max(psi_star_XZ_zoom_map, [], 1);
[val col_ind] =max(val); 
% The maximum value is at [row_ind(col_ind), col_ind]
Xmax_psi=row_ind(col_ind);
Zmax_psi=col_ind;


[val Xmax_psi]=max(max(psi_star_XZ_zoom_map'));
[val Zmax_psi]=max(psi_star_XZ_zoom_map(Xmax_psi,:));

if (psi_sep_rank>=3)
    figure(6);
    if psi_sep_rank==3
        val = final_psi_values(psi_sep_rank:psi_sep_rank);
        Cont_psi_star=contour(X_scale_zoom,Z_scale_zoom,psi_star_XZ_zoom_map',[val val]);
    else
        Cont_psi_star=contour(X_scale_zoom,Z_scale_zoom,psi_star_XZ_zoom_map',final_psi_values(1:psi_sep_rank));
    end
    Ncontours_region3 = find_number_of_contours(Cont_psi_star);
    if DISPLAY_OUTPUTS==1
        pause(0.1);
    else
        hold on;
        %Cont_psi_star_object=contourcs(X_scale_zoom,Z_scale_zoom,psi_star_XZ_zoom_map',final_psi_values(1:psi_sep_rank+1));
    end

else
    Ncontours_region3=0;
end
%hold on


index_data=1;

n_rank=0;
psi_value_prev=0;
Nbegin=1;
N_separatrix_13=0;
matching_index=0;
integ_value_prev=0;
missing_psi_value_rank=0;

Nlength=MIN_LENGTH+1;


Ncontours=Ncontours_region3;
Ncontours_begin=1;
n_rank=0;
Nlength_data=zeros(1,1);

if (psi_sep_rank>1)
    following_contours_to_get_halves;
%     n_rank=n_rank-1;
    Ncontours_region3=n_rank;
end

N_separatrix_13=n_rank;


% disp('pause')
% pause



% Then getting the contours in region 1
% Leaving again separatrix outside
Npsi_values=round(xPsih_zero+1);
Psih_final_recalc=max(psi_star_final./(1+0.15*(0:Nradial-1)/(Nradial-1)),0);
defining_psi_value_final_values;
% psi_sep_rank=psi_sep_rank+1;

if (f>2)
%     delta_rx=max(1.1+0.1*f,1.2);
%     delta_rx=min(delta_rx,2.6);
    delta_rx=0.1;
    define_RZ_mask_rx_rank;

    psi_star_XZ_zoom_map=psi_star_XZ_zoom_map_copy+3*(1-mask_XZ_zoom_map_reconnection)*abs(min_Psih);
    [val row_ind] =min(psi_star_XZ_zoom_map, [], 1);
    [val col_ind] =min(val);
    % The min value is at [row_ind(col_ind), col_ind]
    Xmin_psi=row_ind(col_ind);
    Zmin_psi=col_ind;
    
%     [val col_ind]=max(psi_star_XZ_zoom_map(Xmin_psi,:));
%     
%     dist_max_core=sqrt((Xmax_psi-Xmin_psi)^2+(Zmax_psi-Zmin_psi)^2);
%     dist=abs(col_ind-Zmin_psi);
%     dist_max_core=max(dist_max_core,dist);
%     
%     % removing points too far away from the core
%     for (x=2:NZ_zoom-1)
%         for (z=2:NZ_zoom-1)
%             dist=sqrt((x-Xmin_psi)^2+(z-Zmin_psi)^2);
%             if(dist>dist_max_core)
%                 mask_XZ_zoom_map_reconnection(x,z)=0;
%             end
%         end
%     end
    %psi_star_XZ_zoom_map=min(psi_star_XZ_zoom_map,psi_limit13);
else
    delta_rx=0;
    define_RZ_mask_rx_rank;
end

psi_star_XZ_zoom_map=psi_star_XZ_zoom_map_copy+3*(1-mask_XZ_zoom_map_reconnection)*abs(min_Psih);

if psi_limit13>0
    init_psi_values_region1=psi_sep_rank;
    if (f>1)
        for (x=1:NZ_zoom)
            for (z=1:NZ_zoom)
                if (mask_XZ_zoom_map_reconnection(x,z)==1)
                    if (psi_star_XZ_zoom_map(x,z)>psi_limit13)
                        psi_star_XZ_zoom_map(x,z) = 3*abs(min_Psih);
                    end
                end
            end
        end
    end
end


clear hc
clear Cont_psi_star

% Trying to get as close to separatrix as possible
Ncontours_region1=0;

if (psi_limit13>0)&&(f<EXPULSION_FRAME_NUMBER+1)
    
    if DISPLAY_OUTPUTS==1
        figure(6);
        [Cont_psi_star hc]=contour(X_scale_zoom,Z_scale_zoom,psi_star_XZ_zoom_map',final_psi_values(init_psi_values_region1:end));
        pause(0.1);
        Ncontours_region1 = size(get(hc,'Children'),1);
    else
        figure(6);
        Cont_psi_star=contour(X_scale_zoom,Z_scale_zoom,psi_star_XZ_zoom_map',final_psi_values(init_psi_values_region1:end));
        Ncontours_region1 = find_number_of_contours(Cont_psi_star);
        pause(0.1);
    end

end

% Total number of contours region 3 + region 1 
Ncontours=Ncontours_region3+Ncontours_region1;
disp(strcat('Ncontours = ',num2str(Ncontours)));

if psi_limit13>0
    Ncontours_begin=Ncontours_region3+1;
    following_contours_to_get_halves;
    if SAVE_LOG_FILE==1
        fprintf(fid, '--------------------------\n');
    end
    Ncontours=n_rank-1;
    Ncontours_region1=Ncontours-Ncontours_region3;
end


% adding specific contours for separatrix
% add a specific contour for the separatrix at r=rx
% only if the interpolation allows for sharp gradients
% and large enough frame number

% if (strcmp(DT_INTERPOLATION_METHOD,'quadratic'))&&((f>round(0.2*SMALL_FRAME_NUMBER))||(f<2))
if (strcmp(DT_INTERPOLATION_METHOD,'quadratic'))
    
    Ncontours=Ncontours+1;
    n_rank=Ncontours;
    NB_DL=4000;
    DL_DOMEGA=(2*pi)/(NB_DL);
    Nlength=NB_DL;
    
    omega=0;
    for(s=1:NB_DL)
        omega=(NP-1)*((s-1)*DL_DOMEGA)/(2*pi)+1;
        x_value_map(n_rank,s)=interp2((1:NP),(1:Nradial),X_PR_map',omega,rx_precise,'*linear');
        z_value_map(n_rank,s)=interp2((1:NP),(1:Nradial),Z_PR_map',omega,rx_precise,'*linear');
    end
    
    [ Nlength_half_first Nlength_half Nlength_end ] = find_precise_halves_positions( Nlength, x_value_map(n_rank,:), z_value_map(n_rank,:), Line_ref_X, Line_ref_Z, MIN_LENGTH_HALF, MIN_LENGTH_HALF );
    Nstart_data(n_rank)=Nlength_half_first;
    Nhalf_data(n_rank)=Nlength_half;
    Nlength_data(n_rank)=NB_DL;
    Nend_data(n_rank)=Nlength_end;
    
    final_psi_values(end+1)=interp1(1:Nradial,psi_star_initial,rx_precise);
    
end



% shifting data so that first data point is on axis
% and that next x axis intersection point is at "Nhalf" position in array

for (n=1:Ncontours)
    Nlength=Nlength_data(n);
    Nlength_half_first=Nstart_data(n);
    Nlength_half=Nhalf_data(n);
    Nlength_end=Nend_data(n);
    DNlength_end=Nlength-Nlength_end;
    if (Nlength_half_first~=-1)
        clear buffer
        
        buffer=x_value_map(n,1:Nlength_half_first);
        x_value_map(n,1:Nlength-Nlength_half_first)=x_value_map(n,Nlength_half_first+1:Nlength);
        x_value_map(n,Nlength-Nlength_half_first+1:Nlength)=buffer;
            
        clear buffer
        buffer=z_value_map(n,1:Nlength_half_first);
        z_value_map(n,1:Nlength-Nlength_half_first)=z_value_map(n,Nlength_half_first+1:Nlength);
        z_value_map(n,Nlength-Nlength_half_first+1:Nlength)=buffer;

        
        Nhalf_data(n)=Nlength_half-Nlength_half_first;
    elseif (Nlength_end~=-1)&&(DNlength_end>0)
        clear buffer
        buffer=x_value_map(n,Nlength_end+1:Nlength);
        x_value_map(n,DNlength_end+1:Nlength)=x_value_map(n,1:Nlength_end);
        x_value_map(n,1:DNlength_end)=buffer;   
                    
        clear buffer
        buffer=z_value_map(n,Nlength_end+1:Nlength);
        z_value_map(n,DNlength_end+1:Nlength)=z_value_map(n,1:Nlength_end);
        z_value_map(n,1:DNlength_end)=buffer;

        Nhalf_data(n)=Nlength_half+DNlength_end;

    end
    
end


% inserting extra precise values to be exactly on reference axis
Line_ref_zoom_X=interp1((1:(Nradial-1)*2+1),Line_ref_X,(0:0.5:(Nradial-1)*2)+1,'cubic');
Line_ref_zoom_Z=interp1((1:(Nradial-1)*2+1),Line_ref_Z,(0:0.5:(Nradial-1)*2)+1,'cubic');

HORIZONTAL_LIKE_REF=0;
if (phi<(pi/4))||((phi>3*pi/4)&&(phi<5*pi/4))||(phi>(7*pi/4))
    HORIZONTAL_LIKE_REF=1;
end

LENGTH_DATA_MAX=size(x_value_map,2);

for (n_rank=1:Ncontours)
    
    
    Nlength=Nlength_data(n_rank);
    Nlength_half=Nhalf_data(n_rank);
    %dist_values=sqrt((line_ref_pot_X-x_value_map(n_rank,1)).^2+(line_ref_pot_Z-z_value_map(n_rank,1)).^2);

    rank_axis_begin=eval_coords_ref_axis_pot(x_value_map(n_rank,1),z_value_map(n_rank,1),Line_ref_zoom_X,Line_ref_zoom_Z);
    rank_axis_end=eval_coords_ref_axis_pot(x_value_map(n_rank,Nlength),z_value_map(n_rank,Nlength),Line_ref_zoom_X,Line_ref_zoom_Z);
    rank_axis_middle=eval_coords_ref_axis_pot(x_value_map(n_rank,Nlength_half),z_value_map(n_rank,Nlength_half),Line_ref_zoom_X,Line_ref_zoom_Z);
        
    %replacing values of the contour
    x_value_map(n_rank,1)=Line_ref_zoom_X(rank_axis_begin);
    z_value_map(n_rank,1)=Line_ref_zoom_Z(rank_axis_begin);
    
    x_value_map(n_rank,Nlength)=Line_ref_zoom_X(rank_axis_end);
    z_value_map(n_rank,Nlength)=Line_ref_zoom_Z(rank_axis_end);
    
    x_value_map(n_rank,Nlength_half)=Line_ref_zoom_X(rank_axis_middle);
    z_value_map(n_rank,Nlength_half)=Line_ref_zoom_Z(rank_axis_middle);
 
    
    
    if (HORIZONTAL_LIKE_REF==1)
        x_value_map(n_rank,1)=0.5*(x_value_map(n_rank,2)+x_value_map(n_rank,Nlength-1));
        %z_value_map(n_rank,1)=0.5*(z_value_map(n_rank,2)+z_value_map(n_rank,Nlength+1));
        
        x_value_map(n_rank,Nlength)=0.5*(x_value_map(n_rank,2)+x_value_map(n_rank,Nlength-1));
        %z_value_map(n_rank,Nlength+2)=0.5*(z_value_map(n_rank,2)+z_value_map(n_rank,Nlength+1));
        
        x_value_map(n_rank,Nlength_half)=0.5*(x_value_map(n_rank,Nlength_half-1)+x_value_map(n_rank,Nlength_half+1));

    else
        %x_value_map(n_rank,1)=0.5*(x_value_map(n_rank,2)+x_value_map(n_rank,Nlength+1));
        z_value_map(n_rank,1)=0.5*(z_value_map(n_rank,2)+z_value_map(n_rank,Nlength-1));
        
        %x_value_map(n_rank,Nlength+2)=0.5*(x_value_map(n_rank,2)+x_value_map(n_rank,Nlength+1));
        z_value_map(n_rank,Nlength)=0.5*(z_value_map(n_rank,2)+z_value_map(n_rank,Nlength-1));
        
        z_value_map(n_rank,Nlength_half)=0.5*(z_value_map(n_rank,Nlength_half-1)+z_value_map(n_rank,Nlength_half+1));

    end
    
    if (Nlength<LENGTH_DATA_MAX)
        x_value_map(n_rank,Nlength+1:LENGTH_DATA_MAX)=x_value_map(n_rank,Nlength+1:LENGTH_DATA_MAX).*0+x_value_map(n_rank,Nlength);
        z_value_map(n_rank,Nlength+1:LENGTH_DATA_MAX)=z_value_map(n_rank,Nlength+1:LENGTH_DATA_MAX).*0+z_value_map(n_rank,Nlength);
    end
%    Nlength_data(n_rank)=Nlength;
%     if (Nstart_data(n_rank)~=-1)
%         Nstart_data(n_rank)=Nstart_data(n_rank);
%     end
%     if (Nhalf_data(n_rank)~=-1)
%         Nhalf_data(n_rank)=Nhalf_data(n_rank);
%     end
%     if (Nend_data(n_rank)==Nlength)
%         Nend_data(n_rank)=Nlength_data(n_rank);
%     end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% disp('pause')
% pause


