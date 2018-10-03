Nbegin=1;
psi_value=-1;
Nbegin_rank=1;
% maximum space between two points to beconsidered to be from the same
% contour
DLMAX=4;

for (n=Ncontours_begin:Ncontours)
    if SAVE_LOG_FILE==1
        fprintf(fid, ' n = %d \n',n);
    end
    
    NL_RATIO=abs(round(NB_PHI/4)-abs(abs(phi_rank-round(NB_PHI/2))-round(NB_PHI/4)))/round(NB_PHI/4);
    if SAVE_LOG_FILE==1
        fprintf(fid, '--------------------------\n');
    end
    [psi_value Nlength_global nb_pieces contour_data_x contour_data_z] = describe_contour(Nbegin, Cont_psi_star, DLMAX);

%     psi_value=Cont_psi_star(1,Nbegin);
%     x_value_prev=Cont_psi_star(1,Nbegin);
%     z_value_prev=Cont_psi_star(2,Nbegin);
    
    Nlength=Nlength_global;
    if (nb_pieces==-1)
        disp('Excluding a contour because of broken pieces ... ')
    end
    
    if (nb_pieces~=-1) && (Nlength>MIN_LENGTH)

        if SAVE_LOG_FILE==1
            fprintf(fid, ' n_rank = %d \n',n_rank);
        end

        n_rank=n_rank+1;

        s_offset=0;
        
        Nlength_data(n_rank)=Nlength;
        psi_value_list(n_rank)=psi_value;

        
      
        NOFFSET=max(round(0.12*NL_RATIO*Nlength),MIN_LENGTH_HALF);
        NOFFSET=min(NOFFSET,round(0.22*Nlength)-1);
        [ Nlength_half_first Nlength_half Nlength_end ] = find_precise_halves_positions( Nlength, contour_data_x, contour_data_z, Line_ref_X, Line_ref_Z, NOFFSET ,MIN_LENGTH_HALF);
            
        if SAVE_LOG_FILE==1
            %fprintf(fid, ' Nlength_half_first = %d \n',Nlength_half_first);
            fprintf(fid, ' Nlength_half = %d \n',Nlength_half);
            %fprintf(fid, ' Nlength_end = %d \n',Nlength_end);
            fprintf(fid, ' Nlength = %d \n',Nlength);
        end
          
        Nstart_data(n_rank)=Nlength_half_first;
        Nhalf_data(n_rank)=Nlength_half;
        Nend_data(n_rank)=Nlength_end;
            
 
             
        x_value_map(n_rank,(1:Nlength))=contour_data_x;
        z_value_map(n_rank,(1:Nlength))=contour_data_z;
              
    end
    Nbegin=Nbegin+Nlength+1; 

end



% if (Nlength<=MIN_LENGTH)
%     n_rank=n_rank-1;
% end
    
Ncontours=n_rank;
