
% smooth integrand map

%integ_element_XZ_map=smooth_avg_small_map(integ_element_XZ_map);
%integ_element_XZ_map=smooth_avg_map(integ_element_XZ_map);


% integrate along contour

index_data=1;
s=1;



clear dl_X_data dl_Z_data;
clear E_potential_X_data E_potential_Z_data E_potential_XZ_data

for (n=1:Ncontours)
    Nlength=Nlength_data(n);
    dl_X_data(n,1)=0.5*(x_value_map(n,2)-x_value_map(n,1));
    dl_Z_data(n,1)=0.5*(z_value_map(n,2)-z_value_map(n,1)); 
    for(s=2:Nlength-1)
        dl_X_data(n,s)=0.5*(x_value_map(n,s+1)-x_value_map(n,s-1));
        dl_Z_data(n,s)=0.5*(z_value_map(n,s+1)-z_value_map(n,s-1));        
    end
    dl_X_data(n,Nlength)=0.5*(x_value_map(n,Nlength)-x_value_map(n,Nlength-1));
    dl_Z_data(n,Nlength)=0.5*(z_value_map(n,Nlength)-z_value_map(n,Nlength-1)); 
end

data_length=sum(Nlength_data);
E_potential_XZ_data=zeros(1,data_length);
E_potential_X_data=zeros(1,data_length);
E_potential_Z_data=zeros(1,data_length);

index_data=1;

for(n=1:Ncontours)
%for (n=24:38)
%     disp('-------------------------------');
%     disp(n);
    
    integ_value=0;
    integ_value_prev=0;
    inc_Epot_prev=0;
    Epot_value_prev=0;
    
    x_value_prev=x_value_map(n,1);
    z_value_prev=z_value_map(n,1);
    
    Nlength=Nlength_data(n);
    Nlength_half=Nhalf_data(n);


    for (s=1:Nlength_half)
        E_potential_integral_avg_dl;
    end
    if CATCH_BACK_POTENTIAL_INTEGRATION_TO_ZERO==1
        clear E_correction
        delta_Epot=0.5*Epot_value_prev;
        if SAVE_LOG_FILE==1
            fprintf(fid, '--------------------------\n');
            fprintf(fid, ' n_rank = %d \n',n);
            fprintf(fid, ' delta_Epot = %d \n',delta_Epot);
        end
        E_potential_XZ_data(index_data-Nlength_half:index_data-1)=E_potential_XZ_data(index_data-Nlength_half:index_data-1)-delta_Epot;
        E_correction=-(0:Nlength_half-1)*(2*delta_Epot/(Nlength_half))+delta_Epot;
        E_potential_XZ_data(index_data-Nlength_half:index_data-1)=E_potential_XZ_data(index_data-Nlength_half:index_data-1)+E_correction;
    end
    
    Epot_value_prev=0;
    
    for (s=Nlength_half+1:Nlength)
        E_potential_integral_avg_dl;
    end
    if CATCH_BACK_POTENTIAL_INTEGRATION_TO_ZERO==1
        clear E_correction
        delta_Epot=0.5*Epot_value_prev;
        Nlength_half=Nlength-Nlength_half+1;
        E_potential_XZ_data(index_data-Nlength_half:index_data-1)=E_potential_XZ_data(index_data-Nlength_half:index_data-1)-delta_Epot;
        E_correction=-(0:Nlength_half-1)*(2*delta_Epot/(Nlength_half))+delta_Epot;
        E_potential_XZ_data(index_data-Nlength_half:index_data-1)=E_potential_XZ_data(index_data-Nlength_half:index_data-1)+E_correction;
    end
end

if SAVE_LOG_FILE==1
    fprintf(fid, '--------------------------\n');
end


E_potential_RZ_map=griddata(E_potential_X_data,E_potential_Z_data,E_potential_XZ_data,XX,ZZ,'cubic');
E_potential_RZ_map=E_potential_RZ_map';
E_potential_RZ_map(isnan(E_potential_RZ_map))=0;

% E_potential_XZ_zoom_map=griddata(E_potential_X_data,E_potential_Z_data,E_potential_XZ_data,XX_zoom,ZZ_zoom,'cubic');
% E_potential_XZ_zoom_map=E_potential_XZ_zoom_map';
% E_potential_XZ_zoom_map(isnan(E_potential_XZ_zoom_map))=0;


%E_potential_XZ_zoom_map=conv2(E_potential_XZ_zoom_map,[1 1 1 ; 1 1 1 ;1 1 1]);
