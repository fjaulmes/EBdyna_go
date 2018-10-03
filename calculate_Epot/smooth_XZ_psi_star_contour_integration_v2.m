
% smooth integrand map
%    integ_element_XZ_map=smooth_avg_map(integ_element_XZ_map);

% integrate along contour

index_data=1;


clear x_value_data z_value_data dl_X_data dl_Z_data integ_value_data s_half_data s_length_data
clear E_potential_X_data E_potential_Z_data E_potential_XZ_data

for (n=1:Ncontours)
    Nlength_half=Nhalf_data(n);
    Nlength=Nlength_data(n);
    
    index_data=1;
    x_value_data(n,1)=x_value_map(n,1);
    z_value_data(n,1)=z_value_map(n,1);
    xpos=round(x_value_map(n,1)/DX_zoom+mid_X_zoom);
    zpos=round(z_value_map(n,1)/DX_zoom+mid_Z_zoom);
    dl_X_data(n,1)=0.5*(x_value_map(n,2)-x_value_map(n,Nlength));
    dl_Z_data(n,1)=0.5*(z_value_map(n,2)-z_value_map(n,Nlength)); 
    integ_value_data(n,1)=integ_element_XZ_map(xpos,zpos);
 
    for(s=3:2:Nlength_half-2)
        index_data=index_data+1;
        x_value_data(n,index_data)=x_value_map(n,s);
        z_value_data(n,index_data)=z_value_map(n,s);
        xpos=round(x_value_map(n,s)/DX_zoom+mid_X_zoom);
        zpos=round(z_value_map(n,s)/DX_zoom+mid_Z_zoom);
        integ_value_data(n,index_data)=integ_element_XZ_map(xpos,zpos);
        dl_X_data(n,index_data)=(x_value_map(n,s+1)-x_value_map(n,s-1));
        dl_Z_data(n,index_data)=(z_value_map(n,s+1)-z_value_map(n,s-1));        
    end
    
    if (s==Nlength_half-2)
        s_prev=max(Nlength_half-1,1);
    else
        s_prev=max(Nlength_half-2,1);
    end
    
    s=Nlength_half;
    index_data=index_data+1;
    x_value_data(n,index_data)=x_value_map(n,s);
    z_value_data(n,index_data)=z_value_map(n,s);    
    xpos=round(x_value_map(n,s)/DX_zoom+mid_X_zoom);
    zpos=round(z_value_map(n,s)/DX_zoom+mid_Z_zoom);
    integ_value_data(n,index_data)=integ_element_XZ_map(xpos,zpos);
    dl_X_data(n,index_data)=x_value_map(n,s+1)-x_value_map(n,s_prev);
    dl_Z_data(n,index_data)=z_value_map(n,s+1)-z_value_map(n,s_prev);
    s_half_data(n)=index_data;
    
    for(s=Nlength_half+2:2:Nlength-2)
        index_data=index_data+1;
        x_value_data(n,index_data)=x_value_map(n,s);
        z_value_data(n,index_data)=z_value_map(n,s);
        xpos=round(x_value_map(n,s)/DX_zoom+mid_X_zoom);
        zpos=round(z_value_map(n,s)/DX_zoom+mid_Z_zoom);
        integ_value_data(n,index_data)=integ_element_XZ_map(xpos,zpos);
        dl_X_data(n,index_data)=(x_value_map(n,s+1)-x_value_map(n,s-1));
        dl_Z_data(n,index_data)=(z_value_map(n,s+1)-z_value_map(n,s-1));        
    end

    if (s==Nlength_half-2)
        s_prev=Nlength-1;
    else
        s_prev=Nlength-2;
    end
    
    s=Nlength;
    index_data=index_data+1;
    x_value_data(n,index_data)=x_value_map(n,s);
    z_value_data(n,index_data)=z_value_map(n,s);    
    xpos=round(x_value_map(n,s)/DX_zoom+mid_X_zoom);
    zpos=round(z_value_map(n,s)/DX_zoom+mid_Z_zoom);
    integ_value_data(n,index_data)=integ_element_XZ_map(xpos,zpos);
    dl_X_data(n,index_data)=(x_value_map(n,s)-x_value_map(n,s_prev));
    dl_Z_data(n,index_data)=(z_value_map(n,s)-z_value_map(n,s_prev));
    s_length_data(n)=index_data;
    
end

data_length=sum(s_length_data);
E_potential_XZ_data=zeros(1,data_length);
E_potential_X_data=zeros(1,data_length);
E_potential_Z_data=zeros(1,data_length);

index_data=1;

for (n=1:Ncontours)
%     disp('-------------------------------');
%     disp(n);
    
    integ_value=0;
    integ_value_prev=0;
    inc_Epot_prev=0;
    Epot_value_prev=0;
    
    Nlength=s_length_data(n);
    Nlength_half=s_half_data(n);


    
    for (s=1:Nlength_half)
        E_potential_integral_avg_dl_v2;
    end
    
    if CATCH_BACK_POTENTIAL_INTEGRATION_TO_ZERO==1
        delta_Epot=0.5*Epot_value_prev;
        clear E_correction
        E_potential_XZ_data(index_data-Nlength_half:index_data-1)=E_potential_XZ_data(index_data-Nlength_half:index_data-1)-delta_Epot;
        E_correction=-(0:Nlength_half-1)*(2*delta_Epot/(Nlength_half))+delta_Epot;
        E_potential_XZ_data(index_data-Nlength_half:index_data-1)=E_potential_XZ_data(index_data-Nlength_half:index_data-1)+E_correction;
    end
    

    
    Epot_value_prev=0;
    
    for (s=Nlength_half+1:Nlength)
        
        E_potential_integral_avg_dl_v2;
        
    end

    if CATCH_BACK_POTENTIAL_INTEGRATION_TO_ZERO==1
        delta_Epot=0.5*Epot_value_prev;
        Nlength_half=Nlength-Nlength_half+1;
        clear E_correction
        E_potential_XZ_data(index_data-Nlength_half:index_data-1)=E_potential_XZ_data(index_data-Nlength_half:index_data-1)-delta_Epot;
        E_correction=-(0:Nlength_half-1)*(2*delta_Epot/(Nlength_half))+delta_Epot;
        E_potential_XZ_data(index_data-Nlength_half:index_data-1)=E_potential_XZ_data(index_data-Nlength_half:index_data-1)+E_correction;
    end
    %end
    psi_value_prev=psi_value;
end





E_potential_RZ_map=griddata(E_potential_X_data,E_potential_Z_data,E_potential_XZ_data,XX,ZZ,'cubic');
E_potential_RZ_map=E_potential_RZ_map';
E_potential_RZ_map(isnan(E_potential_RZ_map))=0;

% E_potential_XZ_zoom_map=griddata(E_potential_X_data,E_potential_Z_data,E_potential_XZ_data,XX_zoom,ZZ_zoom,'cubic');
% E_potential_XZ_zoom_map=E_potential_XZ_zoom_map';
% E_potential_XZ_zoom_map(isnan(E_potential_XZ_zoom_map))=0;


%E_potential_XZ_zoom_map=conv2(E_potential_XZ_zoom_map,[1 1 1 ; 1 1 1 ;1 1 1]);
