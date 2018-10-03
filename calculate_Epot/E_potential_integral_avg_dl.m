
% Integration calculation along psi star contour

        x_value=x_value_map(n,s);
        z_value=z_value_map(n,s);
        
        x_pos=round(x_value/DX_zoom+mid_X_zoom);
        z_pos=round(z_value/DZ_zoom+mid_Z_zoom);
        integ_value=integ_element_XZ_map(x_pos,z_pos);
%        integ_value=interp2(X_scale_zoom,Z_scale_zoom,integ_element_XZ_map,x_value,z_value);
        
        dl_value_X=dl_X_data(n,s);
        dl_value_Z=dl_Z_data(n,s);
        dl_value=sqrt(dl_value_X^2+dl_value_Z^2);
        
        inc_Epot=dl_value*integ_value;
        
        %index_value_map(n,s)=index_data;
        
        %dl_data(index_data)=dl_value;
        
        E_potential_X_data(index_data)=x_value;
        E_potential_Z_data(index_data)=z_value;
        E_potential_XZ_data(index_data)=Epot_value_prev+inc_Epot;
        Epot_value_prev=E_potential_XZ_data(index_data);

        integ_value_prev=integ_value;
        index_data=index_data+1;