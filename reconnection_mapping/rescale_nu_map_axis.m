



if(omega_cont13==1)  
    nu_min_cont13=interp1(scale_r(1:NR),nu_axis(1:NR,1),r_nu_cont13);
    nu_min_rx=interp1(1:NR,nu_axis(:,1),pos_psi_rx);
else
    nu_min_cont13=interp1(scale_r(1:NR),nu_axis(1:NR,2),r_nu_cont13);
    nu_min_rx=interp1(1:NR,nu_axis(:,2),pos_psi_rx);
end

nu_min=min(nu_min_cont13,nu_min_rx);
nu_axis=nu_axis-nu_min;


pos_psi_middle=round(0.5*(pos_psi_rx+pos_r_nu_cont13));
if (omega_cont13==Nomega)
    pos_psi_middle=round(0.5*pos_psi_rx);
end

if nu_axis(ceil(pos_psi_rx),1)<nu_axis(pos_psi_middle,1)
    
    nu_max=max(nu_axis(pos_r_nu_cont13:ceil(pos_psi_rx),1));
    nu_axis=nu_max-nu_axis;
    nu_max_rx=interp1(1:NR,nu_axis(:,1),pos_psi_rx);
    
    nu_axis=(nu_axis/nu_max_rx)*area_region3;
    
else

    nu_max_cont13=interp1(scale_r(1:NR),nu_axis(1:NR,2),r_nu_cont13);
    nu_max_rx=interp1(1:NR,nu_axis(:,2),pos_psi_rx);
    nu_max=max(nu_max_cont13,nu_max_rx);
    
    nu_max_rx=interp1(1:NR,nu_axis(:,2),pos_psi_rx);
    
    nu_axis=(nu_axis/nu_max_rx)*area_region3;
    
end


[eps nu_min_pos]=min(nu_axis(1:ceil(pos_psi_rx),:));

if (omega_cont13==Nomega)
    nu_max=interp1(scale_r(1:NR),nu_axis(1:NR,2),r_nu_cont13,'pchip');
else
    nu_max=interp1(scale_r(1:NR),nu_axis(1:NR,1),r_nu_cont13,'pchip');
end



