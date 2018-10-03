

    
%nu_min=interp1(scale_X(1:NR),nu(:,omega_cont13),pos_r_nu_cont13-1);
%nu_min=interp1(scale_X(1:NR),nu(1:NR,1),pos_psi_rx);
nu_min_cont13=interp1(1:NR,nu(1:NR,omega_cont13),x_nu_cont13);
nu_min_rx=interp1(1:NR,nu(:,omega_cont13),pos_psi_rx);
nu_min=min(nu_min_cont13,nu_min_rx);
nu=nu-nu_min;
%nu=max(nu,0);
%nu=sqrt(nu);

if (omega_cont13==1)
    pos_psi_middle=round(0.5*(pos_psi_rx+pos_r_nu_cont13));
else
    pos_psi_middle=round(0.5*pos_psi_rx);
end

if nu(xMixing_radius,1)<nu(pos_psi_middle,1)
    
    nu_max=max(nu(1:pos_psi_rx_round,1));
    nu=nu_max-nu;
    nu_max_rx=interp1(1:NR,nu(:,1),pos_psi_rx);
    
    nu=(nu/nu_max_rx)*area_region3;
    
else

    nu_max_cont13=interp1(1:NR,nu(1:NR,omega_cont13),x_nu_cont13);
    nu_max_rx=interp1(1:NR,nu(:,omega_cont13),pos_psi_rx);
    nu_max=max(nu_max_cont13,nu_max_rx);
    
    %nu=nu_max-nu;
    nu_max_rx=interp1(1:NR,nu(:,omega_cont13),pos_psi_rx);
    
    nu=(nu/nu_max_rx)*area_region3;
    
end

%nu_max=area_region3;


[eps nu_min_pos]=min(nu(1:pos_psi_rx_round,:));
%nu_real_max=interp1(1:pos_psi_rx,nu(1:pos_psi_rx,1),x_nu_cont13,'cubic');
%nu_max=interp1(1:NR,nu(1:NR,1),pos_psi_rx,'cubic');


nu_max_cont13=interp1(1:NR,nu(1:NR,omega_cont13),x_nu_cont13);
nu_max_rx=interp1(1:NR,nu(:,omega_cont13),pos_psi_rx);
nu_max=max(nu_max_cont13,nu_max_rx);
nu_max_inf=min(nu_max_cont13,nu_max_rx);


nu=(nu_max/nu_max_inf)*nu;

%for (omega=1:Nomega-1)
%for (r=1:nu_min_pos(omega))
%     nu(r,omega)=nu(r,omega)*nu_max/nu_real_max;
%end
%end


