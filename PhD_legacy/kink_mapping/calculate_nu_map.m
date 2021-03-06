% Deriving the expression of flux contours in the bean shaped zone


NU_OFFSET=1e-4;
    
der_Psih=interp1(1:NR-1,dPsih(1:NR-1),pos_psi_rx,'cubic');
% a1=a1_coef*(1/(r_nu_cont13_initial))*interp1(1:NR-1,dPsih(1:NR-1),1);
a2=a2_coef*(1/(rx))*abs(der_Psih);
% a2=1

a1=a2;


nu=zeros(NR,Nomega);
z1=zeros(NR,Nomega);
z2=zeros(NR,Nomega);


for (r=1:min(round(pos_psi_rx),NR))
    r_value=scale_r(r);
    rv2=r_value^2;
    for (omega=1:Nomega)
        omega_value=scale_omega(omega);
        z1(r,omega)=a1*((rv2-rx^2)*cos(alpha)+rx*(rx+r_value*cos(omega_value))*sin(alpha));
        z2(r,omega)=a2*(rx^2-rv2);
    end
end

den=z1+z2+NU_OFFSET;
nu=z1.*z2./(den);
nu(isnan(nu))=0;

% if(omega_cont13==1)
%     for (r=1:pos_r_nu_cont13-2)
%         nu(r,:)=nu(pos_r_nu_cont13-1,:);
%     end
% end

% for (r=pos_psi_rx_round+2:NR)
%     nu(r,:)=interp1(1:NR,nu(:,:),pos_psi_rx+1);
% end



