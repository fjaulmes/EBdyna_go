function eps=evaluate_a1_coef_nu_axis(a1_coef,a2_coef,dPsih,Psih,psi_region1,ksi0,area_region3,x_nu_cont13,r_nu_cont13,x_nu_cont13_initial,r_nu_cont13_initial,omega_cont13,pos_psi_rx,rx,scale_r,scale_omega,NR,Nomega,alpha,nu_values,psi3_nu,Psih_limit13)

%%
pos_psi_rx_round=ceil(pos_psi_rx);

% Deriving the expression of flux contours in the bean shaped zone

der_Psih=interp1(1:NR-1,dPsih(1:NR-1),pos_psi_rx);
a1=a1_coef*(1/(r_nu_cont13_initial))*interp1(1:NR-1,dPsih(1:NR-1),x_nu_cont13_initial);
a2=a2_coef*(1/(rx))*abs(der_Psih);



run('calculate_nu_axis');

pos_r_nu_cont13=round(x_nu_cont13);
pos_r_nu_cont13_initial=round(x_nu_cont13_initial);
pos_shif_axis=pos_r_nu_cont13_initial-pos_r_nu_cont13;
%ksi0<0
r13_value=rx+2*ksi0;

run('rescale_nu_map_axis');
Dnu=(nu_max)/(size(nu_values,2)-1);

for (n=1:size(nu_values,2))
    nu_values(n)=(n-1)*Dnu;
end

for (r=1:NR)
    psi_star_3(r,1)=interp1(nu_values,psi3_nu,nu_axis(r,1),'pchip');
    psi_star_3(r,2)=interp1(nu_values,psi3_nu,nu_axis(r,2),'pchip');
end



der_psi=zeros(1,NR);

if (r13_value>=0)
    for(r=1:pos_r_nu_cont13-1)

        psi_star_3(r,1)=Psih(r+pos_shif_axis);
    end
else
    for(r=pos_r_nu_cont13+1:pos_psi_rx_round)

        psi_star_3(r,2)=psi_region1(r,Nomega);
    end
%    for(r=1:pos_r_nu_cont13_initial)

%        psi_star_3(pos_r_nu_cont13+pos_r_nu_cont13_initial+r,2)=Psih(r+1);
%    end
end

for(r=pos_psi_rx_round:NR)
    psi_star_3(r,1)=Psih(r);
    psi_star_3(r,2)=Psih(r);
end

if(omega_cont13==1)    
    for(r=2:pos_psi_rx_round+1)
        Dr_value=scale_r(r+1)-scale_r(r-1);
        der_psi(r)=(psi_star_3(r+1,1)-psi_star_3(r-1,1))/Dr_value;
    end
    psi_limit13=abs(interp1(scale_r(1:NR),psi_star_3(1:NR,1),r_nu_cont13,'pchip'));
else    
    for(r=2:pos_psi_rx_round+1)
        Dr_value=scale_r(r+1)-scale_r(r-1);
        der_psi(r)=-(psi_star_3(r+1,2)-psi_star_3(r-1,2))/Dr_value;
    end
    psi_limit13=abs(interp1(scale_r(1:NR),psi_star_3(1:NR,2),r_nu_cont13,'pchip'));
end

for(r=1:pos_r_nu_cont13)
    if (pos_shif_axis>0)
        der_psi(r)=dPsih(r+pos_shif_axis+1);
%    else
%        der_psi(r)=dPsih(r-pos_shif_axis+1);
    end
end

for(r=pos_psi_rx_round:NR-1)
    der_psi(r)=dPsih(r);
end

deriv_psi_limit13=abs(interp1(scale_r(1:NR-1),der_psi(1:NR-1),r_nu_cont13,'pchip'));
deriv_psi_limit13_initial=abs(interp1(scale_r(1:NR-1),dPsih(1:NR-1),r_nu_cont13_initial,'pchip'));

eps=abs(deriv_psi_limit13-deriv_psi_limit13_initial);
%eps=(deriv_psi_limit13-deriv_psi_limit13_initial)^2;

%for(r=2:pos_psi_rx+1)
%    der_psi(r)=(psi_star_3(r+1,1)-psi_star_3(r-1,1))/(2*Dr);
%end
    
%deriv_psi_limit13=abs(interp1(scale_r(1:pos_psi_rx+1),der_psi(1:pos_psi_rx+1),rx,'pchip'));
%deriv_psi_limit13_initial=abs(interp1(scale_r(1:NR),dPsih(1:NR),rx,'pchip'));

%eps=eps+(deriv_psi_limit13-deriv_psi_limit13_initial)^2;
