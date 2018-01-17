function eps=evaluate_a2_coef_nu_axis(a1_coef,a2_coef,dPsih,Psih,psi_region1,sign_ksi0,area_region3,x_nu_cont13,r_nu_cont13,x_nu_cont13_initial,r_nu_cont13_initial,omega_cont13,pos_psi_rx,rx,scale_r,scale_omega,NR,Nomega,alpha,nu_values,psi3_nu,Psih_limit13)

% Deriving the expression of flux contours in the bean shaped zone


der_Psih=interp1(1:NR-1,dPsih(1:NR-1),pos_psi_rx);

a1=a1_coef*(1/(r_nu_cont13_initial))*interp1(1:NR-1,dPsih(1:NR-1),x_nu_cont13_initial);
a2=a2_coef*(1/(rx))*abs(der_Psih);



run('calculate_nu_axis');


pos_r_nu_cont13=round(x_nu_cont13);
pos_r_nu_cont13_initial=round(x_nu_cont13_initial);
pos_shif_axis=pos_r_nu_cont13_initial-pos_r_nu_cont13;
run('rescale_nu_map_axis');
Dnu=(nu_max)/(size(nu_values,2)-1);

for (n=1:size(nu_values,2))
    nu_values(n)=(n-1)*Dnu;
end

for (r=1:NR)
    psi_star_3(r,1)=interp1(nu_values,psi3_nu,nu_axis(r,1),'cubic');
    psi_star_3(r,2)=interp1(nu_values,psi3_nu,nu_axis(r,2),'cubic');
end



der_psi=zeros(1,NR);

if (pos_shif_axis>0)
    for(r=1:pos_r_nu_cont13-1)

        psi_star_3(r,1)=Psih(r+pos_shif_axis+1);
    end
end

pos_psi_rx_round=ceil(pos_psi_rx);
for(r=pos_psi_rx_round:NR)
    psi_star_3(r,1)=Psih(r);
    psi_star_3(r,2)=Psih(r);
end


for(r=pos_psi_rx_round:NR-1)
    der_psi(r)=dPsih(r);
end


for(r=2:pos_psi_rx_round)
    Dr_value=scale_r(r+1)-scale_r(r-1);
    der_psi(r)=(psi_star_3(r+1,1)-psi_star_3(r-1,1))/Dr_value;
end
    
deriv_psi_limit32=abs(interp1(scale_r(1:NR-1),der_psi(1:NR-1),rx,'*linear'));
deriv_psi_limit32_initial=abs(interp1(scale_r(1:NR-1),dPsih(1:NR-1),rx,'*linear'));

eps=abs(deriv_psi_limit32-deriv_psi_limit32_initial);
%eps=(deriv_psi_limit32-deriv_psi_limit32_initial)^2;

