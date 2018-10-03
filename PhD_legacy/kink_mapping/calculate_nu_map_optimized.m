% Deriving the expression of flux contours in the bean shaped zone

%NU_OFFSET=0;

if (f<=round(0.4*f_final))
    f_range=(0.12*f+4)/4;
    coef_ratio_sat=1.25+0.002*f;
elseif (f<=round(0.5*f_final))
    f_range=((0.12/PRECISE_MESH)*f+4)/4;
    coef_ratio_sat=1.27+0.005*f;
elseif (f<=round(0.9*f_final))
    f_range=((0.12/PRECISE_MESH)*f+4)/4;
    coef_ratio_sat=1.27+0.008*f;
else
    f_range=((0.25/PRECISE_MESH)*f+4)/4;
%     f_range=(0.05*f_final+18)/4;
    coef_ratio_sat=1.27+0.011*f;
end

a1_coef_min=max(0+0.7*f_range,0.5);
a1_coef_max=min(0.8+4.0*f_range,3e2);
%A1_INF_FACTOR=0.5;

a2_coef_min=max(0+0.6*f_range,0.5);
a2_coef_max=min(0.6+4.0*f_range,3e2);

a1_coef_max=max(a1_coef_max,a2_coef_max);


%FRAME_LIMIT_SUP_FOR_CONTINUITY=5;

if ((f>=FRAME_LIMIT_INF_FOR_CONTINUITY)&&(f<FRAME_LIMIT_SUP_FOR_CONTINUITY))
    [a1_coef eps]=fminbnd(@(a1_coef) evaluate_a1_coef_nu_axis(a1_coef_prev,a2_coef_prev,dPsih,Psih,psi_region1,sign_ksi0,area_region3,x_nu_cont13,r_nu_cont13,x_nu_cont13_initial,r_nu_cont13_initial,omega_cont13,pos_psi_rx,rx,scale_r,scale_omega,NR,Nomega,alpha,nu_values,psi3_nu,Psih_limit13),a1_coef_min,a1_coef_max,...
        optimset('TolX',1e-1,'Display','off'));
    a1_coef_raw=a1_coef;
    %[a2_coef eps]=fminbnd(@(a2_coef) evaluate_a2_coef_nu_axis(a1_coef,a2_coef_prev,dPsih,Psih,psi_region1,sign_ksi0,area_region3,x_nu_cont13,r_nu_cont13,x_nu_cont13_initial,r_nu_cont13_initial,omega_cont13,pos_psi_rx,rx,scale_r,scale_omega,NR,Nomega,alpha,nu_values,psi3_nu,Psih_limit13),a2_coef_min,a2_coef_max);
    %a1_coef=max(a1_coef,a2_coef);
%     a1_coef=0.25*a1_coef+0.75*a1_coef_prev;
%     a2_coef=0.5*(a2_coef+a2_coef_prev);

    %[a1_coef eps]=fminbnd(@(a1_coef) evaluate_a1_coef_nu_axis(a1_coef,a2_coef,dPsih,Psih,psi_region1,sign_ksi0,area_region3,x_nu_cont13,r_nu_cont13,x_nu_cont13_initial,r_nu_cont13_initial,omega_cont13,pos_psi_rx,rx,scale_r,scale_omega,NR,Nomega,alpha,nu_values,psi3_nu,Psih_limit13),a2_coef,1.1*a1_coef_max);
    [a2_coef eps]=fminbnd(@(a2_coef) evaluate_a2_coef_nu_axis(a1_coef,a2_coef,dPsih,Psih,psi_region1,sign_ksi0,area_region3,x_nu_cont13,r_nu_cont13,x_nu_cont13_initial,r_nu_cont13_initial,omega_cont13,pos_psi_rx,rx,scale_r,scale_omega,NR,Nomega,alpha,nu_values,psi3_nu,Psih_limit13),a2_coef_min,a2_coef_max,...
        optimset('TolX',1e-1,'Display','off'));
    if (f<=round(0.5*f_final))
        a1_coef=max(a1_coef,a2_coef);
    end
    
    a1_coef=max(a1_coef,a1_coef_raw);
%     a2_coef=0.5*(a2_coef+a2_coef_prev);
    if (f<=round(0.5*f_final))
        a2_coef=0.5*a2_coef+0.5*a2_coef_prev;
    else
        a2_coef=0.5*a2_coef+0.5*a2_coef_prev;
    end
%a1_coef=a2_coef;
%     [a1_coef eps]=fminbnd(@(a1_coef) evaluate_a1_coef_nu_axis(a1_coef,a2_coef,dPsih,Psih,psi_region1,sign_ksi0,area_region3,x_nu_cont13,r_nu_cont13,x_nu_cont13_initial,r_nu_cont13_initial,omega_cont13,pos_psi_rx,rx,scale_r,scale_omega,NR,Nomega,alpha,nu_values,psi3_nu,Psih_limit13),a1_coef_min,a1_coef_max);
%     
    [a1_coef eps]=fminbnd(@(a1_coef) evaluate_a1_coef_nu_axis(a1_coef,a2_coef,dPsih,Psih,psi_region1,sign_ksi0,area_region3,x_nu_cont13,r_nu_cont13,x_nu_cont13_initial,r_nu_cont13_initial,omega_cont13,pos_psi_rx,rx,scale_r,scale_omega,NR,Nomega,alpha,nu_values,psi3_nu,Psih_limit13),0.9*a1_coef_raw,a1_coef_max,...
    optimset('TolX',1e-2,'Display','off'));
%     [a2_coef eps]=fminbnd(@(a2_coef) evaluate_a2_coef_nu_axis(a1_coef,a2_coef,dPsih,Psih,psi_region1,sign_ksi0,area_region3,x_nu_cont13,r_nu_cont13,x_nu_cont13_initial,r_nu_cont13_initial,omega_cont13,pos_psi_rx,rx,scale_r,scale_omega,NR,Nomega,alpha,nu_values,psi3_nu,Psih_limit13),a2_coef_min,a2_coef_max,...    
%     optimset('TolX',1e-11,'Display','off'));

    if (f<=round(0.5*f_final))
        a1_coef=0.5*a1_coef+0.5*a1_coef_prev;
    else
        a1_coef=0.5*a1_coef+0.5*a1_coef_prev;
    end
%     a2_coef=min(a1_coef,a2_coef);

        
%     if abs(a1_coef-a2_coef)<0.35
%         [a1_coef eps]=fminbnd(@(a1_coef) evaluate_a1_coef_nu_axis(a1_coef,a2_coef,dPsih,Psih,psi_region1,sign_ksi0,area_region3,x_nu_cont13,r_nu_cont13,x_nu_cont13_initial,r_nu_cont13_initial,omega_cont13,pos_psi_rx,rx,scale_r,scale_omega,NR,Nomega,alpha,nu_values,psi3_nu,Psih_limit13),1.01*a2_coef,1.8*a1_coef_max,...
%             optimset('TolX',1e-11,'Display','off'));
%         %[a2_coef eps]=fminbnd(@(a2_coef) evaluate_a2_coef_nu_axis(a1_coef,a2_coef,dPsih,Psih,psi_region1,sign_ksi0,area_region3,x_nu_cont13,r_nu_cont13,x_nu_cont13_initial,r_nu_cont13_initial,omega_cont13,pos_psi_rx,rx,scale_r,scale_omega,NR,Nomega,alpha,nu_values,psi3_nu,Psih_limit13),a2_coef_min,1.2*a2_coef_max);
%     end
%     a1_coef=max(a1_coef,a2_coef);
%     if abs(a1_coef-a2_coef)<0.1
%         [a1_coef eps]=fminbnd(@(a1_coef) evaluate_a1_coef_nu_axis(a1_coef,a2_coef,dPsih,Psih,psi_region1,sign_ksi0,area_region3,x_nu_cont13,r_nu_cont13,x_nu_cont13_initial,r_nu_cont13_initial,omega_cont13,pos_psi_rx,rx,scale_r,scale_omega,NR,Nomega,alpha,nu_values,psi3_nu,Psih_limit13),a2_coef,1.3*a1_coef_max);
%         %[a2_coef eps]=fminbnd(@(a2_coef) evaluate_a2_coef_nu_axis(a1_coef,a2_coef,dPsih,Psih,psi_region1,sign_ksi0,area_region3,x_nu_cont13,r_nu_cont13,x_nu_cont13_initial,r_nu_cont13_initial,omega_cont13,pos_psi_rx,rx,scale_r,scale_omega,NR,Nomega,alpha,nu_values,psi3_nu,Psih_limit13),a2_coef_min,1.3*a2_coef_max);
%     end
%     a1_coef=max(a1_coef,a2_coef);

    
    %[a2_coef eps]=fminbnd(@(a2_coef) evaluate_a2_coef_nu_axis(a1_coef,a2_coef,dPsih,Psih,psi_region1,sign_ksi0,area_region3,x_nu_cont13,r_nu_cont13,x_nu_cont13_initial,r_nu_cont13_initial,omega_cont13,pos_psi_rx,rx,scale_r,scale_omega,NR,Nomega,alpha,nu_values,psi3_nu,Psih_limit13),a2_coef_min,1.2*a2_coef_max);
    [a1_coef eps]=fminbnd(@(a1_coef) evaluate_a1_coef_nu_axis(a1_coef,a2_coef,dPsih,Psih,psi_region1,sign_ksi0,area_region3,x_nu_cont13,r_nu_cont13,x_nu_cont13_initial,r_nu_cont13_initial,omega_cont13,pos_psi_rx,rx,scale_r,scale_omega,NR,Nomega,alpha,nu_values,psi3_nu,Psih_limit13),0.8*a1_coef_raw,a1_coef_max,...
    optimset('TolX',1e-3,'Display','off'));
    [a2_coef eps]=fminbnd(@(a2_coef) evaluate_a2_coef_nu_axis(a1_coef,a2_coef,dPsih,Psih,psi_region1,sign_ksi0,area_region3,x_nu_cont13,r_nu_cont13,x_nu_cont13_initial,r_nu_cont13_initial,omega_cont13,pos_psi_rx,rx,scale_r,scale_omega,NR,Nomega,alpha,nu_values,psi3_nu,Psih_limit13),0.8*a2_coef,a2_coef_max,...
        optimset('TolX',1e-3,'Display','off'));
    
    a1_coef=0.9*a1_coef+0.1*a1_coef_prev;
    a2_coef=0.9*a2_coef+0.1*a2_coef_prev;
%     [a2_coef eps]=fminbnd(@(a2_coef) evaluate_a2_coef_nu_axis(a1_coef,a2_coef,dPsih,Psih,psi_region1,sign_ksi0,area_region3,x_nu_cont13,r_nu_cont13,x_nu_cont13_initial,r_nu_cont13_initial,omega_cont13,pos_psi_rx,rx,scale_r,scale_omega,NR,Nomega,alpha,nu_values,psi3_nu,Psih_limit13),0.8*a2_coef,a2_coef_max,...
%         optimset('TolX',1e-5,'Display','off'));
    
    if abs(a2_coef/a1_coef)>coef_ratio_sat
        disp('coef ratio saturated a1')
        a1_coef=(1/coef_ratio_sat)*a2_coef;
    end
    
    if abs(a1_coef/a2_coef)>coef_ratio_sat
        disp('coef ratio saturated a2')
        a2_coef=(1/coef_ratio_sat)*a1_coef;
    end
    
%     if (a1_coef/a2_coef)>2
%         [a1_coef eps]=fminbnd(@(a1_coef) evaluate_a1_coef_nu_axis(a1_coef,a2_coef,dPsih,Psih,psi_region1,sign_ksi0,area_region3,x_nu_cont13,r_nu_cont13,x_nu_cont13_initial,r_nu_cont13_initial,omega_cont13,pos_psi_rx,rx,scale_r,scale_omega,NR,Nomega,alpha,nu_values,psi3_nu,Psih_limit13),a2_coef,a1_coef_max);
%         [a2_coef eps]=fminbnd(@(a2_coef) evaluate_a2_coef_nu_axis(a1_coef,a2_coef,dPsih,Psih,psi_region1,sign_ksi0,area_region3,x_nu_cont13,r_nu_cont13,x_nu_cont13_initial,r_nu_cont13_initial,omega_cont13,pos_psi_rx,rx,scale_r,scale_omega,NR,Nomega,alpha,nu_values,psi3_nu,Psih_limit13),a2_coef_min,a2_coef_max);
%     end
%     if f>0.5*(f_final-f_initial)+f_initial
%         [a2_coef eps]=fminbnd(@(a2_coef) evaluate_a2_coef_nu_axis(a1_coef,a2_coef,dPsih,Psih,psi_region1,sign_ksi0,area_region3,x_nu_cont13,r_nu_cont13,x_nu_cont13_initial,r_nu_cont13_initial,omega_cont13,pos_psi_rx,rx,scale_r,scale_omega,NR,Nomega,alpha,nu_values,psi3_nu,Psih_limit13),a2_coef_min,1.15*a2_coef_max);
%     end
    %[a2_coef eps]=fminbnd(@(a2_coef) evaluate_a2_coef_nu_axis(a1_coef,a2_coef,dPsih,Psih,psi_region1,sign_ksi0,area_region3,x_nu_cont13,r_nu_cont13,x_nu_cont13_initial,r_nu_cont13_initial,omega_cont13,pos_psi_rx,rx,scale_r,scale_omega,NR,Nomega,alpha,nu_values,psi3_nu,Psih_limit13),a2_coef_min,a2_coef_max);

else
    a1_coef=a1_coef_prev;
    a2_coef=a2_coef_prev;
end
    
der_Psih=interp1(1:NR-1,dPsih(1:NR-1),pos_psi_rx);
a1=a1_coef*(1/(r_nu_cont13_initial))*interp1(1:NR-1,dPsih(1:NR-1),x_nu_cont13_initial);
a2=a1_coef*(1/(rx))*abs(der_Psih);

disp('nu profile optimized for derivative continuity...');
disp('a1_coef = ');
disp(a1_coef);
disp('a2_coef = ');
disp(a2_coef);
% disp('a1 = ');
% disp(a1);
% disp('a2 = ');
% disp(a2);





