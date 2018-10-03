NBFRAME=101
time_evol=zeros(1,NBFRAME);
time_exponential_evol=zeros(1,NBFRAME);
ksi0_evol=zeros(1,NBFRAME);
rx_evol=zeros(1,NBFRAME)+r_value_q1_mean;

max_psi3=max(Psih_final);

ksi0_prev=0;

psi_star_2D_evol=zeros(NBFRAME,NR,Nomega);

time0_offset=1/(NBFRAME-1);
f_rank=1;
psi_star_2D_evol(1,:,:)=psi_star_2D;
% rx_evol(1)=r_value_q1_mean;

f_rank=2;
psi_star_2D_evol(2,:,:)=psi_star_2D;
% rx_evol(2)=r_value_q1_mean;



rPsih_zero=interp1(1:NR,scale_r,xPsih_zero);
r_q1=interp1(1:NR,scale_r,xPsih_max);
rmix=rPsih_zero;

   
Psih_max=max(Psih);
psi_contour_values=Psih_max*([0:22])/22;


f_initial=2;
f_final=199;

f_size=f_final-f_initial+1;


Df=(f_final*0.01-time0_offset)/(f_final-f_initial+2);

time_exponential_evol=(0:NBFRAME-1)*time0_offset;

% up to 2cm dispalcement
exponential_ksi0_evol=0.0008*exp(time_exponential_evol*log(100));
exponential_ksi0_evol(1)=0;

rmix=rPsih_zero;
rx=r_q1;

pos_rmix=xPsih_zero;
pos_psi_rx=xPsih_max;
pos_psi_rx_round=round(pos_psi_rx);

area_region3=interp1((1:NR),surf_flux_precise(1:NR),pos_psi_rx,'cubic');
omega_cont13=0

%alpha=atan(2);
for (f=f_initial:f_final)
    %%
    psi_star_2D=zeros(NR,Nomega);

    f_rank=f;
    %close(6);
    
    disp('----------------');
    clear S3 S3_filt;
    clear psi3_nu;
    
    if (f<10)
        frame_name='0';
    else
        frame_name='';
    end
    frame_name=strcat(frame_name,num2str(f));
    psi_rx=interp1((floor(xPsih_max):ceil(xMixing_radius)),Psih(floor(xPsih_max):ceil(xMixing_radius)),pos_psi_rx);
    x_nu_cont13_initial=interp1(Psih(1:ceil(xPsih_max)),(1:ceil(xPsih_max)),psi_rx,'cubic');
    r_nu_cont13_initial=interp1(Psih(1:ceil(xPsih_max)),scale_r(1:ceil(xPsih_max)),psi_rx,'cubic');
    
%     alpha=(Df*(f_final-f+2)*(1-atan(2)/pi)+atan(2)/pi)*pi;
    alpha=((1/f_final)*(f_final+(f_initial-f))*(1-atan(2)/pi)+atan(2)/pi)*pi;
	Dalpha=(pi-alpha)/pi
%     ksi0=-0.5*rx*tan(alpha);
%     sign_ksi0=sign(ksi0);

    rx_evol(f_rank)=rx;
%     ksi0_evol(f_rank)=abs(ksi0);
    disp('rx=');
    disp(rx);
    
    % Deriving the expression of flux contours in the bean shaped zone

    a1_coef=1;
    a2_coef=1;
    calculate_nu_map;
    nu_max=max(max(nu));
    nu=nu_max-nu;
    
    N_NU=100;
    Dnu=(nu_max)/(N_NU-1);
    
    for (n=1:N_NU)
        nu_values(n)=(n-1)*Dnu;
    end
   
  
    %-----------------------------------------
    % Helical flux in the magnetic island
    %-----------------------------------------    
    min_psi3=0;
    max_psi3=psi_star_max;
    pos_psi_rx_round=round(pos_psi_rx);
    calculate_nu_surfaces;
    
    for (r=1:NR)
        psi_region3(r,:)=interp1(nu_values,psi3_nu,nu(r,:),'cubic');
    end
    psi_region3(isnan(psi_region3))=0;    

    
    for (r=1:NR)
        r_value=scale_r(r);
        for (omega=1:Nomega)
            omega_value=scale_omega(omega);
            if(r_value < rx)
                psi_star_2D(r,omega)=psi_region3(r,omega);
            elseif (r_value == rx)
                psi_star_2D(r,omega)=0.5*(psi_region3(r,omega)+Psih(r));
            else
                psi_star_2D(r,omega)=Psih(r);
            end
        end
    end
    
    psi_star_2D_Xder=gradient(psi_star_2D(1:round((Dalpha+0.05)*xPsih_max)-1,1),scale_r(1:round((Dalpha+0.05)*xPsih_max)-1));
    
    ksi0_pos=interp1(psi_star_2D_Xder,1:round((Dalpha+0.05)*xPsih_max)-1,0,'linear');

    if f>f_initial
        f_rank
		ksi0=interp1(1:NR,scale_r,ksi0_pos,'cubic')
		sign_ksi0=sign(ksi0);
		ksi0_evol(f_rank)=ksi0;
        time_evol(f_rank)=interp1(exponential_ksi0_evol,time_exponential_evol,ksi0,'cubic');
	else
		ksi0=0
    end
    
    psi_star_2D_evol(f_rank,:,:)=psi_star_2D;
    
    % Create grids and convert polar coordinates to rectangular
    [THETA,RR] = meshgrid(scale_omega,scale_r);
    [XX,YY] = pol2cart(THETA,RR);
    
    
    figure(3)
    set(gca,'fontsize',20)
    %surfc(XX,YY,psi_star_2D,'edgecolor','none');view(0,90);
    contour(XX,YY,psi_star_2D,psi_contour_values,'linewidth',2);
    
    axis xy equal;
    %imagesc(scale_r,scale_omega,Psi_star_2D)
    grid on;
    pause(0.22);
    
%     F=getframe;
%     [im,map] = frame2im(F);    %Return associated image data
%     
%     imwrite(im,filename,'bmp');    

    
    disp('...done...');
    %pause

    
end

disp('-----')
disp('fixing time point 2 ')
disp('-----')
%ksi0_evol(4)=0.5*(ksi0_evol(3)+ksi0_evol(5));
%ksi0_evol(3)=0.5*(ksi0_evol(2)+ksi0_evol(4));
%time_evol(4)=3*(time_evol(5))/4;
%time_evol(3)=2*(time_evol(5))/4;
time_evol(f_initial)=(time_evol(f_initial+1))/2;


% figure(8);grid on;hold on;
% 
% set(gca,'FontSize',18);
% 
% plot(radial_r_value, P_initial_profile,'b--');
% plot(radial_r_value, P_final_profile,'r');
% 
% set( findobj(gcf,'Type','line'),'LineWidth',1.5);
% ylabel('pressure (Pa)')
% xlabel('r (m)')
% xlim([0 1]);
% 
% legend('before collapse','afer collapse');

