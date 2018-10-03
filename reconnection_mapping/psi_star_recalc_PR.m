
FIRST_PHASE_FACTOR=1.4

% psi_star_initial;

%pause(1);
% close all


%r1=interp1(scale_X,scale_r_values(xPsih_max));
x_initial_psimax_pos=interp1(scale_r,(1:NR),r_q1);
max_psi3=max(Psih_final);

ksi0_prev=0;
a1_coef_prev=1.2;
a2_coef_prev=1;
%total_psi3_flux(1)=0;

f_initial=1
f_final=round((xPsih_zero-xPsih_max))-2
f_stop_gradient_continuity=f_final+5;

f_size=f_final-f_initial+3;

psi_star_2D_evol=zeros(301,NR,Nomega);
time_evol=zeros(1,f_size);
ksi0_evol=zeros(1,f_size);
rx_evol=zeros(1,f_size);

time0_offset=1/101;
f_rank=1;
psi_star_2D_evol(1,:,:)=psi_star_2D;
rx_evol(1)=r_q1;
time_evol(f_rank)=0;

f_rank=2;
psi_star_2D_evol(2,:,:)=psi_star_2D;
rx_evol(2)=r_q1;
time_evol(f_rank)=(f_rank-1)*time0_offset;

% f_rank=3;
% psi_star_2D_evol(3,:,:)=psi_star_2D;
% rx_evol(3)=r_q1;
% rx_prev=r_q1;
% time_evol(f_rank)=f_rank*time0_offset;

time0_offset=(f_rank-1)*time0_offset

rPsih_zero=interp1(1:NR,scale_r,xPsih_zero);
rmix=rPsih_zero;

volume1_prev=2*pi*(R0+X_axis)*interp1((1:NR),surf_flux_precise(1:NR),xPsih_max,'pchip')
volume2_prev=2*pi*(R0+X_axis)*interp1((1:NR),surf_flux_precise(1:NR),xPsih_max,'pchip')
    
Psih_max=max(Psih);
psi_contour_values=Psih_max*([0:22])/22;

% figure(1);
% hold on;
% grid on;
% set(gca,'FontSize',18);
% 
% %if (f-2*round(0.5*f))==0
%     plot([-(size_r-1:-1:0)  (1:size_r-1)]/(xPsih_zero-1),[psi_star_2D(size_r:-1:1,Nomega) ; psi_star_2D(2:size_r,1)],'r','LineWidth',2);
% %end
% xlabel('r/r_{mix}')
% ylabel('\Psi _*')
% axis([-1 1 -0.005 0.02])

%alpha=atan(2);
f_transition=f_final

% from Xray tomography

%poly_evol=[-0.1 -0.12 0.3 0.2 0]
%time_scale_poly=(0:1000)/1000;
%poly_evol_values=polyval(poly_evol,time_scale_poly);

%rescaled to half of the time values
%add the initial offset to avoid an initial bump
%time_scale_poly=0.5*time_scale_poly+time0_offset;

for (f=1:f_final)
    %%
    psi_star_2D=zeros(NR,Nomega);
    
    f_rank=f+2;
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
%     alpha=0.01*(f)*pi;
%     filename='reconnection_movie\t0';
%     filename=strcat(filename,frame_name,'.bmp');
%     disp(filename);
%     time_alpha_value=0.01*(f-f_initial)/FIRST_PHASE_FACTOR;
%     alpha=pi*time_alpha_value
%     time_value=FIRST_PHASE_FACTOR*alpha/pi+time0_offset
%     ksi0=0.5*tan(alpha)*rx
    
    pos_psi_rx=ceil(x_initial_psimax_pos+f);
    rx=interp1(1:NR,scale_r,pos_psi_rx,'pchip');
    psi_rx=interp1((floor(xPsih_max):ceil(xMixing_radius)),Psih(floor(xPsih_max):ceil(xMixing_radius)),pos_psi_rx,'pchip');
    x_nu_cont13_initial=interp1(Psih(1:ceil(xPsih_max)),1:ceil(xPsih_max),psi_rx,'pchip');
    r_nu_cont13_initial=interp1(Psih(1:ceil(xPsih_max)),scale_r(1:ceil(xPsih_max)),psi_rx,'pchip');
    ksi0=r_nu_cont13_initial-rx;
    sign_ksi0=sign(ksi0);

    rx_evol(f_rank)=rx;
    ksi0_evol(f_rank)=abs(ksi0);
    disp('rx=');
    disp(rx);
    alpha=atan(2*abs(ksi0)/rx);
    %time_value=interp1(poly_evol_values,time_scale_poly,abs(ksi0))
    time_value=FIRST_PHASE_FACTOR*alpha/pi+time0_offset
    time_evol(f_rank)=time_value;
    disp('time=');
    disp(time_evol(f_rank));
    
    x_psi_final_rx=interp1(Psih_final(1:xMixing_radius), scale_X_precise(1:xMixing_radius),psi_rx,'pchip');
    %r_psi_final_rx=interp1((1:xMixing_radius),scale_r(1:xMixing_radius),x_psi_final_rx);

    pos_psi_rx_round=ceil(pos_psi_rx);
    
    
    psi_region1=zeros(NR,Nomega);
    map_region1=zeros(NR,Nomega);
    psi_region2=zeros(NR,Nomega);
    psi_region3=zeros(NR,Nomega);
    
    %-----------------------------------------
    % Helical flux in the moving core region
    %-----------------------------------------
    r_nu_cont13=rx-2*abs(ksi0);

    if (r_nu_cont13>=0)
        x_nu_cont13=interp1(scale_r,(1:NR),r_nu_cont13,'pchip')
        pos_r_nu_cont13=round(x_nu_cont13);
        omega_cont13=1;
    else
        r_nu_cont13=abs(r_nu_cont13);
        x_nu_cont13=interp1(scale_r,(1:NR),r_nu_cont13,'pchip')
        pos_r_nu_cont13=round(x_nu_cont13);
        omega_cont13=Nomega;
    end
    
    for (r=1:pos_psi_rx)
        for (omega=1:Nomega)
            r_value=scale_r(r);
            omega_value=scale_omega(omega);
            r_calc_region1_sq=(((r_value)*cos(omega_value)-(ksi0))^2+(r_value^2)*sin(omega_value)^2);
            r_calc_region1=sqrt(r_calc_region1_sq);
            if (r_calc_region1>rmix)
                r_calc_region1=rmix;
            end
            psi_region1(r,omega)=interp1(scale_r,Psih(1:NR),r_calc_region1,'pchip');
            if (r_calc_region1_sq > (rx-abs(ksi0)+0.5*Dr_avg(round(x_nu_cont13)+1))^2)
                map_region1(r,omega)=0;
            elseif (r_calc_region1_sq > (rx-abs(ksi0))^2)
                map_region1(r,omega)=2;
            else
                map_region1(r,omega)=1;
            end
        end
    end
    
    Psih_limit23=psi_rx;
    


    
    %x_nu_cont13_initial=x_psi_limit13-abs(ksi0)/Dr;
    
    %r_nu_cont13_initial=rx-abs(ksi0);
    %x_nu_cont13_initial=r_nu_cont13_initial/Dr+1;
    pos_r_nu_cont13_initial=round(x_nu_cont13_initial);
    
    calculate_island_area;
    volume1=2*pi*(R0+X_axis)*interp1((1:NR),surf_flux_precise(1:NR),x_nu_cont13_initial,'pchip')
    volume2=2*pi*(R0+X_axis)*interp1((1:NR),surf_flux_precise(1:NR),pos_psi_rx,'pchip')    
    volume3=2*pi*(R0+X_axis)*area_region3;
    Bstar1=interp1((1:NR),Bstar_initial_ext(1:NR),x_nu_cont13_initial,'pchip')
    Bstar2=interp1((1:NR),Bstar_initial_ext(1:NR),pos_psi_rx,'pchip')    
    Bstar3=interp1((1:NR),Bstar_final_ext(1:NR),x_psi_final_rx,'pchip')    
    delta_volume1=abs(volume1-volume1_prev);
    delta_volume2=abs(volume2-volume2_prev);
    gamma=delta_volume2/delta_volume1;
    P1=interp1((1:NR),P_ini_precise(1:NR),x_nu_cont13_initial,'pchip');
    P2=interp1((1:NR),P_ini_precise(1:NR),pos_psi_rx,'pchip');
    P3=P2/(1+1/gamma)+P1/(1+gamma);
    PI1=interp1((1:NR),PI_ini_precise(1:NR),x_nu_cont13_initial,'pchip');
    PI2=interp1((1:NR),PI_ini_precise(1:NR),pos_psi_rx,'pchip');
    PI3=PI2/(1+1/gamma)+PI1/(1+gamma);
    PE1=interp1((1:NR),PE_ini_precise(1:NR),x_nu_cont13_initial,'pchip');
    PE2=interp1((1:NR),PE_ini_precise(1:NR),pos_psi_rx,'pchip');
    PE3=PE2/(1+1/gamma)+PE1/(1+gamma);
    volume1_evol(f)=volume1;
    volume2_evol(f)=volume2;
    volume3_evol(f)=volume3;
    Bstar1_evol(f)=Bstar1;
    Bstar2_evol(f)=Bstar2;
    Bstar3_evol(f)=Bstar3;
	P1_evol(f)=P1;
    P2_evol(f)=P2;
    P3_evol(f)=P3;
	PI1_evol(f)=PI1;
    PI2_evol(f)=PI2;
    PI3_evol(f)=PI3;
	PE1_evol(f)=PE1;
    PE2_evol(f)=PE2;
    PE3_evol(f)=PE3;
    P_final_profile_raw(f)=P3;
    pos_P_final_profile(f)=x_psi_final_rx;
    

    a1_coef=1.2;
    a2_coef=1.2;


    calculate_nu_map;
    
    
    %determining N_NU arbitrary steps
    %for scaling the flux in nu scale
    N_NU=12+f;
    nu_values=zeros(1,N_NU);
    
    rescale_nu_map_alt;
    
    Dnu=(nu_max)/(N_NU-1);
    
    for (n=1:N_NU)
        nu_values(n)=(n-1)*Dnu;
    end
    
    for (omega=1:Nomega)
        psi_region2(1:NR,omega)=Psih(1:NR);
    end

    
    % continuity of psi
    
    Psih_limit13=psi_rx;
    deriv_psi_limit13=abs(interp1((1:xMixing_radius),dPsih(1:xMixing_radius),x_nu_cont13_initial,'pchip'));
    min_psi3=Psih_limit13;
    
    
    
    %-----------------------------------------
    % Helical flux in the expanding magnetic island
    %-----------------------------------------
    FRAME_LIMIT_INF_FOR_CONTINUITY=4;
    FRAME_LIMIT_SUP_FOR_CONTINUITY=f_stop_gradient_continuity;
    
    if ((f>=FRAME_LIMIT_INF_FOR_CONTINUITY)&&(f<FRAME_LIMIT_SUP_FOR_CONTINUITY))
        calculate_nu_surfaces;
        calculate_nu_map_optimized;
        calculate_nu_map;
        rescale_nu_map_alt;
        Dnu=(nu_max)/(N_NU-1);
        
        for (n=1:N_NU)
            nu_values(n)=(n-1)*Dnu;
        end
%     else
%         a1_coef=1.1*a1_coef_prev
%         a2_coef=1.1*a2_coef_prev;
    end
    
    calculate_nu_surfaces;


%     run('calculate_nu_map_optimized');    
%     run('calculate_nu_map');
%     run('rescale_nu_map_alt');
%     Dnu=(nu_max)/(N_NU-1);
%     
%     for (n=1:N_NU)
%         nu_values(n)=(n-1)*Dnu;
%     end

    

%    if (f>f_initial)
        for (r=1:NR)
            psi_region3(r,:)=interp1(nu_values,psi3_nu,nu(r,:),'pchip');
        end
%    else
%        psi_region3=psi_region3+psi_star_max;
%    end



    psi_region3(isnan(psi_region3))=0;
    
    
    for (r=1:NR)
        r_value=scale_r(r);
        for (omega=1:Nomega)
            omega_value=scale_omega(omega);
            if (map_region1(r,omega)==1)
                psi_star_2D(r,omega)=psi_region1(r,omega);
            elseif (map_region1(r,omega)==2)
                 psi_star_2D(r,omega)=0.8*psi_region1(r,omega)+0.2*psi_region3(r,omega);
            else
                if(r <= pos_psi_rx)
                    psi_star_2D(r,omega)=psi_region3(r,omega);
%                 elseif (r == pos_psi_rx-1)
%                     psi_star_2D(r,omega)=0.6*psi_region3(r,omega)+0.4*Psih(r);
%                 elseif (r == pos_psi_rx)
%                     psi_star_2D(r,omega)=0.5*psi_region3(r,omega)+0.5*Psih(r);
                else
                    psi_star_2D(r,omega)=Psih(r);
                end
            end
        end
    end
    
    %run('adjust_separatrix_psi_values')
    
    
    psi_star_2D_evol(f_rank,:,:)=psi_star_2D;
    
    % Create grids and convert polar coordinates to rectangularfigure(3)
    figure(3);
    set(gca,'FontSize',16);

    [THETA,RR] = meshgrid(scale_omega,scale_r);
    [XX,YY] = pol2cart(THETA,RR);
    
    
    %figure(3)
    %surfc(XX,YY,psi_star_2D,'edgecolor','none');view(0,90);
    contour(XX,YY,psi_star_2D,psi_contour_values);
    xlabel('r (m)');
    ylabel('r (m)');
    
    axis xy equal;
    max_psi_star=max(max(psi_star_2D));
    %imagesc(scale_r,scale_omega,psi_star_2D,[max_psi_star-0.05 max_psi_star])

    grid on;
    if (mod(f,2) == 0)
        display_psi_star_profile;
    end
    pause(0.22);
    
%     F=getframe;
%     [im,map] = frame2im(F);    %Return associated image data
%     
%     imwrite(im,filename,'bmp');
    
    disp('...done...');
    %pause
    rx_prev=rx;
    a1_coef_prev=min(a1_coef,3e2);
    a2_coef_prev=min(a2_coef,3e2);
    ksi0_prev=ksi0;
    
    volume1_prev=volume1;
    volume2_prev=volume2;
    volume3_prev=volume3;
end
%%
time1_offset=time_value

FRAME_LIMIT_SUP_FOR_CONTINUITY=f

% for next script
psi_star_2D_transition=psi_star_2D;
psi_star_2D_expulsion=psi_star_2D;


%close all
size_profile=ceil(xPsih_zero);
P_final_profile=P_initial_profile;
P_final_profile(1:size_profile)=interp1(pos_P_final_profile,P_final_profile_raw,1:size_profile,'pchip');


FILENAME=strcat(DATA_FOLDER,'psi_profiles_kadomtsev.mat')
save (FILENAME,'-append','plasma_beta_tot','P_initial_profile','P_final_profile','P0');

