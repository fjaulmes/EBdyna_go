
% close all



r_nu_cont13=rx;
omega_cont13=Nomega;
rPsih_zero=interp1((1:NR),scale_r,xPsih_zero);

f_initial=round(xPsih_zero-xPsih_max)-1;



f_final=297;
f_size=f_final-f_initial+1;
Df=(0.97-time0_offset)/(f_final-f_initial+2);

psi_star_2D_evol_final=zeros(f_size,NR,Nomega);

rmix=rPsih_zero;
rx=rmix;
ksi0=-rmix;

pos_rmix=xPsih_zero;
pos_psi_rx=pos_rmix;
x_nu_cont13=pos_psi_rx;
x_nu_cont13_initial=2;
r_nu_cont13=rx;
r_nu_cont13_initial=2*Dr_avg(round(pos_psi_rx));

area_region3=interp1((1:NR),surf_flux_precise(1:NR),pos_rmix,'pchip');


%alpha=atan(2);
for (f=f_initial:f_final)
    
    %%
    psi_star_2D=zeros(NR,Nomega);

    f_rank=f+1;
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
    alpha=(Df*(f-f_initial+1)*(1-atan(2)/pi)+atan(2)/pi)*pi;
%     filename='reconnection_movie\t0';
%     filename=strcat(filename,frame_name,'.bmp');
%     disp(filename);
    
    time_value=((1-time1_offset)/(pi-atan(2)))*(alpha-atan(2))+time1_offset;
    time_evol(f_rank)=time_value;
    disp('time=');
    disp(time_evol(f_rank));
    
    rx_evol(f_rank)=rx;
    ksi0_evol(f_rank)=abs(ksi0);
    disp('rx=');
    disp(rx);
    
    pos_psi_rx_round=ceil(pos_psi_rx);

    %psi_rx=interp1(scale_X(1:xPsih_zero),Psih(1:xPsih_zero),pos_psi_rx);
    psi_rx=0;
    psi_rank_rx=interp1(Psih_values,[1:NPsih],psi_rx);
    x_psi_final_rx=xPsih_zero;
    
    %diff_in_psi_rx=abs(Psih(round((rx-ksi0)/DX))-Psih(round(rx/DX)));
    
    
    %psi value at limit position rx
    %psi_limit_23=psi_region1(pos_psi_rx,half_Nomega);
    
    min_psi_region1=min(min(psi_region1));
    max_psi_region1=max(max(psi_region1));
    %Dpsi_region1=max_psi_region1-min_psi_region1;
    
    Psih_limit13=0;
    Psih_limit23=0;
    
     
    % Deriving the expression of flux contours in the bean shaped zone

    a1_coef=a1_coef_prev+0.05*(1-a1_coef_prev);
    a2_coef=a2_coef_prev+0.05*(1-a2_coef_prev);
    calculate_nu_map;
    
    nu_min_rx=interp1(1:NR,nu(:,omega_cont13),pos_psi_rx);
    nu_min=nu_min_rx;
    nu=nu-nu_min;
    nu_max=max(nu(1:pos_psi_rx_round,1));
    nu=nu_max-nu;
    nu_max_rx=interp1(1:NR,nu(:,1),pos_psi_rx);
    nu=(nu/nu_max_rx)*area_region3;
    nu_max=interp1(1:NR,nu(:,omega_cont13),pos_psi_rx);
    
    N_NU=round(0.25*NR);
    Dnu=(nu_max)/(N_NU-1);
    nu_values=((1:N_NU)-1)*Dnu;

   


%     if (f<FRAME_LIMIT_SUP_FOR_CONTINUITY)
%         calculate_nu_surfaces;
%         calculate_nu_map_optimized;
%         calculate_nu_map;
%         
%         nu_min_rx=interp1(1:NR,nu(:,omega_cont13),pos_psi_rx);
%         nu_min=nu_min_rx;
%         nu=nu-nu_min;
%         nu_max=max(nu(1:pos_psi_rx_round,1));
%         nu=nu_max-nu;
%         nu_max_rx=interp1(1:NR,nu(:,1),pos_psi_rx);
%         nu=(nu/nu_max_rx)*area_region3;
%         nu_max=interp1(1:NR,nu(:,omega_cont13),pos_psi_rx);
%         
%         %determining N_NU arbitrary steps
%         %for scaling the flux in nu scale
%         N_NU=500;
%         Dnu=(nu_max)/(N_NU-1);
%         
%         for (n=1:N_NU)
%             nu_values(n)=(n-1)*Dnu;
%         end
%     else
        a1_coef
        a2_coef
%     end
    
    
    %-----------------------------------------
    % Helical flux in the magnetic island
    %-----------------------------------------    
    min_psi3=0;
    max_psi3=psi_star_max;
    
    calculate_nu_surfaces;
    
    %psi3_nu=interp1(scale_X(1:round(x_psi_final_rx)+1),Psih_final(1:round(x_psi_final_rx)+1),nu_values);
    %psi3_nu=interp1(scale_X(1:x_psi_final_rx),Psih_final(1:x_psi_final_rx),nu_values);

    for (r=1:NR)
        %r_value=r*Dr;
        %for (omega=1:Nomega)
            %nu_val=nu(r,omega);
        psi_region3(r,:)=interp1(nu_values,psi3_nu,nu(r,:),'pchip');
        %end
    end
    psi_region3(isnan(psi_region3))=0;    

    % figure(5);hold on;
    % plot(psi3);
    
    % total_psi3_flux(f_rank)=0;
    % total_psi_flux(f_rank)=0;
    
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
    
%     figure(2);
%     
%     hold on;
%     plot(scale_r(1:pos_psi_rx+2),psi_star_2D(1:pos_psi_rx+2,1));
    
    % ignoring the transition frames with the push on the other side
    % this looks not very physical
%     if (f<=f_initial+2)
%         psi_star_2D_transition=psi_star_2D;
%     elseif (f==f_initial+3)
%         psi_star_2D_evol_final(f-f_initial-2,:,:)=0.25*psi_star_2D+0.75*psi_star_2D_expulsion;
%         psi_star_2D_evol_final(f-f_initial-1,:,:)=0.5*(psi_star_2D+psi_star_2D_expulsion);
%         psi_star_2D_evol_final(f-f_initial,:,:)=0.75*psi_star_2D+0.25*psi_star_2D_expulsion;
%         psi_star_2D_evol_final(f-f_initial+1,:,:)=psi_star_2D;
%     elseif (f==f_initial+3)
%         psi_star_2D_evol_final(f-f_initial,:,:)=0.5*(psi_star_2D+psi_star_2D_evol_final(f-f_initial,:,:));
%         psi_star_2D_evol_final(f-f_initial+1,:,:)=psi_star_2D;
%     end
    if (f==f_initial)
%         psi_star_2D_transition=0.5*(psi_star_2D+psi_star_2D_expulsion);
%         psi_star_2D_evol_final(f-f_initial+1,:,:)=psi_star_2D_transition;
        psi_star_2D_evol_final(f-f_initial+1,:,:)=psi_star_2D;
    elseif (f==f_initial+1)
%         psi_star_2D_evol_final(f-f_initial,:,:)=0.5*(psi_star_2D+psi_star_2D_transition);
        psi_star_2D_evol_final(f-f_initial+1,:,:)=psi_star_2D;
    else
        psi_star_2D_evol_final(f-f_initial+1,:,:)=psi_star_2D;
    end
    
    % Create grids and convert polar coordinates to rectangular
    [THETA,RR] = meshgrid(scale_omega,scale_r);
    [XX,YY] = pol2cart(THETA,RR);
    
    
    figure(3)
    %surfc(XX,YY,psi_star_2D,'edgecolor','none');view(0,90);
    contour(XX,YY,psi_star_2D,psi_contour_values);
    
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
    a1_coef_prev=a1_coef;
    a2_coef_prev=a2_coef;
    
end

psi_star_2D_evol(f_initial+1:f_final+1,:,:)=psi_star_2D_evol_final(:,:,:);
rx_evol(f_initial+1:f_final+1)=rmix;
ksi0_evol(f_initial+1:f_final+1)=rmix;
%save psi_star_evol.mat psi_star_2D_evol 


%%
figure(8);grid on;hold on;

set(gca,'FontSize',18);

%P_final_profile(end+1)=P_final_profile(end);
P_final_profile_result=P_initial_profile;

P_final_profile_result(1:round(xPsih_zero/PRECISE_MESH))=interp1(scale_r,P_final_profile(1:length(scale_r)),radial_r_value(1:round(xPsih_zero/PRECISE_MESH)),'pchip');

plot(radial_r_value, P_initial_profile,'b--');
plot(radial_r_value, P_final_profile_result,'r');

P_final_profile=P_final_profile_result;

set( findobj(gcf,'Type','line'),'LineWidth',1.5);
ylabel('pressure (Pa)')
xlabel('r (m)')
xlim([0 1]);

legend('before collapse','afer collapse');

