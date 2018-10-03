

Nvalues=zeros(length(radial_bins),1);
Tvalues=zeros(length(radial_bins),1);
Pvalues=zeros(length(radial_bins),1);

SKIP_EDGE=1

for psi_pos=psi_pos_inf:psi_pos_sup
    
    disp('-------------------------------------------------')
    Nvalues(psi_pos)=mean(N_profile_radial_ini(radial_bins(psi_pos)-radial_bin_half_size:radial_bins(psi_pos)+radial_bin_half_size));
    
    Ti_psi=mean(Ti_profile(radial_bins(psi_pos)-radial_bin_half_size+SKIP_EDGE:radial_bins(psi_pos)+radial_bin_half_size+SKIP_EDGE)) %eV
% 	vthe=sqrt(2*eV*Te_psi/me);
    Pvalues(psi_pos)=Ti_psi*Nvalues(psi_pos);
    
    % lowering edge values to keep reasonable amount confined
%     if psi_pos>=(psi_pos_sup-30)
%         Ti_psi=Ti_psi*(0.999-(0.75-0.025*(psi_pos_sup-psi_pos))^4.5);
%     end
%     Tvalues(psi_pos)=Ti_psi;
    if psi_pos>=(psi_pos_sup-20)
        Ti_psi=Ti_psi*(0.9999-(0.4-0.02*(psi_pos_sup-psi_pos))^2.8);
    end
    Tvalues(psi_pos)=Ti_psi;

        
    vthi=sqrt(2*eV*Ti_psi/mHe);

    
    for n=1:N_energy_bins+1
        f_D_speed_dist(psi_pos,n)=((mHe/(2*pi*Ti_psi*eV))^(3/2))*(4*pi*speed_D_range(n)^2)*exp(-(mHe*speed_D_range(n)^2)/(2*Ti_psi*eV));
        f_D_E_dist(psi_pos,n)=eV*energy_bin_size*(2*sqrt(energy_D_range(n)*eV/pi).*(1/(Ti_psi*eV))^(3/2)).*exp(-energy_D_range(n)/Ti_psi);
    end
    
end
Nvalues=Pvalues./Tvalues;
% for psi_pos=psi_pos_sup-11:psi_pos_sup
%     Nvalues(psi_pos)=(1+0.001*(psi_pos-psi_pos_sup+11))*Nvalues(psi_pos);
% end
% Nvalues(end)=1.00025*Nvalues(end-1);

% for psi_pos=psi_pos_inf:psi_pos_sup
%     Ti_psi=Tvalues(psi_pos);
%     % lowering edge values to keep reasonable amount confined
%     if psi_pos>=(psi_pos_sup-5)
%         Ti_psi=Ti_psi*(0.99-(0.75-0.15*(psi_pos_sup-psi_pos))^2);
%     end
%     Tvalues(psi_pos)=Ti_psi;
% end
% 
% EDGE_CORR=1.05;
% SIZE_EDGE_CORR=6
% edge_corr_values=linspace(1,EDGE_CORR,SIZE_EDGE_CORR).^4;
% % edge_corr_values=edge_corr_values*EDGE_CORR/edge_corr_values(end);
% 
% Nvalues(end-SIZE_EDGE_CORR+1:end)=Nvalues(end-SIZE_EDGE_CORR+1:end).*edge_corr_values';
%%
figure;
hold on
% plot(C_density_profile)
plot(1:Nradial,Ti_profile);
plot(radial_bins,Tvalues)


% Nvalues=Pvalues./Tvalues;

%%
figure;
hold on
% plot(C_density_profile)
plot(1:Nradial,N_profile_radial_ini);
plot(radial_bins,Nvalues)

%%
figure;
hold on
% plot(C_density_profile)
plot(1:Nradial,N_profile_radial_ini.*Ti_profile);
plot(radial_bins,Nvalues.*Tvalues)