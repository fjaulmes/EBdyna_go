

% *************************************************************
% calculate volume of each flux surface
% *************************************************************
JPHI_XZ_map(find(isnan(JPHI_XZ_map)))=0;

tor_flux_profile=zeros(1,Nradial);
surf_flux=zeros(1,Nradial);
volume_flux=zeros(1,Nradial);
I_flux_map=Z_psi_fit_up*0;
I_flux=zeros(1,Nradial);
I_flat_flux=zeros(1,Nradial);
radial_r_value_flux=zeros(1,Nradial);


%DS=(DX)^2;

for (n=2:Nradial)
    %I_flux(n)=0;
    for(x=X1_Nradial(n):X2_Nradial(n))
        %Z_max=round(0.5*(Z_psi(n,x)+Z_psi_fit(n,x)))+mid_Z+1;
        Z_max_up=Z_psi_fit_up(n,x);
        Z_max_down=Z_psi_fit_down(n,x);
        DS=(Z_max_up-Z_max_down)*DX;
        surf_flux(n)=surf_flux(n)+DS;
        for z=floor(Z_max_down/DX+mid_Z):ceil(Z_max_up/DX+mid_Z)
        tor_flux_profile(n)=tor_flux_profile(n)+Bphi_XZ_map(x,z)*DX*DX;
        end
        %surf_flux(n)=surf_flux(n)+DS;
    end
end

scale_tor_x=[X1_Nradial(n) X2_Nradial(n)];
size_tor_X=X2_Nradial(n)-X1_Nradial(n)+1;
volume_tor=zeros(Nradial,size_tor_X);
%volume_tor_diff=zeros(Nradial,size_tor_X);

for (n=2:Nradial)
    for(x=X1_Nradial(n):X2_Nradial(n))
        Z_max_up=Z_psi_fit_up(n,x);
        Z_max_down=Z_psi_fit_down(n,x);
        DS=(Z_max_up-Z_max_down)*DX;
        volume_tor(n,x)=2*pi*Rpos(x)*DS;
    end
end

for (n=2:Nradial)
    volume_tor_diff(n,:)=volume_tor(n,:)-volume_tor(n-1,:);
end



if psi_global<0
surf_flux_pol=polyfit(psi_profile+psi_global,surf_flux,8);
surf_flux_fit=polyval(surf_flux_pol,psi_profile+psi_global);
else
surf_flux_pol=polyfit(psi_scale-psi_global,surf_flux,8);
surf_flux_fit=polyval(surf_flux_pol,psi_scale-psi_global);
end
if (surf_flux_fit(2)<surf_flux(2))
    surf_flux_fit=surf_flux_fit+(surf_flux(2)-surf_flux_fit(2));
end
surf_flux_fit(1)=surf_flux(1);
surf_flux_fit(2)=0.5*(surf_flux_fit(3)+surf_flux_fit(1));

surf_flux=surf_flux_fit;




disp('Volume of vaccuum vessel = ');
disp(sum(volume_tor(end,:)));





% % adjust toroidal flux calculation
% for (r=2:Nradial)
%     tor_flux_profile(r)=tor_flux_profile(r)/surf_flux(r);
% end
