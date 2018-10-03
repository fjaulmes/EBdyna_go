Z_psi_fit_up=Z_psi_fit;

% *************************************************************
% calculate volume of each flux surface
% *************************************************************

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
        %I_flux_map(n,x)=sum(JPHI_XZ_map(x,round(Z_max_down/DX)+mid_Z:round(Z_max_up/DX)+mid_Z))*DX*DX;
        %I_flat_flux(n)=I_flat_flux(n)+JPHI_flat_pressure_XZ_map(x,mid_Z)*DS;
        % toroidal flux calculation
        for z=floor(Z_max_down/DX+mid_Z):ceil(Z_max_up/DX+mid_Z)
            tor_flux_profile(n)=tor_flux_profile(n)+Bphi_XZ_map(x,z)*DX*DX;
        end
        %surf_flux(n)=surf_flux(n)+DS;
    end
end
for (n=2:Nradial)
    I_flux(n)=sum(I_flux_map(n,:));
end
scale_tor_x=[min(X1_Nradial) max(X2_Nradial)];
size_tor_X=max(X2_Nradial)-min(X1_Nradial)+1;
volume_tor=zeros(Nradial,size_tor_X);
volume_tor_diff=volume_tor*0;
%volume_tor_diff=zeros(Nradial,size_tor_X);

figure
hold on
for (n=2:Nradial)
    for(x=min(X1_Nradial):max(X2_Nradial))
        Z_max_up=Z_psi_fit_up(n,x);
        Z_max_down=Z_psi_fit_down(n,x);
        DS=(Z_max_up-Z_max_down)*DX;
        volume_tor(n,x-min(X1_Nradial)+1)=2*pi*Rpos(x)*DS;
    end
    plot(X1_Nradial(n):X2_Nradial(n),Z_psi_fit_up(n,X1_Nradial(n):X2_Nradial(n)));
    plot(X1_Nradial(n):X2_Nradial(n),Z_psi_fit_down(n,X1_Nradial(n):X2_Nradial(n)));
    
end


for (n=2:Nradial)
    volume_tor_diff(n,:)=volume_tor(n,:)-volume_tor(n-1,:);
end


if SIGN_CO_CURRENT_FIELD>0
    surf_flux_pol=polyfit(psi_profile+psi_global,surf_flux,4);
    surf_flux_fit=polyval(surf_flux_pol,psi_profile+psi_global);
else
    surf_flux_pol=polyfit(psi_profile-psi_global,surf_flux,4);
    surf_flux_fit=polyval(surf_flux_pol,psi_profile-psi_global);
end

%matching with the original first value

surf_flux_fit=surf_flux_fit-surf_flux_fit(2)+surf_flux(2);

surf_flux=surf_flux_fit-surf_flux_fit(1);
%surf_flux(1)=0;

radial_r_value_flux=sqrt(surf_flux/pi);

volume_flux=2*pi*surf_flux.*(R0+xi_psi_Nradial);



disp('Volume of vaccuum vessel = ');
disp(sum(volume_tor(end,:)));

