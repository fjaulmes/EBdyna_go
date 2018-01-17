% reset_data_analysis_environment

filename=strcat(DATA_FOLDER,'plasma_current.mat');
load(filename);
filename=strcat(DATA_FOLDER,'volume_flux_geometry.mat');
load(filename);
filename=strcat(DATA_FOLDER,'B_fields.mat');
load(filename,'BpolX_initial_map');
load(filename,'BpolZ_initial_map');

Bpol_XZ_map=sqrt(BpolX_initial_map.^2+BpolZ_initial_map.^2);

mag_energ_pol=zeros(1,Nradial);
LI_profile=zeros(1,Nradial);
li_profile=zeros(1,Nradial);




for (n=2:Nradial)
    for(x=X1_Nradial(n):X2_Nradial(n))
        Z_max_up=Z_psi_fit_up(n,x);
        Z_max_down=Z_psi_fit_down(n,x);
        DS=(Z_max_up-Z_max_down)*DX;
        for z=floor(Z_max_down/DX+mid_Z):ceil(Z_max_up/DX+mid_Z)
            mag_energ_pol(n)=mag_energ_pol(n)+Bpol_XZ_map(x,z)^2*DX*DX*2*pi*(X_scale(x)+R0);
        end
    end
end

mag_energ_pol=mag_energ_pol/(2*mu0);

LI_profile=2*mag_energ_pol./I_flux.^2;

li_profile=2*LI_profile/(mu0*R0);

n=psi_rank_q1
figure(1)
hold on
grid on
plot(X1_Nradial(n):X2_Nradial(n),Z_psi_fit_up(n,X1_Nradial(n):X2_Nradial(n)))
plot(X1_Nradial(n):X2_Nradial(n),Z_psi_fit_down(n,X1_Nradial(n):X2_Nradial(n)))

a1=0.5*(X_scale(X2_Nradial(n))-X_scale(X1_Nradial(n)))
b1=0.5*(max(Z_psi_fit_up(n,:))-min(Z_psi_fit_down(n,:)))

% modified elongation correction from Porcelli
kappa1=(b1/a1)^2


filename=strcat(DATA_FOLDER,'plasma_current.mat');
save(filename,'-append','mag_energ_pol','li_profile');

filename=strcat(DATA_FOLDER,'volume_flux_geometry.mat');
save(filename,'-append','a1','b1','kappa1');