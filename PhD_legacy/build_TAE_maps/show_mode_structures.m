figure(3)

% close all
clear phi_val theta_val pot_s
pos_r1=interp1(radial_r_value_flux,1:Nradial,r1)
pos_r2=interp1(radial_r_value_flux,1:Nradial,r2)
% pos_r3=interp1(radial_r_value_flux,1:Nradial,r3)


for offset=1:2:size_r_TAE-1
    posTAE_value=pTAE_inf+offset
    qcheck=q_initial_profile(posTAE_value)
    Phi_theta_phi_map= squeeze(Epot_map_phi(1:end-1,:,posTAE_value-pTAE_inf+1));
%     Phi_theta_phi_map = squeeze(Epot_map_phi(1:end-1,:,posTAE_value-pTAE_inf+1));
    for n=2:nTAE
        Phi_theta_phi_map=[Phi_theta_phi_map ; squeeze(Epot_map_phi(1:end-1,:,posTAE_value-pTAE_inf+1)) ];
    end
    for n=1:nTAE
        Phi_theta_phi_map=[Phi_theta_phi_map ; squeeze(Epot_map_phi(1:end-1,:,posTAE_value-pTAE_inf+1)) ];
    end
    size_phi=size(Phi_theta_phi_map,1)-1;
    phi_scale=((0:size_phi)/size_phi)*4*pi;
%     imagesc(phi_scale,theta,Phi_theta_phi_map');axis xy;hold on
%     plot(phi_val,theta_val,'g.','linewidth',4)
%%

    for s=1:size_phi
        phi_val(s)=4*pi*(s-1)/(size_phi+1);
        theta_val(s)=mod(phi_val(s)*(1/qcheck),2*pi);
        pot_s(s)=interp1(theta,Phi_theta_phi_map(s,:),theta_val(s));
    end
    phi_left=round((qTAE-1)*(size_phi-1));
    for s=1:phi_left
        phi_val(s+size_phi)=4*pi+4*pi*(s-1)/(size_phi+1);
        theta_val(s+size_phi)=mod(phi_val(s+size_phi)*(1/qcheck),2*pi);
        pot_s(s+size_phi)=interp1(theta,Phi_theta_phi_map(s,:),theta_val(s+size_phi));
    end
%     phi_left=round((qTAE-1)*(size_phi-1))
%     for s=1:phi_left
%         phi_val(s+2*size_phi)=4*pi+2*pi*(s-1)/(size_phi+1);
%         theta_val(s+2*size_phi)=mod(phi_val(s+2*size_phi)*(1/qTAE),2*pi);
%         pot_s(s+2*size_phi)=interp1(theta,Phi_theta_phi_map(s,:),theta_val(s+2*size_phi));
%     end
    plot(phi_val,pot_s,'color',[1 1 1]*offset/size_r_TAE);
hold on
%%
    title(num2str(f))
    pause
end
close all