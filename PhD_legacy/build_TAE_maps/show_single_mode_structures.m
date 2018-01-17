close all
clear phi_val theta_val pot_s
pos_r1=interp1(radial_r_value_flux,1:Nradial,r1)
pos_r2=interp1(radial_r_value_flux,1:Nradial,r2)
pos_r3=interp1(radial_r_value_flux,1:Nradial,r3)

m_adjsut_profile=zeros(Nradial,1);
posTAE=round(interp1(q_initial_profile,1:Nradial,qTAE))
for p=1:NP
q_PR_map(p,:)=q_initial_profile;
end
flc_inc=sqrt((X_PR_map.^2+Z_PR_map.^2)./q_PR_map.^2+Rpos_PR_map.^2);

for f=1:NB_FRAME
    time=NB_OSCILLATIONS*((f-1)/(NB_FRAME))*2*pi/omega_TAE;
    for phi_rank=1:NB_PHI
        for r=2:pTAE_sup
            r_value=radial_r_value_flux(r);
            Rvalues=Rpos_PR_map(:,r);
            vAvalues=vA_PR_map(:,r);
            
            phi_nTAE=2*pi*(phi_rank-1)/(NB_PHI-1);
            
            if UNIFORM_TAE_PROPAGATION==0
                omega_values=kTAE*vAvalues';
            else
                omega_values=omega_TAE;
                vAvalues=vA_TAE;
            end
            
            vAmin=min(vAvalues);
            vAmax=max(vAvalues);
%             Rmin=min(Rvalues);
%             Rmax=max(Rvalues);
            mtheta=0*theta;
%             mtheta=qTAE*(nTAE+0.5-omega_TAE*(Rvalues./vAvalues));
%             mtheta(find(Rvalues>=R0))=0.5+mtheta1-0.5*(Rvalues(find(Rvalues>=R0))-R0)/(Rmax-R0);
%             mtheta(find(Rvalues<R0))=0.5+mtheta1+0.5*(Rvalues(find(Rvalues<R0))-R0)/(Rmin-R0);
%             mtheta(find(vAvalues>=vA_TAE))=0.5+mtheta1+0.5*(vAvalues(find(vAvalues>=vA_TAE))-vA_TAE)/(vAmax-vA_TAE);
%             mtheta(find(vAvalues<vA_TAE))=0.5+mtheta1-0.5*(vAvalues(find(vAvalues<vA_TAE))-vA_TAE)/(vAmin-vA_TAE);
%              mtheta=mtheta';
            flc_s=flc_inc(:,r)*DPHI*qTAE;
            flc_s(1)=0.5*flc_s(1);
            for sval=2:NP
                flc_s(sval)=flc_s(sval-1)+flc_inc(sval,r)*DPHI*qTAE;
            end
            s_coord=flc_s';
            
            
            % ksi_mm1_t=cos(-phi_nTAE+mtheta0*(theta+time*OMEGA_VA/(mtheta0)));
            options = optimset('TolFun',1e-6);
            kpllm=(nTAE-mtheta1/q_initial_profile(r))/R0;
%             m_adjust=fminsearch(@(x) check_ksi_continuity(s_coord,vA_TAE,omega_TAE,theta,qTAE,nTAE,x), (r-posTAE)/posTAE ,options);
            m_guess=q_initial_profile(r)*nTAE;
            m_adjust_guess=0.5*(m_guess-mtheta1);
            m_adjust=fminsearch(@(x) check_ksi_continuity(s_coord,kpllm,theta,mtheta1,m_adjust_guess,x),m_adjust_guess,options);
%             elseif r<posTAE
%                 m_adjust=fzero(@(x) check_ksi_continuity(s_coord,kpllm,theta,mtheta1,x), 0.3,0.5,options);
%             elseif r<pos_r2
%                 m_adjust=fzero(@(x) check_ksi_continuity(s_coord,kpllm,theta,mtheta1,x), 0.3,1.0,options);
%             elseif r<pos_r3
%                 m_adjust=fzero(@(x) check_ksi_continuity(s_coord,kpllm,theta,mtheta1,x), 0.7,1.3,options);
%             else
%                 m_adjust=fminbnd(@(x) check_ksi_continuity(s_coord,kpllm,theta,mtheta1,x), 1.2,2.2,options);
            m_adjust_profile(r)=m_adjust;
%             M_theta=(q_initial_profile(r)*(nTAE)+m_adjust)*theta-kpllm.*s_coord;
             M_theta=(mtheta1+m_adjust)*theta-kpllm.*s_coord;
           
            %             m_adjust=0.1;
%             eps_ksi=1;
%             while eps_ksi>0.005
%                 M_theta=(qTAE*(nTAE)+m_adjust)*theta-(omega_TAE./vAvalues').*flc_s';
%                 mod_M_theta=mod(M_theta,2*pi);              
%                 m_adjust=m_adjust+0.001;
%                 eps_ksi=abs(mod_M_theta(1)-mod_M_theta(end));
%             end
            ksi_m_t=cos(phi_nTAE-M_theta-time*omega_values);
            % ksi_mp1_t=cos(-phi_nTAE+mtheta2*(theta-time*OMEGA_VA/(mtheta2)));
            % ksi_mp2_t=cos(-phi_nTAE+mtheta3*(theta-time*OMEGA_VA/(mtheta3)));
            
            % iksi_mm1_t=nTAE*sin(-phi_nTAE+mtheta0*(theta+time*OMEGA_VA/(mtheta0)));
            iksi_m_t=-nTAE*sin(phi_nTAE-mtheta.*theta-time*omega_values);
            % iksi_mp1_t=nTAE*sin(-phi_nTAE+mtheta2*(theta-time*OMEGA_VA/(mtheta2)));
            % iksi_mp2_t=nTAE*sin(-phi_nTAE+mtheta3*(theta-time*OMEGA_VA/(mtheta3)));
            
            % Aksi_mm1_t=-sin(mtheta0*(theta-2*pi*time*OMEGA_VA/mtheta0));
            % Aksi_m_t=-sin(mtheta1*(theta-2*pi*time*OMEGA_VA/mtheta1));
            % Aksi_mp1_t=-sin(mtheta2*(theta+2*pi*time*OMEGA_VA/mtheta2));
            % Aksi_mp2_t=-sin(mtheta3*(theta+2*pi*time*OMEGA_VA/mtheta3));
            global_width=exp(-((r_value-rTAE)^2)/(2*TAE_WIDTH*(r2-r1)));
            
            Phi_PR_map(:,r)=Phi0*global_width*(ksi_m_t);
            iPhi_PR_map(:,r)=Phi0*global_width*(iksi_m_t);
            
        end
        Epot_map_phi(phi_rank,:,:)=Phi_PR_map(:,pTAE_inf:pTAE_sup);
    end
    posTAE_value=posTAE
    Phi_theta_phi_map = squeeze(Epot_map_phi(1:end-1,:,posTAE_value-pTAE_inf+1));
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
        theta_val(s)=mod(phi_val(s)*(1/qTAE),2*pi);
        pot_s(s)=interp1(theta,Phi_theta_phi_map(s,:),theta_val(s));
    end
    phi_left=round((qTAE-1)*(size_phi-1));
    for s=1:phi_left
        phi_val(s+size_phi)=4*pi+4*pi*(s-1)/(size_phi+1);
        theta_val(s+size_phi)=mod(phi_val(s+size_phi)*(1/qTAE),2*pi);
        pot_s(s+size_phi)=interp1(theta,Phi_theta_phi_map(s,:),theta_val(s+size_phi));
    end
%     phi_left=round((qTAE-1)*(size_phi-1))
%     for s=1:phi_left
%         phi_val(s+2*size_phi)=4*pi+2*pi*(s-1)/(size_phi+1);
%         theta_val(s+2*size_phi)=mod(phi_val(s+2*size_phi)*(1/qTAE),2*pi);
%         pot_s(s+2*size_phi)=interp1(theta,Phi_theta_phi_map(s,:),theta_val(s+2*size_phi));
%     end
    plot(phi_val,pot_s)
%%
    title(num2str(f))
    pause(0.2)
end