
phi_rank
m_adjust_profile=zeros(Nradial,1);
m_guess=q_initial_profile(pTAE_inf-1)*nTAE;
% if (DOM_M<2)
%     m_adjust_profile1(pTAE_inf-1)=nTAE-mtheta1;
%     m_adjust_profile2(pTAE_inf-1)=nTAE-mtheta2-DOM_M;
% else
%     m_adjust_profile1(pTAE_inf-1)=nTAE-mtheta1+1;
%     m_adjust_profile2(pTAE_inf-1)=nTAE-mtheta2;
% end

% m_adjust_profile1(pTAE_inf-1)=-0.0+DOM_M;
% m_adjust_profile2(pTAE_inf-1)=-0.5;
% m_adjust_profile1(pTAE_inf-1)=-0.4+DOM_M;
% m_adjust_profile2(pTAE_inf-1)=1.5;
m_adjust_profile1(pTAE_inf-1)=0.4+DOM_M;
m_adjust_profile2(pTAE_inf-1)=0.4;


for r=pTAE_inf:pTAE_sup
    r_value=radial_r_value_flux(r);
    Rvalues=Rpos_PR_map(:,r);
%     vAvalues=vA_PR_map(:,r);
    vAvalues=vA_TAE;
    calculate_temporal_mode_structures;
    
    %     sb_width1=exp(-((r_value-r0)^2)/(2*MODE_WIDTH*(r3-r0)));
    %     sb_width2=exp(-((r_value-r3)^2)/(2*MODE_WIDTH*(r3-r0)));
    m_width1=coupling_m_mp1*exp(-((r_value-r1)^2)/(coupling_m_mp1*MODE_WIDTH*(r2-r1)));
    m_width2=exp(-((r_value-r2)^2)/(MODE_WIDTH*(r2-r1)));
    global_width=exp(-((r_value-rTAE)^2)/(TAE_WIDTH*(r2-r1)));
    
    Phi_PR_map(:,r)=Phi0*global_width*(m_width1*ksi_m_t+m_width2*ksi_mp1_t);
    iPhi_PR_map(:,r)=nTAE*Phi0*global_width*(m_width1*iksi_m_t+m_width2*iksi_mp1_t);
    dPhi_PR_map(:,r)=-(omega_TAE)*Phi0*global_width*(m_width1*iksi_m_t+OPP_W*m_width2*iksi_mp1_t);
    if UNIFORM_TAE_PROPAGATION==0
        A_PR_map(:,r)=(1/omega_values')*Phi0*global_width*(kpllm1*m_width1*ksi_m_t+OPP_W*kpllm2*m_width2*ksi_mp1_t);
        iA_PR_map(:,r)=(nTAE/omega_values')*Phi0*global_width*(kpllm1*m_width1*iksi_m_t+OPP_W*kpllm2*m_width2*iksi_mp1_t);
        dA_PR_map(:,r)=-Phi0*global_width*(kpllm1*m_width1*iksi_m_t+kpllm2*m_width2*iksi_mp1_t);
        omega_PR_map(:,r)=omega_values';
    else
        A_PR_map(:,r)=(1/omega_TAE)*Phi0*global_width*(kpllm1*m_width1*ksi_m_t+OPP_W*kpllm2*m_width2*ksi_mp1_t);
        iA_PR_map(:,r)=(nTAE/omega_TAE)*Phi0*global_width*(kpllm1*m_width1*iksi_m_t+OPP_W*kpllm2*m_width2*iksi_mp1_t);
        dA_PR_map(:,r)=-Phi0*global_width*(kpllm1*m_width1*iksi_m_t+kpllm2*m_width2*iksi_mp1_t);
    end
    
end

psi_star_PR_map(:,1:pTAE_sup)=-A_PR_map(:,1:pTAE_sup).*Rpos_PR_map(:,1:pTAE_sup).*Btor_PR_map(:,1:pTAE_sup)./Btot_PR_map(:,1:pTAE_sup);
