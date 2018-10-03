
%kpllm1=(nTAE-mtheta1./alphas_q_value)/R0;
%kpllm1=kTAE;
kpllm1=interp1(1:Nradial,kpllm1_profile,alphas_psi,'*cubic');
%m_adjust=interp1(1:Nradial,m_adjust_profile1,alphas_psi,'*cubic');
m_adjust=0.5;
M_theta1=(mtheta1+m_adjust).*alphas_theta-kpllm1.*alphas_s_value;

%kpllm2=(nTAE-mtheta2./alphas_q_value)/R0;
%kpllm2=kTAE;
kpllm2=interp1(1:Nradial,kpllm2_profile,alphas_psi,'*cubic');
%m_adjust=interp1(1:Nradial,m_adjust_profile2,alphas_psi,'*cubic');
m_adjust=0.5;
M_theta2=(mtheta2+m_adjust).*alphas_theta-kpllm2.*alphas_s_value;

%m_width1=coupling_m_mp1*exp(-((alphas_r_value-r1).^2)/(sqrt(coupling_m_mp1)*MODE_WIDTH*(r2-r1)));
%m_width2=exp(-((alphas_r_value-r2).^2)/(MODE_WIDTH*(r2-r1)));
%global_width=exp(-((alphas_r_value-rTAE).^2)/(TAE_WIDTH*(r2-r1)));
    m_width1=coupling_m_mp1*exp(-((alphas_r_value-r1).^2)/(MODE_WIDTH*(r2-r1))^2);
    m_width2=exp(-((alphas_r_value-r2).^2)/(MODE_WIDTH*(r2-r1))^2);
    global_width=exp(-((alphas_r_value-rTAE).^2)/(TAE_WIDTH*(r2-r1))^2);

ksi_m_t=cos(nTAE*alphas_pos_phi_wrapped-M_theta1-time*omega_TAE);
ksi_mp1_t=cos(nTAE*alphas_pos_phi_wrapped-M_theta2-OPP_W*time*omega_TAE);

iksi_m_t=-sin(nTAE*alphas_pos_phi_wrapped-M_theta1-time*omega_TAE);
iksi_mp1_t=-sin(nTAE*alphas_pos_phi_wrapped-M_theta2-OPP_W*time*omega_TAE);


ksi_m_t=(TAE_amplitude*Phi0*coupling_m_mp1)*global_width.*m_width1.*ksi_m_t;
ksi_mp1_t=(TAE_amplitude*Phi0)*global_width.*m_width2.*ksi_mp1_t;
iksi_m_t=(TAE_amplitude*Phi0*coupling_m_mp1)*global_width.*m_width1.*iksi_m_t;
iksi_mp1_t=(TAE_amplitude*Phi0)*global_width.*m_width2.*iksi_mp1_t;


alphas_iEpot=iksi_m_t+iksi_mp1_t;
alphas_iA=(1/omega_TAE)*(kpllm1.*iksi_m_t+kpllm2.*iksi_mp1_t);

alphas_Epot=ksi_m_t+ksi_mp1_t;
% alphas_grad_Phi=nTAE*(iksi_m_t+iksi_mp1_t);
% alphas_dt_Epot=-omega_TAE*(iksi_m_t+OPP_W*iksi_mp1_t);

alphas_Apll_value=(1/omega_TAE)*(kpllm1.*ksi_m_t+kpllm2.*ksi_mp1_t);

% the multiplication by -Rpos is not useful here
% alphas_grad_psi_star=(nTAE/omega_TAE)*(kpllm1.*iksi_m_t+kpllm2.*iksi_mp1_t);
% alphas_dt_psi_star=-(kpllm1.*iksi_m_t+OPP_W*kpllm2.*iksi_mp1_t);


alphas_Epot(OUTER_PART)=0;
alphas_Apll_value(OUTER_PART)=0;
alphas_iEpot(OUTER_PART)=0;
alphas_iA(OUTER_PART)=0;
alphas_grad_Phi(OUTER_PART)=0;
alphas_dt_Epot(OUTER_PART)=0;
alphas_grad_psi_star(OUTER_PART)=0;
alphas_dt_psi_star(OUTER_PART)=0;