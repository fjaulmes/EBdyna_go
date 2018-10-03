
s_coord=flc_s(:,r)';
omega_values=kTAE*vAvalues';


% mtheta0t=0*theta+mtheta0;
% mtheta1t=0*theta+mtheta1;
% mtheta2t=0*theta+mtheta2;
% mtheta3t=0*theta+mtheta3;

% s_coord=flc_s(:,r)';
% 
% m_adjust=0.11;
% options = optimset('TolFun',1e-3);
% m_adjust=fminbnd(@(x) check_ksi_continuity(s_coord,vAvalues,omega_TAE,theta,qTAE,nTAE,x), -0.5,0.6,options);
% M_theta=(qTAE*(nTAE)+m_adjust)*theta-(omega_TAE./vAvalues').*s_coord;

m_guess=q_initial_profile(r)*nTAE;
% m_adjust_guess=0.8*(m_guess-mtheta1);
% kpllm1=(omega_TAE./vAvalues');
%kpllm1_coef=min(max((1-(radial_r_value_flux(r)-r1)/(rTAE-r1)),0),1);
%kpllm1=(nTAE-mtheta1/q_initial_profile(r))/R0;
kpll_TAE=2*pi/(2*length_fl_profile(r));
kpllm1=OPP_W*kpll_TAE;
% kpllm1=kTAE;

%kpllm1=(kpllm1_coef)*kpllm1+(1-kpllm1_coef)*kpll_TAE;
% kpllm1=kpll_TAE;
kpllm1_profile(r)=kpllm1;
% m_adjust=(m_guess-mtheta1);
% m_adjust=fminbnd(@(x) check_ksi_continuity(s_coord,kpllm1,theta,mtheta1,x),m_adjust_profile1(r-1)-0.3,m_adjust_profile1(r-1)+0.3,options);
% if abs(m_adjust_guess>0.1)
% m_adjust=fminsearch(@(x) check_ksi_continuity(s_coord,kpllm1,theta,mtheta1,m_adjust_guess,x),m_adjust_profile1(r-1),options);
% else
% m_adjust=(m_guess-mtheta1);
% end
m_adjust=0.5;
m_adjust_profile1(r)=m_adjust;
%M_theta1=(mtheta1+m_adjust)*theta+kpllm1.*s_coord;
M_theta1=(mtheta1)*theta;

% m_adjust_guess=0.8*(m_guess-mtheta2);
% kpllm2=(omega_TAE./vAvalues');
%kpllm2_coef=min(max(((radial_r_value_flux(r)-rTAE)/(r2-rTAE)),0),1);
% kpllm2=kTAE;
kpllm2=kpll_TAE;

%kpllm2=(nTAE-mtheta2/q_initial_profile(r))/R0;
%kpll_TAE=2*pi/(2*length_fl_profile(r));
%kpllm2=(kpllm2_coef)*kpllm2+(1-kpllm2_coef)*kpll_TAE;
kpllm2_profile(r)=kpllm2;

% m_adjust=(m_guess-mtheta2);
% if abs(m_adjust_guess>0.1)
% m_adjust=fminbnd(@(x) check_ksi_continuity(s_coord,kpllm2,theta,mtheta2,x),m_adjust_profile2(r-1)-0.3,m_adjust_profile2(r-1)+0.3,options);
% else
% m_adjust=(m_guess-mtheta2);
% end
m_adjust=-0.5;
m_adjust_profile2(r)=m_adjust;
%M_theta2=(mtheta2+m_adjust)*theta+kpllm2.*s_coord;
M_theta2=(mtheta2)*theta;


% ksi_mm1_t=cos(phi_nTAE-mtheta0t.*theta-time*omega_values);
ksi_m_t=cos(M_theta1-phi_nTAE+time*omega_values);
ksi_mp1_t=cos(M_theta2-phi_nTAE+time*omega_values);
% ksi_mp2_t=cos(phi_nTAE-mtheta3t.*theta+time*omega_values);

% iksi_mm1_t=-nTAE*sin(phi_nTAE-mtheta0t.*theta-time*OMEGA_VA);
iksi_m_t=-sin(M_theta1-phi_nTAE+time*omega_values);
iksi_mp1_t=-sin(M_theta2-phi_nTAE+time*omega_values);
% iksi_mp2_t=-nTAE*sin(phi_nTAE-mtheta3t.*theta+time*OMEGA_VA);

% Aksi_mm1_t=-sin(mtheta0*(theta-2*pi*time*OMEGA_VA/mtheta0));
% Aksi_m_t=-sin(mtheta1*(theta-2*pi*time*OMEGA_VA/mtheta1));
% Aksi_mp1_t=-sin(mtheta2*(theta+2*pi*time*OMEGA_VA/mtheta2));
% Aksi_mp2_t=-sin(mtheta3*(theta+2*pi*time*OMEGA_VA/mtheta3));