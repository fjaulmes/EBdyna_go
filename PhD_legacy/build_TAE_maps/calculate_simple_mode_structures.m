

phi_nTAE=2*pi*(phi_rank-1)/(NB_PHI-1);

if UNIFORM_TAE_PROPAGATION==0
    omega_values=kTAE*vAvalues';
else
    omega_values=omega_TAE;
end

% vAmin=min(vAvalues);
% vAmax=max(vAvalues);
% % Rmin=min(Rvalues);
% % Rmax=max(Rvalues);
% mtheta=0*theta;
% % mtheta(find(Rvalues>=R0))=0.5+mtheta1-0.5*(Rvalues(find(Rvalues>=R0))-R0)/(Rmax-R0);
% % mtheta(find(Rvalues<R0))=0.5+mtheta1+0.5*(Rvalues(find(Rvalues<R0))-R0)/(Rmin-R0);
% mtheta(find(vAvalues>=vA_TAE))=0.5+mtheta1+0.5*(vAvalues(find(vAvalues>=vA_TAE))-vA_TAE)/(vAmax-vA_TAE);
% mtheta(find(vAvalues<vA_TAE))=0.5+mtheta1-0.5*(vAvalues(find(vAvalues<vA_TAE))-vA_TAE)/(vAmin-vA_TAE);
% 

s_coord=flc_s(:,r)';


% ksi_mm1_t=cos(-phi_nTAE+mtheta0*(theta+time*OMEGA_VA/(mtheta0)));
m_adjust=0.11;
eps_ksi=1;
options = optimset('TolFun',1e-3);
m_adjust=fminbnd(@(x) check_ksi_continuity(s_coord,vAvalues,omega_TAE,theta,qTAE,nTAE,x), -0.5,0.6,options);
M_theta=(qTAE*(nTAE)+m_adjust)*theta-(omega_TAE./vAvalues').*s_coord;
% while eps_ksi>0.005
%     M_theta=(qTAE*(nTAE)+m_adjust)*theta-(omega_TAE./vAvalues').*s_coord;
%     mod_M_theta=mod(M_theta,2*pi);
%     m_adjust=m_adjust+0.001;
%     eps_ksi=abs(mod_M_theta(1)-mod_M_theta(end));
% end
m_adjust_profile(r)=m_adjust;

ksi_m_t=cos(phi_nTAE-M_theta-time*omega_values);

% ksi_mm1_t=cos(-phi_nTAE+mtheta0*(theta+time*OMEGA_VA/(mtheta0)));
% ksi_mp1_t=cos(-phi_nTAE+mtheta2*(theta-time*OMEGA_VA/(mtheta2)));
% ksi_mp2_t=cos(-phi_nTAE+mtheta3*(theta-time*OMEGA_VA/(mtheta3)));

% iksi_mm1_t=nTAE*sin(-phi_nTAE+mtheta0*(theta+time*OMEGA_VA/(mtheta0)));
iksi_m_t=-nTAE*sin(phi_nTAE-M_theta-time*omega_values);
% iksi_mp1_t=nTAE*sin(-phi_nTAE+mtheta2*(theta-time*OMEGA_VA/(mtheta2)));
% iksi_mp2_t=nTAE*sin(-phi_nTAE+mtheta3*(theta-time*OMEGA_VA/(mtheta3)));
dksi_m_t=omega_values.*sin(phi_nTAE-M_theta-time*omega_values);

% Aksi_mm1_t=-sin(mtheta0*(theta-2*pi*time*OMEGA_VA/mtheta0));
% Aksi_m_t=-sin(mtheta1*(theta-2*pi*time*OMEGA_VA/mtheta1));
% Aksi_mp1_t=-sin(mtheta2*(theta+2*pi*time*OMEGA_VA/mtheta2));
% Aksi_mp2_t=-sin(mtheta3*(theta+2*pi*time*OMEGA_VA/mtheta3));