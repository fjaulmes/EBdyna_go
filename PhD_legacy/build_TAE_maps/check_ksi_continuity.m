function [ eps_ksi ] = check_ksi_continuity(s_coord,kpll,theta,mtheta, m_adjust )
% for optimization of poloidal profile according to parallel wave number

M_theta=(mtheta+m_adjust)*theta-kpll.*s_coord;
% M_theta=(qTAE*(nTAE)+m_adjust)*theta;
mod_M_theta=mod(M_theta,2*pi);
% if abs(m_adjust-m_guess)<0.5
    
% eps_ksi=(mod_M_theta(1)-mod_M_theta(end))^2+(cos(M_theta(1))-cos(M_theta(end)))^2;
eps_ksi=(cos(M_theta(1))-cos(M_theta(end)))^2;
% else
%     eps_ksi=1000;
% end

end

