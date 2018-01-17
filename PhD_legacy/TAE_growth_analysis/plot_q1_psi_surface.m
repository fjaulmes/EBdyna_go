psi_q1=interp1(q_initial_profile,psi_scale,1);
eps=1-q_initial_profile(psi_rank_q1);
sign_dist=sign(eps);
if (sign_dist ~= 0)
    dist_psi=abs(eps/(q_initial_profile(psi_rank_q1)-q_initial_profile(psi_rank_q1+sign_dist)));
else
    dist_psi=0;
end

Z_psi_q1=Z_psi_fit_up(psi_rank_q1,:)*(1-dist_psi)+Z_psi_fit_up(psi_rank_q1+sign_dist,:)*(dist_psi);
Z_psi_q1_down=Z_psi_fit_down(psi_rank_q1,:)*(1-dist_psi)+Z_psi_fit_down(psi_rank_q1+sign_dist,:)*(dist_psi);
% xi_value_q1=xi_psi_fit(psi_rank_q1)*(1-dist_psi)+xi_psi_fit(psi_rank_q1+sign_dist)*dist_psi;
X1_q1=X1_Nradial(psi_rank_q1)*(1-dist_psi)+X1_Nradial(psi_rank_q1+sign_dist)*dist_psi;
X1_q1_rank=round(X1_q1);
X2_q1=X2_Nradial(psi_rank_q1)*(1-dist_psi)+X2_Nradial(psi_rank_q1+sign_dist)*dist_psi;
X2_q1_rank=round(X2_q1);
Z_psi_q1(1:X1_q1_rank)=Z_axis;
Z_psi_q1(X2_q1_rank:Nradial)=Z_axis;
Z_psi_q1_down(1:X1_q1_rank)=Z_axis;
Z_psi_q1_down(X2_q1_rank:Nradial)=Z_axis;

r_value_q1_mean=radial_r_value_flux(psi_rank_q1)

hold on;
grid on

plot(X_scale(X1_q1_rank:X2_q1_rank),Z_psi_q1(X1_q1_rank:X2_q1_rank),'k--','linewidth',3);
plot(X_scale(X1_q1_rank:X2_q1_rank),Z_psi_q1_down(X1_q1_rank:X2_q1_rank),'k--','linewidth',3);