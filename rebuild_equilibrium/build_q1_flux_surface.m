


% *************************************************************
% defining the safety factor map
% *************************************************************

disp('... Shape of q=1 flux surface ...');


q_map=q_XZ_map;

Z_psi_q1=zeros(1,NR)+Z_axis;
Z_psi_q1_down=zeros(1,NR)+Z_axis;


% *************************************************************
% position of q=1 surface
% *************************************************************



q1=ones(1,Nradial);
[eps Nq1]=min(abs(q1-q_profile));
psi_rank_q1=Nq1;
% of course q=1 should not be the separatrix!
if(psi_rank_q1>=Nradial)
	psi_rank_q1=Nradial-1;
end
eps=1-q_profile(psi_rank_q1);
sign_dist=sign(eps);
if (sign_dist ~= 0)
    dist_psi=abs(eps/(q_profile(psi_rank_q1)-q_profile(psi_rank_q1+sign_dist)));
else
    dist_psi=0;
end
psi_q1=psi_profile(psi_rank_q1)*(1-dist_psi)+psi_profile(psi_rank_q1+sign_dist)*dist_psi;


Z_psi_q1=Z_psi_fit(psi_rank_q1,:)*(1-dist_psi)+Z_psi_fit(psi_rank_q1+sign_dist,:)*(dist_psi);
Z_psi_q1_down=Z_psi_fit_down(psi_rank_q1,:)*(1-dist_psi)+Z_psi_fit_down(psi_rank_q1+sign_dist,:)*(dist_psi);
% xi_value_q1=xi_psi_fit(psi_rank_q1)*(1-dist_psi)+xi_psi_fit(psi_rank_q1+sign_dist)*dist_psi;
X1_q1=X1_Nradial(psi_rank_q1)*(1-dist_psi)+X1_Nradial(psi_rank_q1+sign_dist)*dist_psi;
X1_q1_rank=round(X1_q1);
X2_q1=X2_Nradial(psi_rank_q1)*(1-dist_psi)+X2_Nradial(psi_rank_q1+sign_dist)*dist_psi;
X2_q1_rank=round(X2_q1);
Z_psi_q1(1:X1_q1_rank)=Z_axis;
Z_psi_q1(X2_q1_rank:NR)=Z_axis;
Z_psi_q1_down(1:X1_q1_rank)=Z_axis;
Z_psi_q1_down(X2_q1_rank:NR)=Z_axis;

% pos_xi_value_q1=round(xi_value_q1/DX)+1;
% X1_center=pos_xi_value_q1+2*NX+1;
Z1_center=2*NX+1;
r_value_q1_mean=radial_r_value_flux(psi_rank_q1)

figure(2);

%%
hold on;
grid on

plot(X_scale(X1_q1_rank:X2_q1_rank),Z_psi_q1(X1_q1_rank:X2_q1_rank),'k','linewidth',2);
plot(X_scale(X1_q1_rank:X2_q1_rank),Z_psi_q1_down(X1_q1_rank:X2_q1_rank),'k','linewidth',2);


