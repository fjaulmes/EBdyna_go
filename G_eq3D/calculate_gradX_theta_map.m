

% corrected theta maps for complete interpolation
QNB_THETA=round(0.25*NB_THETA);
HQNB_THETA=round(0.5*QNB_THETA);
theta_low_XZsmall_map=theta_XZsmall_map;
theta_up_XZsmall_map=theta_XZsmall_map;
theta_up_XZsmall_map(find(theta_XZsmall_map<QNB_THETA*DTHETA))=theta_XZsmall_map(find(theta_XZsmall_map<QNB_THETA*DTHETA))+2*pi;
theta_low_XZsmall_map(find(theta_XZsmall_map>(NB_THETA-QNB_THETA-2)*DTHETA))=theta_XZsmall_map(find(theta_XZsmall_map>(NB_THETA-QNB_THETA-2)*DTHETA))-2*pi;


gradX_theta_map=zeros(size_X,size_Z);

for (x=3:size_X-2)
    for (z=3:size_Z-2)
%         if theta_XZsmall_map(x,z)<HQNB_THETA*DTHETA
%             grad_theta_Z(x,z)=(1/12)*(-theta_up_XZsmall_map(x,z+2)+theta_up_XZsmall_map(x,z-2))+(2/3)*(theta_up_XZsmall_map(x,z+1)-theta_up_XZsmall_map(x,z-1));
%         elseif theta_XZsmall_map(x,z)>(NB_THETA-HQNB_THETA-2)*DTHETA
%             grad_theta_Z(x,z)=(1/12)*(-theta_low_XZsmall_map(x,z+2)+theta_low_XZsmall_map(x,z-2))+(2/3)*(theta_low_XZsmall_map(x,z+1)-theta_low_XZsmall_map(x,z-1));
        if (theta_XZsmall_map(x,z)<HQNB_THETA*DTHETA)||(theta_XZsmall_map(x,z)>(NB_THETA-HQNB_THETA-2)*DTHETA)
            gradX_theta_map(x,z)=(1/12)*(-theta_low_XZsmall_map(x+2,z)+theta_low_XZsmall_map(x-2,z))+(2/3)*(theta_low_XZsmall_map(x+1,z)-theta_low_XZsmall_map(x-1,z));
        else
            gradX_theta_map(x,z)=(1/12)*(-theta_XZsmall_map(x+2,z)+theta_XZsmall_map(x-2,z))+(2/3)*(theta_XZsmall_map(x+1,z)-theta_XZsmall_map(x-1,z));
        end
    end
end

gradX_theta_map(isnan(gradX_theta_map))=0;

gradX_theta_map=gradX_theta_map/DX;

% figure(7)
% set(gca,'FontSize',26);
% % title('\theta');
% 
% hold on
% grid on
% contour(scale_X,scale_Z,theta_XZsmall_map',(0:0.5:6));axis xy;
% colormap('jet')
% 
% imagesc(scale_X,scale_Z,psi_norm_XZsmall_map')
% colormap('gray')
% 
% contour(scale_X,scale_Z,theta_XZsmall_map',(0:0.5:6));axis xy;
% 
% % contour(scale_X,scale_Z,theta_XZsmall_map',[0 0],'k','LineWidth',4);axis xy
% 
% xlabel('X (m)');ylabel('Z (m)')
% axis equal;
% % 
% % 
figure(8)
set(gca,'FontSize',26);
% title('gradZ(\theta)');

hold on
grid on
colorbar

contour(scale_X,scale_Z,gradX_theta_map',[0 0],'k','LineWidth',5);axis xy
contour(scale_X,scale_Z,gradX_theta_map',(-3:1:3),'--','LineWidth',2);axis xy;
axis equal;
xlim([-0.5 1])
ylim([-0.9 0.9])
xlabel('X (m)');ylabel('Z (m)')

% 
% SMALL_ORBITS=find(CO_PASSING_POP'.*(r_avg>0.34).*(r_avg<0.36));
% plot(X_output_pop(1:144,SMALL_ORBITS(10)),Z_output_pop(1:144,SMALL_ORBITS(10)),'r','linewidth',3)
% SMALL_ORBITS=find(COUNTER_PASSING_POP'.*(r_avg>0.3).*(r_avg<0.32));
% plot(X_output_pop(1:144,SMALL_ORBITS(10)),Z_output_pop(1:144,SMALL_ORBITS(10)),'b--','linewidth',5)
% plot(X_output_pop(1:144,ALL_TRAPPED(90)),Z_output_pop(1:144,ALL_TRAPPED(90)),'g--','linewidth',6)
% % 
