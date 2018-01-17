
% sqrtg_S_avg=mean((sqrtg_B_PR_map(1:Nomega-1,:))./(Btor_PR_map(1:Nomega-1,:)),1);
% Dr_avg=mean((dist_surf_PR_map(1:Nomega-1,:)),1);
% 
% clear nu_surf_radius;
% r_count=1;
% nu_surf_radius(1)=0;
% 
% for (omega=1:Nomega)
%     for (r=1:pos_psi_rx)
%         r_value=scale_r(r);
%         omega_value=scale_omega(omega);
%         if (r<=pos_psi_rx) && (map_region1(r,omega)~=1);
%             Dr1=Dr_avg(r);
%             if (omega==1) || (omega==Nomega)
%                 nu_surf_radius(r_count)=sqrtg_S_avg(r)*Domega*Dr1;
%                 %nu_surf_radius(r_count)=Domega*Dr1;
%             else
%                 Dr2=Dr1;
%                 nu_surf_radius(r_count)=sqrtg_S_avg(r)*Domega*Dr1+sqrtg_S_avg(r)*Domega*Dr2;
%                 %nu_surf_radius(r_count)=Domega*Dr1+Domega*Dr2;
%             end
%             r_count=r_count+1;
%         end
%     end
% end
% 
%area_region3=sum(nu_surf_radius);
area_region1=interp1(1:NR,surf_flux_precise(1,1:NR)',x_nu_cont13_initial,'cubic');
area_region2=interp1(1:NR,surf_flux_precise(1,1:NR)',pos_psi_rx,'cubic');
%area_region3=interp1(1:size_r,surf_flux(1,1:size_r)',x_psi_final_rx,'cubic');

area_region3=area_region2-area_region1
x_psi_final_rx=interp1(surf_flux_precise(1:NR),(1:NR),area_region3,'cubic');