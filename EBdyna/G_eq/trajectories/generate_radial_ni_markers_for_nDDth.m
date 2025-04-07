function [x_ni Ti_m ni_m ni_markers_weight]=generate_radial_ni_markers_for_nDDth(dim,maps,input)
% dim = 
%                 mid_X: 189
%                 mid_Z: 392
%             mid_Zzero: 392
%                    DX: 0.0037
%               scale_X: [1x376 double]
%               scale_Z: [1x782 double]
%                NB_PSI: 129
%                    R0: 0.9002
%                Z_axis: -5.4007e-05
%             mid_Xzero: 189
%                X_axis: 0.0381
%             psi_scale: [1x129 double]
%          scale_psi_mp: [176x1 double]
%                size_X: 376
%                size_Z: 782
%     psi_scale_correct: [1x129 double]
%              psi_norm: [81x1 double]
%             vtor_prof: [81x1 double]
%               ni_prof: [81x1 double]
%               Ti_prof: [81x1 double]
%               Te_prof: [81x1 double]
%               ne_prof: [81x1 double]
%           pn_nDD_bins: [1x28 double]
%            nfast_prof: [81x1 double]
%             nimp_prof: [81x1 double]
%                vessel: [1x1 struct]
%                    DZ: 0.0037
%                DX_inv: 268.5053
%                DZ_inv: 268.5053
%          func_Z_cross: @(R)polyval(par.offset_midplane,R)
% dim.vessel
% ans = 
%     R_inner_vessel: [1x62 double]
%     R_lower_vessel: [1x53 double]
%     R_outer_vessel: [1x37 double]
%     R_upper_vessel: [1x52 double]
%     Z_inner_vessel: [1x62 double]
%     Z_lower_vessel: [1x53 double]
%     Z_outer_vessel: [1x37 double]
%     Z_upper_vessel: [1x52 double]
%         wall_RZmap: [376x782 double]
%         
% maps = 
%               psi_global: 0.0746
%                   psi_XZ: [376x782 double]
%              psi_norm_XZ: [376x782 double]
%             psi_norm1_XZ: [376x782 double]
%          theta_normal_XZ: [376x782 double]
%     theta_phase_shift_XZ: [376x782 double]
%                     B_2D: [376x782x3 double]
%                     
N_stats=dim.LENGTH_NI_STATS;

x_ni=nan(N_stats,3);
% ni_markers_weight=0;

% loop on markers
DIMX_MAX=max(dim.vessel.R_outer_vessel)-dim.R0;
DIMX_MIN=min(dim.vessel.R_inner_vessel)-dim.R0;
DIMZ_MAX=max(dim.vessel.Z_upper_vessel);
DIMZ_MIN=min(dim.vessel.Z_lower_vessel);

AVOID_BOUNDARY=0.95;

X_m=rand(N_stats,1).*(DIMX_MAX-DIMX_MIN)*AVOID_BOUNDARY+DIMX_MIN;
Z_m=rand(N_stats,1).*(DIMZ_MAX-DIMZ_MIN)*AVOID_BOUNDARY+DIMZ_MIN;
phi_m=rand(N_stats,1).*2*pi;

psi_m=interp2(dim.scale_X,dim.scale_Z,maps.psi_norm1_XZ',X_m,Z_m,'*linear');

REJECTED=(psi_m>1);


while sum(REJECTED)>0
    X_m(REJECTED)=rand(sum(REJECTED),1).*(DIMX_MAX-DIMX_MIN)*AVOID_BOUNDARY+DIMX_MIN;
    Z_m(REJECTED)=rand(sum(REJECTED),1).*(DIMZ_MAX-DIMZ_MIN)*AVOID_BOUNDARY+DIMZ_MIN;
    psi_m=interp2(dim.scale_X,dim.scale_Z,maps.psi_norm1_XZ',X_m,Z_m,'*linear');
    REJECTED=((psi_m>1)|isnan(psi_m));
end
psi_m=interp2(dim.scale_X,dim.scale_Z,maps.psi_norm1_XZ',X_m,Z_m,'*cubic');
ni_m=interp1(dim.psi_norm,dim.ni_prof',psi_m,'pchip');
DV_m=interp1(dim.psi_norm,dim.DV_prof',psi_m,'pchip');
Ti_m=interp1(dim.psi_norm,dim.Ti_prof',psi_m,'pchip');

N_METIS_RADIAL_POINTS=21

% ni_markers_weight_raw=ni_m.*DV_m; % crude approximation of local number of markers
% beware 22 is number of points in METIS +1 
ni_weight_raw=histc(psi_m,linspace(0,1,N_METIS_RADIAL_POINTS+1));  % crude approximation of local number of markers
ni_markers_number=interp1(dim.psi_norm(1:N_METIS_RADIAL_POINTS),ni_weight_raw(1:N_METIS_RADIAL_POINTS),psi_m,'pchip');
ni_dvol=ni_m.*DV_m;
ni_markers_weight=ni_dvol./ni_markers_number;
% total number of ion in plasma is known
% ni_markers_weight=trapz(dim.ni_prof(1:N_METIS_RADIAL_POINTS).*dim.DV_prof(1:N_METIS_RADIAL_POINTS)).*ni_markers_weight./sum(ni_markers_weight);



x_ni_cyl = [X_m Z_m phi_m];

x_ni(:,1)=x_ni_cyl(:,1).*cos(x_ni_cyl(:,3));
x_ni(:,2)=-x_ni_cyl(:,1).*sin(x_ni_cyl(:,3));
x_ni(:,3)=x_ni_cyl(:,2);
