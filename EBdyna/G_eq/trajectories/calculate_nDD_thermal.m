function [output]=calculate_nDD_thermal(output,input,time_stamp)
global dim const 

CALC_MARKERS=input.N_job;
    % Bosch-Hale coefficients for sigma of DD from beam target
    A1_NDD =   53701;
    A2_NDD =   330.2700;
    A3_NDD =  -0.1271;
    A4_NDD =   2.9327e-05;
    A5_NDD =  -2.5151e-09;

    
% fields required from initial maps:            dim.x_ni dim.Ti_x_ni dim.ni_markers_weight

Ti_n_thermal=zeros(CALC_MARKERS,1);
ni_dV_n_thermal=zeros(CALC_MARKERS,1);
ndd_val_th=zeros(CALC_MARKERS,1);

x_n_thermal=zeros(CALC_MARKERS,3);
v_n_th=zeros(CALC_MARKERS,3);

vXYZ_MB_Ti_th=zeros(CALC_MARKERS,3);
vXYZ_MB_Ti_th2=zeros(CALC_MARKERS,3);
vXYZ_COM_th=zeros(CALC_MARKERS,3);
v_rel_th=zeros(CALC_MARKERS,3);

SELECT_STATS=randi(dim.LENGTH_NI_STATS,CALC_MARKERS,1);
Ti_n_thermal=dim.Ti_x_ni(SELECT_STATS);
ni_n_thermal=dim.ni_val(SELECT_STATS);
ni_dV_n_thermal=dim.ni_markers_weight(SELECT_STATS);
x_n_thermal=dim.x_ni(SELECT_STATS,:);
sigma_MB_Ti=sqrt(Ti_n_thermal.*(const.eV/const.mBulk));

vXYZ_MB_Ti_th(:,1)	= random('Normal',0,sigma_MB_Ti);
vXYZ_MB_Ti_th(:,2)	= random('Normal',0,sigma_MB_Ti);
vXYZ_MB_Ti_th(:,3)	= random('Normal',0,sigma_MB_Ti);
vXYZ_MB_Ti_th2(:,1)	= random('Normal',0,sigma_MB_Ti);
vXYZ_MB_Ti_th2(:,2)	= random('Normal',0,sigma_MB_Ti);
vXYZ_MB_Ti_th2(:,3)	= random('Normal',0,sigma_MB_Ti);

% COM velocities
vXYZ_COM_th = (const.mBulk.*vXYZ_MB_Ti_th(:,:)+const.mBulk.*vXYZ_MB_Ti_th2(:,:))./(2*const.mBulk);
v_rel_th = sqrt(sum((vXYZ_MB_Ti_th2(:,:) - vXYZ_MB_Ti_th(:,:)).^2,2));
% neutron source in COM referential
theta_iso = 2.*rand(CALC_MARKERS,1)-1;
theta_iso = asin(theta_iso);
alpha_iso =  (2*pi).*rand(CALC_MARKERS,1);
vXYZ_neutrons(:,1)    = cos(theta_iso).*cos(alpha_iso);
vXYZ_neutrons(:,2)    = cos(theta_iso).*sin(alpha_iso);
vXYZ_neutrons(:,3)    = sin(theta_iso);
vXYZ_neutrons(:,:)    = const.v0_nDD.*vXYZ_neutrons(:,:);
% total velocity in laboratory referential [cartesians!]
v_n_th = vXYZ_neutrons(:,:) + vXYZ_COM_th(:,:);

% recalculate the neutron rate for the 2 maxwellians
E_sigv=(0.5).*(0.5*const.mBulk/const.eV)*sum(v_rel_th(:,:).^2,2).*1e-3;
% cross section in center of mass for D
s_E=(A1_NDD+E_sigv.*(A2_NDD+E_sigv.*(A3_NDD+E_sigv.*(A4_NDD+E_sigv.*A5_NDD))));
sv_val=s_E./(E_sigv.*exp(31.397./sqrt(E_sigv)));  % in millibarns (10^-31 m^2)
sv_val=sv_val.*v_rel_th*1e-31;
% half factor from BH formula with Kronecker symbol (Bosch-Hale
% paper) : each ion marker sees the half density
% ndd_val_th=(ni_dV_n_thermal.*ni_n_thermal).*sv_val;       % neutron rate you would have over 1s (absolute, not volumic)
ndd_val_th=0.5.*(ni_dV_n_thermal.*ni_n_thermal).*sv_val;       % neutron rate you would have over 1s (absolute, not volumic)


output.x_n_thermal(:,:,time_stamp)=x_n_thermal(:,:);
output.v_n_thermal(:,:,time_stamp)=v_n_th(:,:);
output.ndd_val_thermal(:,time_stamp)=ndd_val_th(:);
