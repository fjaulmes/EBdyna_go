function [ v Eexch_e Eexch_i Etot_e Etot_i delta_omega] = calculate_collisions(const,par,dim,X_SLOW_ELEC_THRESH,input,...
    Ti,vth_i,ni,psi_gc_value,psi_xi_pos,psi_slope,v_norm,Ekin,v,B,Bfield,Eexch_e,Eexch_i,Etot_e,Etot_i,TESTING_SD)
%CALCULATE_COLLISIONS gives you new velocities according to collision rates
%security on maximum energy loss in one collision!
EKIN_FRAC_LOSS_MAX=0.8;
V_FRAC_LOSS_MAX=sqrt(EKIN_FRAC_LOSS_MAX);

%   Detailed explanation goes here
%             ni=interp1qr(dim.psi_norm',dim.ni_prof',psi_gc_value,psi_xi_pos,psi_slope);
ne=interp1qr(dim.pn,dim.ne_prof,psi_gc_value,psi_xi_pos,psi_slope);
Te=interp1qr(dim.pn,dim.Te_prof,psi_gc_value,psi_xi_pos,psi_slope);
%             vth_e=interp1qr(dim.psi_norm',dim.vthe_prof',psi_gc_value,psi_xi_pos,psi_slope);
%             vth_i=interp1qr(dim.psi_norm',dim.vthi_prof',psi_gc_value,psi_xi_pos,psi_slope);
vth_e=sqrt(3*Te*(const.eV/const.me));


% collision rates without log(lambda) 
nu_e=(1/(4*pi)).*input.Z.^2.*(const.eV^2/const.epsilon0)^2./(input.m*const.me).*ne./(vth_e.^3); % simplified by mf/me already
nu_i=(1/(4*pi)).*input.Z.^2.*(const.eV^2/const.epsilon0)^2./(input.m.^2).*ni;

% log(lambda) : original EBdyna paper, modified with new effective velocity for ions


% RABBIT formulas (without Zbulk)
% mass number in Rabbit seems wrt. to approximately mH
% modified slightly the eV^2 terms
veff_e_sq1=(Te.*const.eV./const.me+v_norm.^2);
veff_i_sq1=(Ti.*const.eV./const.mD+v_norm.^2);
omega_ci_sq=1.74.*(1./const.ABulk).*ni+(9.163e15/const.ABulk^2).*Bfield.^2;
omega_ce_sq=1.74.*(1./const.Ae).*ni+(9.163e15/const.Ae^2).*Bfield.^2;
rmax=((omega_ce_sq./veff_e_sq1)+(omega_ci_sq./veff_i_sq1)).^-0.5;

veff_e_sq2=(3.*Te.*const.eV./const.me+v_norm.^2);
veff_i_sq2=(3.*Ti.*const.eV./const.mD+v_norm.^2);

rmin_cl_Bulk=0.13793.*(input.Z).*(const.ABulk+input.Afast)./(input.Afast.*const.ABulk.*veff_i_sq2);
rmin_cl_e=0.13793.*(input.Z).*(const.Ae+input.Afast)./(input.Afast.*const.Ae.*veff_e_sq2);
rmin_qu_Bulk=1.9121e-8.*(const.ABulk+input.Afast)./(input.Afast.*const.ABulk.*sqrt(veff_i_sq2));
rmin_qu_e=0.9121e-8.*(const.Ae+input.Afast)./(input.Afast.*const.Ae.*sqrt(veff_e_sq2));
rmin_Bulk=max(rmin_cl_Bulk,rmin_qu_Bulk);
rmin_e=max(rmin_cl_e,rmin_qu_e);

log_lambda_Fi=log(rmax./rmin_Bulk);
log_lambda_Fe=log(rmax./rmin_e);

% based on S. Ward paper:
%log_lambda_De=0.5*log(Te./(ne.*const.eV./const.epsilon0+(const.eV/const.me).*Bfield.^2));
%mr_i=(const.mBulk.*input.m)./(const.mBulk+input.m);
%log_lambda_Fi=log_lambda_De-log(const.eV^2)+log(4*pi*const.epsilon0*mr_i.*veff_i_sq);
%log_lambda_Fe=log_lambda_De-log(const.eV^2)+log(4*pi*const.epsilon0*const.me.*veff_e_sq);


% log_lambda_Fe=15.2-0.5*log(ne*1e-20)+log(0.5.*1e-3.*(const.me/const.eV).*veff_e_sq);
%log_lambda_Fi=log_lambda_De-log((const.eV^2)./(4*pi*const.epsilon0*mr_i.*(v_norm).^2));
%log_lambda_Fe=15.2-0.5*log(ne*1e-20)+log(0.5.*1e-3.*(const.me/const.eV).*(max(vth_e,v_norm)).^2);



%lambda_De = sqrt((const.epsilon0.*Te)./(ne.*const.eV));
%b90_D=(const.eV^2)./(4*pi*const.epsilon0*0.5*const.mD.*veff_i_sq);
%b90_e=(const.eV^2)./(4*pi*const.epsilon0*const.me.*veff_e_sq);
%log_lambda_ie_ASCOT=log(lambda_De./b90_e);
%log_lambda_ii_ASCOT=log(lambda_De./b90_D);


if par.REMOVE_SOL_PLASMA
    log_lambda_Fe(ne<1)=0;
    log_lambda_Fi(ne<1)=0;
end

% vth_e_sq=(2*Te*(const.eV/const.me));
% v_crit=sqrt(2*18.66.*Te*(const.eV/input.m));
v_crit=(3*sqrt(pi)/4.*const.me./input.m).^(1/3).*vth_e;  % = sqrt(2*18.66.*Te*(const.eV/input.m));
% 
		   
% v_crit = 5.33*1e4  .* sqrt(Te) .* const.A_i.^-1/3
% E_crit = Te .*( 3 * sqrt(pi)/4 * sqrt(input.m/const.me) * input.m/const.mBulk).^2/3;
% v_crit= sqrt(2.*E_crit./input.m);

%%
% omega_e2 = 1.74.*ne./const.A_e+9.18.*1e15.*(1./const.A_e)^2.*B_field2;
% omega_i2 = 1.74.*ni./const.A_i+9.18.*1e15.*(1./const.A_i)^2.*B_field2;
% vrel_e2  = 9.58*1e7*(Te./const.A_e+2*Ekin./input.A_NB);
% vrel_i2  = 9.58*1e7*(Ti./const.A_i+2*Ekin./input.A_NB);
% rmax=sqrt(1./(omega_e2./vrel_e2+omega_i2./vrel_i2));
% vrel2_e2  = 9.58*1e7*(3*Te./const.A_e+2*Ekin./input.A_NB);
% vrel2_i2  = 9.58*1e7*(3*Ti./const.A_i+2*Ekin./input.A_NB);
% 
% rmin1_e=0.13793.*(input.Z).*(const.A_e+input.A_NB)./(const.A_e.*input.A_NB.*vrel2_e2);
% rmin1_i=0.13793.*(input.Z).*(const.A_i+input.A_NB)./(const.A_i.*input.A_NB.*vrel2_i2);
% rmin2_e=1.9121.*1e-8.*(input.Z).*(const.A_e+input.A_NB)./(const.A_e.*input.A_NB.*sqrt(vrel2_e2));
% rmin2_i=1.9121.*1e-8.*(input.Z).*(const.A_i+input.A_NB)./(const.A_i.*input.A_NB.*sqrt(vrel2_i2));
% rmin_e=max(rmin1_e,rmin2_e);
% rmin_i=max(rmin1_i,rmin2_i);
% log_lambda_De=log(rmax./rmin_e);
% log_lambda_Di=log(rmax./rmin_i);

%%
%securiy
%log_lambda_Di=max(log_lambda_Di,1);

nu_e=nu_e.*log_lambda_Fe;
nu_i=nu_i.*log_lambda_Fi./(v_norm.^3);
% nu_i=nu_i.*log_lambda_Di;

% to get some slowing down in a much shorter time....
if TESTING_SD
    nu_e=nu_e*50;
    nu_i=nu_i*50;
end

% Spreading identical to NUBEAM calculation : using order 0 values of collision freq.
% delta_omega=(v_crit./v_norm).^3;
% delta_v_spread=sqrt((par.NR_FUND_IN_LOOP.*par.dt.*nu_e).*(2*const.eV/input.m).*(Te+Ti.*delta_omega));
% considering only simple ion collision term
%delta_omega=(pi/2)*sqrt(par.dt.*par.NR_FUND_IN_LOOP.*nu_e.*delta_omega);

% we discriminate between region of fast and slow electrons
sqrt_x_v=(v_norm./vth_e);
x_v=sqrt_x_v.^2;
SLOW_ELEC_GROUP=x_v>X_SLOW_ELEC_THRESH;
% if sum(SLOW_ELEC_GROUP)>=1 & TESTING_SD
%     disp('testing SLOW_ELEC_GROUP....')
% end
%             SLOW_ELEC_GROUP=COLL_GROUP&SLOW_ELEC_GROUP;

psi_x=x_v*0;

% beware psi expressions modified wrt. reference (JD Callen Chap. 2)
% because of the xpression usd for nu_e that is also different

% 4/3*sqrt(pi) = 0.752252778
psi_x(~SLOW_ELEC_GROUP)=0.752252778.*(1-3.*x_v(~SLOW_ELEC_GROUP)/5+3.*x_v(~SLOW_ELEC_GROUP).^2/14);  % here all multiplicative terms were in the nu_e fit
% at large velocity (compared with thermal electrons) we need to remove the term in 1/vthe^3 from
% nue (and put back the 1/v^3)
psi_x(SLOW_ELEC_GROUP) =(sqrt_x_v(SLOW_ELEC_GROUP).^-3).*(1-(2*(sqrt_x_v(SLOW_ELEC_GROUP)).*exp(-x_v(SLOW_ELEC_GROUP))./sqrt(pi)).*(1+1./(2.*x_v(SLOW_ELEC_GROUP))-1./(4.*x_v(SLOW_ELEC_GROUP).^2)));

%           psi prime term is negligible for ion -> electrons
%             psi_prime_x=(2*(sqrt_x_v).*exp(-x_v)./sqrt(pi));
%             pmpp_x=(mD/me)*psi_x-psi_prime_x;

nu_e=nu_e.*psi_x;
% nu_e=nu_e.*(1./(1+0.75*(v_norm./vth_e).^3));
% gaussian spread of velocity chosen consistent with our Freidberg electron collision rate
if par.SLOWING_DOWN
    delta_omega=Ti.*(v_crit./v_norm).^3; % intermediate variable (name not consistent)
else
    delta_omega=Ti.*(v_crit./max(v_norm,v_crit)).^3; % otherwise this term is causing chaos in velocity spread
end
if par.CALCULATE_PDEP & par.SLOWING_DOWN
    dv_e=sqrt((par.NR_FUND_IN_LOOP.*par.dt.*nu_e).*(2*const.eV/input.m).*Te);
    dv_i=sqrt((par.NR_FUND_IN_LOOP.*par.dt.*nu_e).*(2*const.eV/input.m).*delta_omega);
    %possible renormalization of the total chaos
%     delta_v_spread=sqrt((par.NR_FUND_IN_LOOP.*par.dt.*nu_e).*(2*const.eV/input.m).*(Te+delta_omega));
%     dv_e=dv_e.*delta_v_spread./(dv_e+dv_i);
%     dv_i=dv_i.*delta_v_spread./(dv_e+dv_i);
    %                 delta_v_spread=dv_e+dv_i;
else
    dv_e=v_norm.*0;
    dv_i=v_norm.*0;
    delta_v_spread=sqrt((par.NR_FUND_IN_LOOP.*par.dt.*nu_e).*(2*const.eV/input.m).*(Te+delta_omega));
end

% corrected NBI ion=> ion frequency (Freidberg)
%             nu_i=interp1qr(dim.psi_norm',dim.nuNBii_v3',psi_gc_value,psi_xi_pos,psi_slope);
%             vth_i=interp1qr(dim.psi_norm',dim.vthi_prof',psi_gc_value,psi_xi_pos,psi_slope);
% Freidberg expression
%             nu_i=nu_i./(v_norm.^3+1.33*vth_i.^3);

% we have a profile for MOMENTUM exchange

%             % for shorter simulation!!! (DEBUG)
%             nu_i=10*nu_i;
%             nu_e=10*nu_e;

sqrt_x_v=(v_norm./vth_i);
x_v=sqrt_x_v.^2;

FAST_ION_GROUP=x_v>X_SLOW_ELEC_THRESH;


psi_x(FAST_ION_GROUP)=1-(2*(sqrt_x_v(FAST_ION_GROUP)).*exp(-x_v(FAST_ION_GROUP))./sqrt(pi)).*(1+1./(2.*x_v(FAST_ION_GROUP))-1./(4.*x_v(FAST_ION_GROUP).^2));
psi_x(~FAST_ION_GROUP)=0.752252778.*(x_v(~FAST_ION_GROUP).^1.5).*(1-3.*x_v(~FAST_ION_GROUP)/5+3.*x_v(~FAST_ION_GROUP).^2/14);

psi_prime_x=(2*(sqrt_x_v).*exp(-x_v)./sqrt(pi));

pmpp_x=(input.m./const.mBulk).*psi_x-psi_prime_x;
%             nu_perp=nu_i.*(psi_x+psi_prime_x-psi_x./(2.*x_v));
%nu_perp=nu_i.*((psi_x+psi_prime_x-psi_x./(2.*x_v))./pmpp_x);
nu_perp=nu_i.*((psi_x+psi_prime_x-psi_x./(2.*x_v)));
% nu_pll=nu_i.*(psi_x./x_v);

% for consistency with formulation at low velocities :  
% beware, energy transfer is twice larger than reference collision rate, 
% BUT when going back o the norm of velocity, the factor 2 is cancelled
nu_i=nu_i.*pmpp_x;
%             plot(x,(psi+psi_prime)./pmpp_x-psi./(2.*x.*pmpp_x))

%             nu_E=nu_i.*(psi_x-psi_prime_x);

% gaussian spread for velocities : can be sped up by
% precomputation of normal distribution values
% 2*sqrt(2*log(2)) = 2.3548 and sqrt(2*log(2))=1.1774

% note: the collision frequency of interest are that of Energy exchange
% see equation (2.97) of JD Callen Chap. 2

if par.CALCULATE_PDEP
    if par.SLOWING_DOWN
        dv_e=par.NR_FUND_IN_LOOP.*par.dt.*v_norm.*nu_e+randn(length(v_norm),1).*dv_e;
        dv_i=par.NR_FUND_IN_LOOP.*par.dt.*v_norm.*nu_i+randn(length(v_norm),1).*dv_i;
        dv_e=max(dv_e,-V_FRAC_LOSS_MAX.*v_norm);
        dv_i=max(dv_i,-V_FRAC_LOSS_MAX.*v_norm);
    end
    delta_v=dv_e+dv_i;
    % energy lost is given to plasma species
    if par.SLOWING_DOWN
        dv_e=0.5*(input.m/const.eV).*(2.*v_norm.*dv_e-dv_e.^2);
        dv_i=0.5*(input.m/const.eV).*(2.*v_norm.*dv_i-dv_i.^2);
%         dv_e=max(dv_e,-EKIN_FRAC_LOSS_MAX.*Ekin);
%         dv_i=max(dv_i,-EKIN_FRAC_LOSS_MAX.*Ekin);
    end
    %if par.CALCULATE_CX
    %    % moderate deposited power according to CX losses
    %    dv_e=dv_e.*CX_val;
    %    dv_i=dv_i.*CX_val;
    %end
	
	% security (randn created a few imaginary values at low velocities)
	dv_e=real(dv_e);
	dv_i=real(dv_i);
	
    Eexch_e=Eexch_e+dv_e;
    Eexch_i=Eexch_i+dv_i;
    Etot_e=Etot_e+dv_e;
    Etot_i=Etot_i+dv_i;
else
    nu_tot=nu_e+nu_i;
    if par.SLOWING_DOWN
        delta_v=par.NR_FUND_IN_LOOP.*par.dt.*v_norm.*nu_tot;
    else
        delta_v=v_norm.*0;
    end
    delta_v=delta_v+randn(length(v_norm),1).*delta_v_spread;
end
%delta_vpll=randn(length(v_norm),1).*sqrt(par.dt.*par.NR_FUND_IN_LOOP.*nu_pll);
  
  
% security (randn created a few imaginary values at low velocities)
delta_v=real(delta_v);
%vpll= dot(B,v,2)./Bfield;
%delta_vpll=randn(length(v_norm),1).*sqrt(par.dt.*par.NR_FUND_IN_LOOP.*nu_pll).*vpll;



%small change along parallel direction
%v=bsxfun(@minus,v,delta_vpll);

v=bsxfun(@times,(v_norm-delta_v)./v_norm,v);


% spread in angle based on normal distribution
delta_omega=(pi/2).*randn(length(v_norm),1).*sqrt(par.dt.*par.NR_FUND_IN_LOOP.*nu_perp);
%delta_omega=(pi/2).*randn(length(v_norm),1).*sqrt(par.dt.*par.NR_FUND_IN_LOOP.*nu_perp);
%delta_omega=(pi/2).*randn(length(v_norm),1).*sqrt(par.dt.*par.NR_FUND_IN_LOOP.*nu_perp)+randn(length(v_norm),1).*sqrt(par.dt.*par.NR_FUND_IN_LOOP.*nu_i);

% security (randn created a few imaginary values at low velocities)
delta_omega=real(delta_omega);

end

