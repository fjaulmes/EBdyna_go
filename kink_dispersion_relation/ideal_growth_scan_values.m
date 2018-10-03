reset_data_analysis_environment
% P_initial_profile=PTOT_profile_interp_ini;
% load('NBI_Phot_data.mat')
load('deltaW_NBI_data.mat', 'deltaW_fast_hat')

% rescaling to account for the xi^2 normalizatio 
% (a bit higher than the average xi value used in the kink displacement calculation)
deltaW_fast_hat=deltaW_fast_hat

Bphi0=Bavg;
r_value_q1_mean=interp1(q_initial_profile,radial_r_value_flux,1)


TE_profile_radial_ini=interp1(psi_bar_rho,TE_profile_interp_ini,(0:Nradial-1)/(Nradial-1));
TI_profile_radial_ini=interp1(psi_bar_rho,TI_profile_interp_ini,(0:Nradial-1)/(Nradial-1));
NE_profile_radial_ini=interp1(psi_bar_rho,NE_profile_interp_ini,(0:Nradial-1)/(Nradial-1));
NI_profile_radial_ini=interp1(psi_bar_rho,NI_profile_interp_ini,(0:Nradial-1)/(Nradial-1));


%
Te_profile=TE_profile_radial_ini*eV;
Ne_profile=NE_profile_radial_ini;
Ne0=NE_profile_radial_ini(1)
Te0=TE_profile_radial_ini(1)*eV

%
Pbulk_profile=(NI_profile_radial_ini.*TI_profile_radial_ini*eV+NE_profile_radial_ini.*TE_profile_radial_ini*eV);

Btot_PR_map=sqrt(BX_PR_map.^2+BZ_PR_map.^2+Btor_PR_map.^2);
derq=gradient(q_initial_profile,radial_r_value_flux);
derP=gradient(P_initial_profile,radial_r_value_flux);
derPi=gradient((NI_profile_radial_ini.*TI_profile_radial_ini*eV),radial_r_value_flux);
derPe=gradient((NE_profile_radial_ini.*TE_profile_radial_ini*eV),radial_r_value_flux);

derTe=gradient(Te_profile,radial_r_value_flux);
derne=gradient(Ne_profile,radial_r_value_flux);
derni=gradient(NI_profile_radial_ini,radial_r_value_flux);

B1=mean(Btot_PR_map(1:end-1,psi_rank_q1))
Ti1=TI_profile_radial_ini(psi_rank_q1)*eV;
Te1=Te_profile(psi_rank_q1);
ne1=Ne_profile(psi_rank_q1)
ni1=NI_profile_radial_ini(psi_rank_q1)
vA1=B1/sqrt(mu0*mD*ni1);
derq1=derq(psi_rank_q1);
omegaA1_tilde=vA1/(sqrt(3)*r_value_q1_mean*R0*derq1);
omegaA1=vA1/R0

omegaci=eV*B1/mD;
%be careful, Pi is half of the total pressure - NBI contribution
Pe0=Ne0*Te0
derPi1=derPi(psi_rank_q1);
derPe1=derPe(psi_rank_q1);
derTe1=derTe(psi_rank_q1);
wdia=(1/(mD*ni1*r_value_q1_mean*omegaci))*derPi1
% wdia=(T1/(mD*n1*r_value_q1_mean*omegaci))*dern1


C0=sqrt(1/(epsilon0*mu0));
omegape=sqrt((ne1*eV^2)/(epsilon0*me))
delta_e=C0/omegape

we=-((1/(mD*ne1*r_value_q1_mean*omegaci))*derPe1+0.71*(1/(mD*r_value_q1_mean*omegaci))*derTe1)

% we=-((T1/(mD*r_value_q1_mean*omegaci))*dern1/n1+0.71*(1/(mD*r_value_q1_mean*omegaci))*derT1)
% we=-((0.71*(1/(mD*r_value_q1_mean*omegaci))*derT1)
% we=-((T1/(mD*r_value_q1_mean*omegaci*n1))*dern1)

Bpol1sq=mean(Bpol_PR_map(1:end-1,psi_rank_q1).^2);
dr_avg=mean(dr_PR_map(1:end-1,:));
dr_avg(1:end-1)=radial_r_value_flux(2:end)-radial_r_value_flux(1:end-1);
integ1=0;
for r=2:psi_rank_q1
    integ1=integ1-radial_r_value_flux(r)^2*dr_avg(r)*(derP(r));
end
beta_pol1=(2*mu0/Bpol1sq)*integ1/r_value_q1_mean^2

%electron ion collision frequency
lambda_D=sqrt(epsilon0*Te_profile./(Ne_profile*eV^2));
Lambda_ei=9*(4*pi/3*Ne_profile.*lambda_D.^3);
log_lambda=log(Lambda_ei);
tau_ei=(3*(2*pi)^1.5)*(sqrt(me)*(Te_profile).^1.5)./(Ne_profile.*log_lambda*eV^2);
tau_ei=(epsilon0/eV)^2*tau_ei;
nu_ei_1=interp1(1:Nradial,1./tau_ei,psi_rank_q1);
Zeff=1.5;
eta1=0.51*Zeff*me*nu_ei_1/(ne1*eV^2)
tau_eta=mu0*r_value_q1_mean^2/eta1;
shear1=r_value_q1_mean*derq1
omegaA1_tilde=vA1/(sqrt(3)*R0*r_value_q1_mean*derq1);
omegaA1_hat=(vA1*r_value_q1_mean*derq1)/(sqrt(3)*R0);
omegaA1_wesson=(vA1*r_value_q1_mean*derq1)/(R0);
tau_H=1/omegaA1_hat;
epsilon_eta=tau_H/tau_eta;

%  gamma Bussac 
deltaq=1-q_initial_profile(1)
r2=interp1(q_initial_profile,radial_r_value_flux,2)
% r2=r_value_q1_mean*(1/deltaq+1)^(1/3)
% beta_coef=(1-(r_value_q1_mean/r2)^2)
% 
% gammaI_nominal=omegaA1*sqrt(3)*pi*(r_value_q1_mean^2/R0^3)*(1/derq1)*deltaq*(3*beta_coef*beta_pol1^2-13/112)
beta_coef=(1-(r_value_q1_mean/r2))
deltaW_ideal_hat=-(1-q_initial_profile(1))*(3*beta_coef*beta_pol1^2-13/112)

% gammaI_nominal=-omegaA1*sqrt(3)*pi*(r_value_q1_mean/R0^2)*(1/derq1)*deltaW_ideal_hat


xvalues=radial_r_value_flux/radial_r_value_flux(psi_rank_q1);
dx_avg=xvalues*0;
dx_avg(2:psi_rank_q1)=xvalues(2:psi_rank_q1)-xvalues(1:psi_rank_q1-1);

% PI_profile=(NI_profile_radial_ini.*TI_profile_radial_ini*eV);
PI_profile=(P_initial_profile-NE_profile_radial_ini.*TE_profile_radial_ini*eV);

%only half of NBI pressure considered contributing to trapped stabilization
PI_profile=0.5*(PI_profile+NI_profile_radial_ini.*TI_profile_radial_ini*eV);

%we now neglect the fast fraction in the trapped contribution
PI_profile=NI_profile_radial_ini.*TI_profile_radial_ini*eV;

cp_integ=0;
Pnorm=PI_profile/PI_profile(1);
for r=2:psi_rank_q1
    cp_integ=cp_integ+((xvalues(r)+xvalues(r-1))/2)^(3/2)*dx_avg(r)*(Pnorm(r));
end

betai0=2*mu0*PI_profile(1)/mean(Btor_PR_map(:,1))^2

deltaW_ko_hat=betai0*(1.5/(6*pi))*cp_integ*(r_value_q1_mean/R0)^(-3/2)

deltaW_el_hat=-3*(li_profile(psi_rank_q1)-0.5)^3*(0.5*(kappa1-1))^2

deltaW_hat_tot=deltaW_ideal_hat+deltaW_ko_hat+deltaW_fast_hat+deltaW_el_hat

gammaI_nominal=-omegaA1*sqrt(3)*pi*(r_value_q1_mean/R0^2)*(1/derq1)*(deltaW_hat_tot)
deltaW_ideal_ksi0=(deltaW_hat_tot)*6*pi^2*r_value_q1_mean^4*Bavg^2/(mu0*R0^3)


gammaI_nominal_wdia=gammaI_nominal/abs(wdia)

%%
% using approximate gamma Bussac can yield a scan in r1 or in pressure
disp('------------------------------------------------------')
PLOT_SIMPLE_DISP=0

if PLOT_SIMPLE_DISP==1
    wdia=1
    % solving the real side for gammaI<wdia/2
    NUMBER_POINTS=20000;
    clear gammaI_low omega_r1 omega_r2
    
    for n=1:NUMBER_POINTS
        gammaI_low(n)=wdia*(n-1)/NUMBER_POINTS-0.5*wdia;
        delta=wdia^2-4*gammaI_low(n)^2;
        omega_r1(n)=(wdia-sqrt(delta))/2;
        omega_r2(n)=(wdia+sqrt(delta))/2;
    end
    
    figure(1)
    set(gca,'fontsize',20)
    hold on
    grid on
    plot(gammaI_low,omega_r1,'b--','linewidth',3)
    
    
    
    % solving the imaginary side for gammaI>wdia/2
    NUMBER_POINTS=40000;
    clear gammaI omega_i1 omega_i2
    
    for n=1:NUMBER_POINTS
        gammaI(n)=0.5*wdia+wdia*(n-1)/NUMBER_POINTS;
        b2=gammaI(n)^2-(wdia^2)/4;
        omega_i1(n)=-sqrt(b2);
        omega_i2(n)=sqrt(b2);
    end
    
    figure(1)
    plot(gammaI,omega_i1,'r','linewidth',3)
    plot(gammaI_low,gammaI_low,'k-.','linewidth',2)
    
    legend('Re(\omega)','Im(\omega)','\omega=\gamma_{I}')
    
    plot(gammaI_low,omega_r2,'b--','linewidth',3)
    plot(gammaI,omega_i2,'r','linewidth',3)
    
    plot(gammaI,omega_i2*0+0.5,'b--','linewidth',3)
    plot(gammaI_low,omega_r2*0,'r','linewidth',3)
    
    plot(gammaI,gammaI,'k-.','linewidth',2)
    
    
    xlabel('\gamma_{I} / \omega_{*i}')
    ylabel('\omega / \omega_{*i}')
    
end

disp('------------------------------------------------------')


%%
deltaW_ideal_ksi0=(deltaW_ideal_hat+deltaW_ko_hat)*6*pi^2*r_value_q1_mean^4*Bavg^2/(mu0*R0^3)

% Drift tearing mode analysis
% gamma_values=(-1.0:0.01:1.0)*abs(wdia);
gamma_values=(-1.5:0.05:1.5)*abs(wdia);

deltaprime=(pi^2*shear1^2*Bavg^2*r_value_q1_mean^2/R0/mu0)/deltaW_ideal_ksi0
gammaDT0_wesson=0.55*deltaprime^(4/5)*tau_H^(-2/5)*tau_eta^(-3/5)

deltaW_ideal_hat_ksi0_values=-gamma_values/(omegaA1*sqrt(3)*pi*(r_value_q1_mean/R0^2)*(1/derq1));
deltaW_ideal_ksi0_values=deltaW_ideal_hat_ksi0_values*6*pi^2*r_value_q1_mean^4*Bavg^2/(mu0*R0^3);
deltaprime_values=(pi^2*shear1^2*Bavg^2*r_value_q1_mean^2/R0/mu0)./deltaW_ideal_ksi0_values;
gammaDT0_wesson_values=0.55*deltaprime_values.^(4/5)*tau_H^(-2/5)*tau_eta^(-3/5);

lambdaH0=-mu0*deltaW_ideal_ksi0/(4*shear1^2*Bavg^2*(r_value_q1_mean/R0)^2)
% lambdaH=gammaI_nominal/omegaA1_hat
% 
% 
lambdaH_values=-mu0*deltaW_ideal_ksi0_values/(4*shear1^2*Bavg^2*(r_value_q1_mean/R0)^2)
% lambdaH_values=gamma_values/omegaA1_hat;
gammaDT0=((epsilon_eta*omegaA1_hat/2)*(we*(we-wdia)/omegaA1_hat^2)^(-1/3))*(2.22/abs(lambdaH0))^(4/3)
% 
gammaDT_values=gammaDT0*(1-(2/3)*((2.22./abs(lambdaH_values)).^(2/3))*(we*(we-wdia)/omegaA1_hat^2)^(1/3));

figure(2)
hold on
plot(gamma_values,gammaDT_values)
plot(gamma_values,gammaDT0_wesson_values,'r--')



%%
%further drift tearing mode analysis
close all
pause(0.5)

gamma_eta_wesson=(tau_H)^(-2/3)*tau_eta^(-1/3);
% Drift tearing mode analysis
% gamma_values=(-1.5:0.05:1.5)*abs(wdia);
delta_WH_wesson=gamma_values/omegaA1_hat;

real_freq_values=(-0.8:0.002:0.8)*abs(wdia);
rposfre0=find(real_freq_values==0);
real_values=real_freq_values/gamma_eta_wesson;
% real_values=real_freq_values*0;

NBREAL=length(real_freq_values)

% delta_WH_wesson=deltaW_ideal_ksi0_values/(shear1^2*Bavg^2*r_value_q1_mean^2/R0/mu0);


% this requires a scan in real part as well

gamma_norm=(0.001:0.002:5.2);
NBIMAG=length(gamma_norm)

gamma_num_values=zeros(length(gamma_norm),NBREAL);
gamma_den_values=zeros(length(gamma_norm),NBREAL);
omega_i1=zeros(NBREAL,1);
omega_i2=zeros(NBREAL,1);
omega_r1=zeros(length(gamma_norm),1);
omega_r2=zeros(length(gamma_norm),1);
omega_r3=zeros(length(gamma_norm),1);
gamma_sol1=zeros(length(gamma_values),1);
gamma_sol2=zeros(length(gamma_values),1);
%%
% gamma_transition_value1=zeros(length(delta_WH_wesson),1);
% gamma_transition_value2=zeros(length(delta_WH_wesson),1);
% freq_transition_value1=zeros(length(delta_WH_wesson),1);
% freq_transition_value2=zeros(length(delta_WH_wesson),1);

for r=1:NBREAL
    r
%     real_values(r)=real_freq_values(r)/gamma_eta_wesson;
    gamma_hat_range=gamma_norm-i*real_values(r);
    % complex values
    gamma_num_values(:,r)=(gamma_complex(0.25*gamma_hat_range.^(1.5)-0.25))';
    gamma_den_values(:,r)=(gamma_complex(0.25*gamma_hat_range.^(1.5)+1.25))';
end

size(gamma_num_values)

% save gamma_values.mat  gamma_num_values gamma_den_values real_values gamma_norm NBREAL NBIMAG


%%
%solve for 0 real frequency first
n_zero_pos=find(gamma_values==0)

for n=1:length(gamma_values)
    no_solution=0;
    gammaI_ratio_wdia=gamma_values(n)/abs(wdia)
    gammaI_ratio_gamma_eta=gamma_values(n)/gamma_eta_wesson
    gamma_hat_range=gamma_norm;
    gamma_hat_range=gamma_hat_range';
    % solve real(RHS)=1
    
    RHS=1/(gamma_values(n)/8/gamma_eta_wesson)*(gamma_hat_range.^(-1.25)).*(gamma_den_values(:,rposfre0))./(gamma_num_values(:,rposfre0));

    [max_value max_pos]=max(real(RHS));
    max_value=max_value-1
    [min_value min_pos]=min(real(RHS));
    min_value=min_value-1
    ini_val=real(RHS(1))-1;
    if sign(max_value*ini_val)<0
        if sign(min_value*ini_val)<0
            int_pos=min(min_pos,max_pos);
        else
            int_pos=max_pos;
        end
    else
        if sign(min_value*ini_val)<0
            int_pos=min_pos;
        else
            %                  disp('no solutions to imaginary part equation');
            %                  min_value
            %                  max_value
            no_solution=1;
        end
    end
%     RHS=(gamma_values(n)/8/gamma_eta_wesson)*(gamma_hat_range.^(1.25)).*(gamma_num_values(:,rposfre0))./(gamma_den_values(:,rposfre0));
    if no_solution==0
        try
            gamma_sol1(n)=fzero(@(x0) solve_real_un( RHS,gamma_norm,x0 ),[min(gamma_norm) gamma_norm(int_pos)]);
        catch
            gamma_sol1(n)=0;
        end
        try
            gamma_sol2(n)=fzero(@(x0) solve_real_un( RHS,gamma_norm,x0 ),[gamma_norm(int_pos) max(gamma_norm)]);
        catch
            gamma_sol2(n)=0;
        end
    else
        gamma_sol1(n)=0;
        gamma_sol2(n)=0;
    end
end


%%
figure(2)
set(gca,'fontsize',20)
hold on
grid on
gamma_sol1(n_zero_pos)=1;
gamma_sol2(n_zero_pos)=1;

plot(gamma_values/abs(wdia),gamma_sol1,'b','linewidth',2)
plot([min(gamma_values) max(gamma_values)]/abs(wdia),[gamma_eta_wesson gamma_eta_wesson]/gamma_eta_wesson,'k--','linewidth',2)
plot(gamma_values/abs(wdia),gamma_values/gamma_eta_wesson,'r--','linewidth',2)

%  plot(gamma_values/abs(wdia),1.8*gammaDT0_wesson_values/gamma_eta_wesson,'r--','linewidth',2)

legend('\gamma','\gamma=\gamma_{\eta0}','\gamma=\gamma_I')


pause(0.1)
xlim([-1.2 1.2])
ylim([-0.5 4])

xl=xlabel('$\gamma_I / |\omega _{*i}|$');
set(xl,'Interpreter','latex');
yl=ylabel('$\gamma  / \gamma _{\eta0}$');
set(yl,'Interpreter','latex');

%%
close all
% delta_WH_wesson=-gamma_values*tau_H;
LHS=gamma_norm;

for n=1:length(delta_WH_wesson)
    %
    gammaI_ratio=gamma_values(n)/abs(wdia)
    gammaI_ratio_gamma_eta=gamma_values(n)/gamma_eta_wesson
    
    % we need to solve the imaginary part as a fnction of the real part
    % first : to do this we set Re (RHS) = 1
    for r=1:NBREAL
        no_solution=0;
        gamma_hat_range=gamma_norm-i*real_values(r);
        gamma_hat_range=gamma_hat_range';
%         RHS=(gamma_values(n)/8/gamma_eta_wesson)*(gamma_hat_range.^(5/4)).*(gamma_num_values(:,r))./(gamma_den_values(:,r));
        RHS=1/(gamma_values(n)/8/gamma_eta_wesson)*(gamma_hat_range.^(-1.25)).*(gamma_den_values(:,r))./(gamma_num_values(:,r));
        [max_value max_pos]=max(real(RHS));
        max_value=max_value-1;
        [min_value min_pos]=min(real(RHS));
         min_value=min_value-1;
       ini_val=real(RHS(1))-1;
        ini_val=real(RHS(1));
        if sign(max_value*ini_val)<0
            if sign(min_value*ini_val)<0
                int_pos=min(min_pos,max_pos);
            else
                int_pos=max_pos;
            end
        else
            if sign(min_value*ini_val)<0
                int_pos=min_pos;
            else
                %disp('no solutions to real part equation')
                no_solution=1;
            end
        end

        if no_solution==0
            try
                omega_i1(r)=fzero(@(x0) solve_real( RHS,gamma_norm,x0 ),[min(gamma_norm) gamma_norm(int_pos)]);
            catch
                omega_i1(r)=0;
            end
            try
                omega_i2(r)=fzero(@(x0) solve_real( RHS,gamma_norm,x0 ),[gamma_norm(int_pos) max(gamma_norm)]);
            catch
                omega_i2(r)=0;
            end
        else
            omega_i1(r)=0;
            omega_i2(r)=0;
        end
    end
    
    
    % then solve Re(RHS) = 1
    
    for index=1:length(gamma_norm)
        %
        no_solution=0;
        gamma_hat_range=gamma_norm(index)-i*real_values;
%         RHS=(gamma_values(n)/8/gamma_eta_wesson)*(gamma_hat_range.^(5/4)).*(gamma_num_values(index,:))./(gamma_den_values(index,:));
        RHS=1/(gamma_values(n)/8/gamma_eta_wesson)*(gamma_hat_range.^(-1.25)).*(gamma_den_values(index,:))./(gamma_num_values(index,:));
        
        % check if function above or below horizontal axis
        [min_value min_pos]=min(imag(RHS));
        [max_value max_pos]=max(imag(RHS));
        if (min_value>0) || (max_value<0)
            no_solution=1
        end
        if no_solution==0
            
            % central symmetry around origin of the function implies
            % splitting the root searching in two
            % rposfre0 =
            [max_value, max_pos]=max(imag(RHS(1:rposfre0)));
            %         max_value=max_value+real_values(r);
            [min_value min_pos]=min(imag(RHS(1:rposfre0)));
            %         min_value=min_value+real_values(r);
            ini_val=imag(RHS(1));
            if sign(max_value*ini_val)<0
                if sign(min_value*ini_val)<0
                    int_pos1=min(min_pos,max_pos);
                else
                    int_pos1=max_pos;
                end
            else
                if sign(min_value*ini_val)<0
                    int_pos1=min_pos;
                else
                    %                  disp('no solutions to imaginary part equation')
                    %                  min_value
                    %                  max_value
                    %                 no_solution=1;
                end
            end
            if max_value==0
                int_pos1=max_pos+1;
            end
            if min_value==0
                int_pos1=min_pos+1;
            end

            
            [max_value, max_pos]=max(imag(RHS(rposfre0-1:end)));
            max_pos=max_pos+rposfre0;
            [min_value min_pos]=min(imag(RHS(rposfre0-1:end)));
            min_pos=min_pos+rposfre0;
            ini_val=imag(RHS(rposfre0-1));
            if sign(max_value*ini_val)<0
                if sign(min_value*ini_val)<0
                    int_pos2=min(min_pos,max_pos);
                else
                    int_pos2=max_pos;
                end
            else
                if sign(min_value*ini_val)<0
                    int_pos2=min_pos;
                else
                    %                  disp('no solutions to imaginary part equation')
                    %                  min_value
                    %                  max_value
                    %                 no_solution=1;
                end
            end
            
            
            
            try
                omega_r1(index)=fzero(@(x0) solve_imag( RHS,real_values,x0 ),[min(real_values) real_values(int_pos1)]);
            catch
                omega_r1(index)=0;
            end
            try
                omega_r2(index)=fzero(@(x0) solve_imag( RHS,real_values,x0 ),[real_values(int_pos1) real_values(int_pos2)]);
            catch
                omega_r2(index)=0;
            end
            try
                omega_r3(index)=fzero(@(x0) solve_imag( RHS,real_values,x0 ),[real_values(int_pos2) max(real_values)]);
            catch
                omega_r3(index)=0;
            end
        else
            omega_r1(index)=0;
            omega_r2(index)=0;
            omega_r3(index)=0;
        end
        
    end
    
    %     complex_solution_no_real_part(n)=omega_i1(201);
    
    
    %     for index=1:length(gamma_norm)
    %         gamma_hat_range=gamma_norm(index)-i*real_values;
    %         RHS=(gamma_values(n)/8/gamma_eta_wesson)*(gamma_hat_range.^(9/4)).*(gamma_num_values(index,:))./(gamma_den_values(index,:));
    %         omega_r1(r)=fzero(@(x0) solve_real( RHS,gamma_norm(index),real_values,x0 ),[-1.0 0.1]);
    %         omega_r2(r)=fzero(@(x0) solve_real( RHS,gamma_norm(index),real_values,x0 ),[-0.1 1.0]);
    %     end
    
    hold on
    grid on
    plot(omega_r1,gamma_norm)
    plot(omega_r2,gamma_norm)
    plot(omega_r3,gamma_norm)
    plot(real_values,omega_i1,'r')
    plot(real_values,omega_i2,'r')
    pause(0.1)
       
    
    % the intersection(s) of the two curves give operating points
    xl=xlabel('$\omega_r / \gamma _\eta$');
    set(xl,'Interpreter','latex');
    yl=ylabel('$\gamma  / \gamma _\eta$');
    set(yl,'Interpreter','latex');
%     pause
%     close all
    
%     pause
%     close all
%     [residual_eps gamma_transition_pos1]=min(abs(RHS(1:end)-LHS(1:end)));
%     gamma_transition_pos1=gamma_transition_pos1;
%     [residual_eps gamma_transition_pos2]=min(abs(RHS(end:-1:gamma_transition_pos1+1)-LHS(end:-1:gamma_transition_pos1+1)));
%     gamma_transition_pos2=length(RHS)-gamma_transition_pos2+1;
%     if gamma_transition_pos2>3
%         [residual_eps gamma_transition_pos1]=min(abs(RHS(2:gamma_transition_pos2-2)-LHS(2:gamma_transition_pos2-2)));
%         gamma_transition_pos1=gamma_transition_pos1+1;
%     end
%     if gamma_transition_pos1==2
%         gamma_transition_pos1=1;
%     end
%     gamma_transition_value1(n)=gamma_norm(gamma_transition_pos1)*gamma_eta_wesson;
%     gamma_transition_value2(n)=gamma_norm(gamma_transition_pos2)*gamma_eta_wesson;
%
end
    
% hold on
% plot(gamma_values,gamma_transition_value1,'r')
% plot(gamma_values,gamma_transition_value2)
