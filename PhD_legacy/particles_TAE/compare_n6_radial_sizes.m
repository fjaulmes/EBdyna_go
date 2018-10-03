load('../data_tokamak/q_profile.mat')

nTAE=6
mtheta0=nTAE-1;
mtheta1=mtheta0+1;
mtheta2=mtheta1+1;
mtheta3=mtheta2+1;
mTAE=0.5*(mtheta1+mtheta2);

q0=mtheta0/nTAE;
q1=mtheta1/nTAE;
q2=mtheta2/nTAE;
q3=mtheta3/nTAE;
qTAE=(mTAE)/nTAE
r0=interp1(q_initial_profile,radial_r_value_flux,q0)
r1=interp1(q_initial_profile,radial_r_value_flux,q1)
r2=interp1(q_initial_profile,radial_r_value_flux,q2)
r3=interp1(q_initial_profile,radial_r_value_flux,q3)
rTAE=interp1(q_initial_profile,radial_r_value_flux,qTAE);
psiTAE=interp1(q_initial_profile,psi_scale,qTAE);
psiposTAE=interp1(q_initial_profile,1:Nradial,qTAE);
psi1=interp1(q_initial_profile,psi_scale,q1);
psi2=interp1(q_initial_profile,psi_scale,q2);

delta_r_small=sqrt((r2-r1)*0.03/2)
delta_r_normal=sqrt((r2-r1)*0.03)
delta_r_m=sqrt((r2-r1)*0.045)
delta_r_l=sqrt((r2-r1)*0.06)

delta_r_values=[delta_r_small delta_r_normal delta_r_m delta_r_l];

WTAE_values=delta_r_values*0;
gamma_values=delta_r_values*0;
omega_values=delta_r_values*0;

load('TAE_data_n6s.mat')
WTAE_values(1)=WTAE_AVG;
omega_values(1)=omega_TAE;

load('TAE_data_n6.mat')
WTAE_values(2)=WTAE_AVG;
omega_values(2)=omega_TAE;

load('TAE_data_n6m.mat')
WTAE_values(3)=WTAE_AVG;
omega_values(3)=omega_TAE;

load('TAE_data_n6l.mat')
WTAE_values(4)=WTAE_AVG;
omega_values(4)=omega_TAE;


load('fewG_alphas_co_evol_TAEn6small.mat')
gamma_values(1)= gamma_value;

load('fewG_alphas_co_evol_TAEn6standard.mat')
gamma_values(2)=gamma_value; 

load('fewG_alphas_co_evol_TAEn6medium.mat')
gamma_values(3)= gamma_value; 

load('fewG_alphas_co_evol_TAEn6large.mat')
gamma_values(4)=gamma_value; 


figure(1)
set(gca,'fontsize',16)
hold on
grid on

plot(delta_r_values,gamma_values./omega_values,'linewidth',3)
xlabel('$$\Delta r$$','Interpreter','latex')
ylabel('\gamma/\omega_{TAE}')


rho_scale=radial_r_value_flux/max(radial_r_value_flux);

figure(2)
set(gca,'fontsize',16)
hold on
grid on

load('TAE_data_n6s.mat')
load('E0001_n6s.mat')
Epot_radial=squeeze(Epot_map_phi(1,1,:));
plot(rho_scale(pTAE_inf:pTAE_sup),Epot_radial,'g','linewidth',2)

load('E0001_n6.mat')
load('TAE_data_n6.mat')
Epot_radial=squeeze(Epot_map_phi(1,1,:));
plot(rho_scale(pTAE_inf:pTAE_sup),Epot_radial,'b','linewidth',2)

load('E0001_n6m.mat')
load('TAE_data_n6m.mat')
Epot_radial=squeeze(Epot_map_phi(1,1,:));
plot(rho_scale(pTAE_inf:pTAE_sup),Epot_radial,'k','linewidth',2)

load('E0001_n6l.mat')
load('TAE_data_n6l.mat')
Epot_radial=squeeze(Epot_map_phi(1,1,:));
plot(rho_scale(pTAE_inf:pTAE_sup),Epot_radial,'r','linewidth',2)


hl=legend(strcat('$$\Delta r=$$',num2str(delta_r_values(1))),strcat('$$\Delta r=$$',num2str(delta_r_values(2))),...
    strcat('$$\Delta r=$$',num2str(delta_r_values(3))),strcat('$$\Delta r=$$',num2str(delta_r_values(4))));
set(hl,'Interpreter','latex');

xlabel('\rho','Interpreter','latex')
ylabel('\Phi (V)')

xlim([0.15 0.6])
