clear all
close all

m=4;

NTHETA=10000;
theta=((0:NTHETA)/NTHETA)*2*pi;


ksi_mm1=cos((m-1)*theta);
ksi_m=cos(m*theta);
ksi_mp1=cos((m+1)*theta);
ksi_mp2=cos((m+2)*theta);

OMEGA_VA=1;
OMEGA_VA3=1/3;

NTIME=1000;
sb_tot=zeros(NTIME,1);
m_tot=zeros(NTIME,1);
TAE_tot=zeros(NTIME,1);

figure(1)


for (t_step=0:NTIME)
    t=4*t_step/NTIME;
    ksi_mm1_t=cos((m-1)*(theta-2*pi*t*OMEGA_VA/(m-1)));
    ksi_m_t=cos((m)*(theta-2*pi*t*OMEGA_VA/(m)));
    ksi_mp1_t=cos((m+1)*(theta+2*pi*t*OMEGA_VA/(m+1)));
    ksi_mp2_t=cos((m+2)*(theta+2*pi*t*OMEGA_VA/(m+2)));
    coupling_sb1=ksi_mm1_t+ksi_mp1_t;
    coupling_sb2=ksi_m_t+ksi_mp2_t;
    coupling_main=ksi_m_t+ksi_mp1_t;
    TAE_shape=ksi_mm1_t+ksi_m_t+ksi_mp1_t+ksi_mp2_t;
    sb_term=mean(coupling_sb1.^2);
    m_term=mean(coupling_main.^2);
    TAE_term=mean(TAE_shape.^2);
    sb_tot(t_step+2)=sb_tot(t_step+1)+sb_term;
    m_tot(t_step+2)=m_tot(t_step+1)+m_term;
    TAE_tot(t_step+2)=TAE_tot(t_step+1)+TAE_term;
%     figure(1)
%     subplot(4,1,1)
%     plot(theta,ksi_mm1_t);
%     subplot(4,1,2)
%     plot(theta,ksi_m_t);
%     subplot(4,1,3)
%     plot(theta,ksi_mp1_t);
%     subplot(4,1,4)
%     plot(theta,ksi_mp2_t);
%     pause(0.05)
    
    figure(1)

    subplot(3,1,1)
    plot(theta,coupling_main);
    ylim([-2 2])
    title('main resonance')
    subplot(3,1,2)
    plot(theta,coupling_sb1);
    ylim([-2 2])
    title('inner sideband')
    subplot(3,1,3)
    plot(theta,TAE_shape);
    ylim([-4 4])
    title('TAE strucutre')
end

figure(3)
hold on
plot(m_tot,'k-.')
plot(sb_tot,'b--')
plot(TAE_tot,'r.')
legend('main coupling','inner sideband','TAE')
xlabel('time')
ylabel('mode energy')
