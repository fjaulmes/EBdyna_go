OMEGAVD_BINS=(-28:0.2:28)*1e5;
OMEGAVD_SCALE=((-28:0.2:27.9)+0.1)*1e5;

hold on
grid on

load('initial_NBI60keV_counter_pre_collapse_all.mat')
load('initial_NBI60keV_counter_precession_stats_all.mat')
[COVD binned]=histc(omega_ik_avg(CO_PASSING),OMEGAVD_BINS);
load('initial_NBI60keV_co_pre_collapse_all.mat')
load('initial_NBI60keV_co_precession_stats_all.mat')
[COUNTERVD binned]=histc(omega_ik_avg(COUNTER_PASSING),OMEGAVD_BINS);

plot(OMEGAVD_SCALE,COVD(1:end-1),'r')
plot(OMEGAVD_SCALE,COUNTERVD(1:end-1),'b--')

legend('counter current (HFS, stab)','co current (LFS, destab)')
xlabel('\omega_{d}')
xlim([-2 2]*1e6)