OMEGAVD_BINS=(-3:0.1:0)*1e4;
OMEGAVD_SCALE=((-3:0.1:-0.05)+0.05)*1e4;

hold on
grid on

[COVD binned]=histc(omega_precess_avg(CO_PASSING),OMEGAVD_BINS);
[COUNTERVD binned]=histc(omega_precess_avg(COUNTER_PASSING),OMEGAVD_BINS);

plot(OMEGAVD_SCALE,COVD(1:end-1),'r')
plot(OMEGAVD_SCALE,COUNTERVD(1:end-1),'b--')

legend('counter current (HFS, stab)','co current (LFS, destab)')
xlabel('\omega_{vD}')