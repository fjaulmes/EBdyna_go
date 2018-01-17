time_stamp=20

[minval Yrs_max]=min(Yrx(2:end-1,time_stamp));
Yrs_max=Yrs_max;
[minval Ycore_max]=min(Ycore(2:end-1,time_stamp));
Ycore_max=Ycore_max;

figure(2);
hold on
plot([0 0],[-1 1],'k--')
plot([-1 1],[0 0],'k--')

plot(Xlinscale(1:Yrs_max),Yrx(1:Yrs_max,time_stamp),'b','linewidth',2)
plot(-Xlinscale(1:Yrs_max),Yrx(1:Yrs_max,time_stamp),'b','linewidth',2)

plot(Xlinscale(1:Ycore_max),Ycore(1:Ycore_max,time_stamp),'r','linewidth',2)
plot(-Xlinscale(1:Ycore_max),Ycore(1:Ycore_max,time_stamp),'r','linewidth',2)

plot([0 Xlinscale(1000)],[Yrx(Yrs_max,time_stamp) Yrx(1000,time_stamp)],'g--','linewidth',2)
plot([0 -Xlinscale(1000)],[Yrx(Yrs_max,time_stamp) Yrx(1000,time_stamp)],'g--','linewidth',2)


xlim([-0.2 0.2])
ylim([-0.23 0.04])
% xlabel('x (m)')
% ylabel('y (m)')
