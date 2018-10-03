tmax=time_stamp-1
figure(1)
hold on

for n=1:11
    plot(pphi_output(1:tmax,n)-pphi_output(1,n));
end


figure(2)
hold on

for n=1:11
    plot(Ekin_gc_output(1:tmax,n)-Ekin_gc_output(1,n));
     plot(Ekin_output(1:tmax,n)-Ekin_output(1,n),'r');
end
% 
% 
% figure(3)
% hold on
% 
% for n=1:20
%     plot(pphi_evol(1:tmax,n)-pphi_evol(1,n));
% end