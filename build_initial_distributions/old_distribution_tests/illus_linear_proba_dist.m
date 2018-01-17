%%
Nsamples=3e6
ftest=zeros(Nsamples,1);
chi=0.5

for(N=1:Nsamples) 
    ftest(N)=my_rand_linear_dist(1,chi,0,1); 
end


%%
close all

v_values=(0:0.001:1);
pdf=-0.5*v_values+1;
cpdf=-0.25*v_values.^2+v_values;
figure(1);
subplot(2,1,1)
grid on
hold on
set(gca,'fontsize',22)
plot([0 1],[1 chi],'r--','linewidth',3)
plot(v_values,cpdf/0.75,'b','linewidth',3)
legend('f','f_{\rm{cdf}}')
xlabel('v')


subplot(2,1,2)
grid on
hold on
set(gca,'fontsize',22)
hist(ftest,20)
set(get(gca,'child'),'FaceColor','g','EdgeColor','b');
plot([0 1],[200000 200000*chi],'r--','linewidth',3)
xlabel('v')
ylabel('N','Rotation',0)
ylim([0 20]*1e4)
