figure(2)
set(gca,'FontSize',22);

hold on;
grid on;

vpll_bin_size=0.05*1e7;
vpll_range=(-1.6:0.05:1.6)*1e7;
vpll_range=vpll_range';

Nvpll0=histc(vpll0(ENERGY_POPULATION_0),vpll_range-0.5*vpll_bin_size);
Nvpll=post_dist_rescale*histc(alphas_vpll(ENERGY_POPULATION),vpll_range-0.5*vpll_bin_size);

Nvpll0=Nvpll0/max(Nvpll0);
Nvpll=Nvpll/max(Nvpll);
% h = findobj(gca,'Type','patch');
% set(h(2),'FaceColor','b','EdgeColor','k');
% set(h(1),'FaceColor','r');

plot(vpll_range,Nvpll0,'b');
plot(vpll_range,Nvpll,'r');

legend('initial',L2);
title('v_{||}');
xlabel('m/s')


figure(3)
hold on
VPLL_POPULATION=find(alphas_vpll>0);
plot(alphas_pos_x(VPLL_POPULATION),alphas_pos_z(VPLL_POPULATION),'r.')
VPLL_POPULATION=find(alphas_vpll<0);
plot(alphas_pos_x(VPLL_POPULATION),alphas_pos_z(VPLL_POPULATION),'b.')