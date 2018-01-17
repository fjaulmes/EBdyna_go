

energy_dist=f_D_energy_percentage(1,:).*energy_D_range*energy_bin_size;
cum_pdf=energy_dist;

for n=2:length(f_D_energy_percentage(1,:))
cum_pdf(n)=energy_dist(n)+cum_pdf(n-1);
end

% plot(energy_D_range,cum_pdf)
plot(energy_D_range,energy_dist)