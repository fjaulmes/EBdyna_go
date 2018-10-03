ZHe=1
mHe=mH
pphiTAE=(mHe/eV)*R0*vA3_TAE-ZHe*psiTAE

NON_EJECT_PART=find(~alphas_ejected);
mean_vpll=mean(vparallel_output(1:tmax,NON_EJECT_PART),1);

VPLL_POP=((mean_vpll>0.5*vA3_TAE).*(mean_vpll<1.5*vA3_TAE));

RED_POP=find(VPLL_POP.*(pphi_output(1,NON_EJECT_PART)>=pphiTAE));
BLUE_POP=find(VPLL_POP.*(pphi_output(1,NON_EJECT_PART)<pphiTAE));

omega_output_non_ej=omega_output(1:tmax,NON_EJECT_PART);
pphi_output_non_ej=pphi_output(1:tmax,NON_EJECT_PART);
% omega_output=9.5*theta_output-8*phipos_output;
% omega_output=mod(omega_output,2*pi);

plot(omega_output(1,RED_POP),pphi_output(1,RED_POP),'r.');

figure(1);
hold on;

for t=1:tmax
figure(1)
clf
hold on
plot(omega_output_non_ej(t,BLUE_POP),pphi_output_non_ej(t,BLUE_POP),'b.');
plot(omega_output_non_ej(t,RED_POP),pphi_output_non_ej(t,RED_POP),'r.');
pause(0.5)
end