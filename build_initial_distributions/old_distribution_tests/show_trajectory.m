figure(1);

imagesc(scale_X,scale_Z,Fmirror_tot_map');axis xy square
brighten(0.4);
hold on;

plot(Xpos_output,Zpos_output,'k');



figure(2);
subplot(2,1,1);
plot(phipos_output);
grid on;

subplot(2,1,2);
plot(vparallel_output);
grid on;
hold on
plot(vparallel_tilde_output,'b--');
plot(vD_output,'r');
