
format compact
shot=z0dinput.shot

cdb = cdb_client();
s_neutrals = cdb.get_signal(['neutron_detection_scintillator/PARTICLE_DIAGNOSTICS_RAW:' num2str(shot)]);
s_neutrals.data = abs(s_neutrals.data);
s_neutrals.data=s_neutrals.data/max(s_neutrals.data);
s_neutrals.data=s_neutrals.data*max(zs.ndd);
try 
s_Prad = cdb.get_signal(['Prad/BOLOMETRY:' num2str(shot)]);
end
% time_axis = cdb.get_signal(['Prad_time_axis/BOLOMETRY:' num2str(shot)]);
% s_wdia = cdb.get_signal(['diamagnet_PP_EnergyBT/MAGNETICS:' num2str(shot)]);
try 
s_wdia = cdb.get_signal(['P_dia_norm/MAGNETICS:' num2str(shot)]);
s_wdia.data = cumtrapz(s_wdia.time_axis.data*1e-3,s_wdia.data);
end
try 
s_wdia = cdb.get_signal('diamagnet_PP_Energy/MAGNETICS:15606');
% s_wdia.data = cumtrapz(s_wdia.time_axis.data,s_wdia.data);
end
s_q0 = cdb.get_signal(['q0/EFIT:' num2str(shot)]);
s_q95 = cdb.get_signal(['q95/EFIT:' num2str(shot)]);
s_q = cdb.get_signal(['q/EFIT:' num2str(shot)]);

[B,A]=localButter(4);
if exist('s_wdia')
s_wdia.data = FiltFiltM(B,A,s_wdia.data);
s_wdia.data = FiltFiltM(B,A,s_wdia.data(1:8:end));
s_wdia.time_axis.data = FiltFiltM(B,A,s_wdia.time_axis.data(1:8:end));
s_wdia.data=mean(2*pi*z0dinput.geo.R.*z0dinput.geo.a*2*pi)*s_wdia.data;
end
if exist('s_Prad')
s_Prad.data = FiltFiltM(B,A,s_Prad.data);
s_Prad.data = FiltFiltM(B,A,s_Prad.data(1:8:end));
s_Prad.time_axis.data = FiltFiltM(B,A,s_Prad.time_axis.data(1:8:end));
end



%%
z0plotflux;
clear xlim ylim

%%
tmin=975
tmax=1200

figure; 

subplot(3,1,1)
set(gca,'fontsize',22)
title('neutron rate [count/s]')
hold on
grid on
if ~isempty(s_neutrals.data)
    plot(s_neutrals.time_axis.data,s_neutrals.data,'b','linewidth',2)
end
plot(zs.temps*1000,zs.ndd,'r','linewidth',2)
legend('exp','METIS')
xlim([tmin tmax]);

subplot(3,1,2)
set(gca,'fontsize',22)
title('radiation [W]')
hold on
grid on
if exist('s_Prad')
plot(s_Prad.time_axis.data,s_Prad.data,'b','linewidth',2)
end
plot(zs.temps*1000,zs.prad,'r','linewidth',2)
legend('exp','METIS')
xlim([tmin tmax]);

subplot(3,1,3)
set(gca,'fontsize',22)
title('plasma w (from dia) [J]')
hold on
grid on
plot(s_wdia.time_axis.data,s_wdia.data,'b','linewidth',2)
plot(zs.temps*1000,zs.wdia,'r','linewidth',2)
plot(zs.temps*1000,(real(z0dinput.cons.pnbi)+imag(z0dinput.cons.pnbi))/100,'g','linewidth',2)

legend('exp','METIS','NBI/100')
xlim([tmin tmax]);


%%

figure
subplot(2,1,1)

set(gca,'fontsize',22)
title('safety factor on axis')
hold on
grid on
%
plot(s_q0.time_axis.data,s_q0.data,'b','linewidth',2)
plot(zs.temps*1000,zs.q0,'r','linewidth',2)

legend('exp','METIS')
xlim([tmin tmax]);


subplot(2,1,2)

set(gca,'fontsize',22)
title('safety factor q95 (x=0.95)')
hold on
grid on
%
plot(s_q95.time_axis.data,s_q95.data,'b','linewidth',2)
plot(zs.temps*1000,zs.q95,'r','linewidth',2)

legend('exp','METIS')
xlim([tmin tmax]);


%%


	h=figure('tag','z0plotnbiq');
 
    clf
    set(h,'defaultaxesfontsize',14,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	    'defaultlinelinewidth',1,'color',[1 1 1])
    zplotprof(gca,post.profil0d.temps,post.profil0d.xli,post.profil0d.qjli,'color','r');
    leg = {'METIS'};
    zplotprof(gca,s_q.time_axis.data*1e-3,s_q.axis1.data,s_q.data,'color','b','marker','o','linestyle','none');
    leg{end+1} ='EFIT';
%     zplotprof(gca,post.zerod.temps* ones(1,2),ones(size(post.zerod.temps,1),2),post.zerod.qeff* ones(1,2),'color','g','marker','*');
%     leg{end+1} ='q_{eff}';
%     try
% 	    zplotprof(gca,s_q.time_axis.data,s_q.axis1.data,s_q.data,'color','c','marker','+','linestyle','none');
% 	    leg{end+1} ='EQUI';
%     end 
    legend(leg);
    xlabel('x (normalized radius)');
    ylabel('safety factor');
