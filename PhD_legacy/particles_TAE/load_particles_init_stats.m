INPUTFILE='initial_MBD_300_precession_stats_all.mat'
INPUTFILE2='initial_MBD_300_pre_collapse_all.mat'

load(INPUTFILE, 'delta_r_avg');
% load(INPUTFILE, 'omega_ae_avg');
load(INPUTFILE, 'r_avg');
load(INPUTFILE, 'vpll_avg');
load(INPUTFILE, 'ALL_PASSING_POP');
load(INPUTFILE, 'CO_PASSING_POP');
load(INPUTFILE2, 'alphas_weight');

Nalphas_simulated=size(pphi_output,2)

r_avg=r_avg(1:Nalphas_simulated);
% omega_ae_avg=omega_ae_avg(1:Nalphas_simulated);
delta_r_avg=delta_r_avg(1:Nalphas_simulated);
ALL_PASSING_POP=ALL_PASSING_POP(1:Nalphas_simulated);
CO_PASSING_POP=CO_PASSING_POP(1:Nalphas_simulated);
vpll_avg=vpll_avg(1:Nalphas_simulated);
alphas_weight=alphas_weight(1:Nalphas_simulated);