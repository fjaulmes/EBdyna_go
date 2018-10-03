% reset_data_analysis_environment
close all

clear all


SAVENAME='AUG30382_2p9_mag_energy_evol.mat'

   
    
load('energy_mag_evol_1.mat');
B2_tot_global_evol=B2_tot_evol;
[Evalue end_ts_begin ]=min(B2_tot_evol);
end_ts_begin=end_ts_begin;
    
for process_number=2:15
    disp('----------------------------------------');
    loadname=strcat('energy_mag_evol_',num2str(process_number),'.mat')
    load(loadname);
    [Evalue end_ts_end ]=min(B2_tot_evol(end_ts_begin:end));
    end_ts_end=end_ts_end+end_ts_begin-1
    
    B2_tot_global_evol(end_ts_begin:end_ts_end)=B2_tot_evol(end_ts_begin:end_ts_end);
    end_ts_begin=end_ts_end;

end
    

load('energy_mag_evol_16.mat');
    B2_tot_global_evol(end_ts_begin:end)=B2_tot_evol(end_ts_begin:end);

