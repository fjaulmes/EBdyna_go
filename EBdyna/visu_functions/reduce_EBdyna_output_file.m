function reduce_EBdyna_output_file(FILENAME)
% clear all

load(FILENAME)

try
load('../../../data_common/wallRZ_CU.mat')
catch
load('../data_common/wallRZ_CU.mat')
end

if ~exist('par')
    error('Please load an EBdyna output first!')
end
if ~exist('hist2d')
    try
    run('../../../startup_EBdyna_go.m')
    catch
    error('Configure by running startup_EBdyna_go...')
    end
end

if ~isfield(par,'FILENAME')
    warning('Filename is undefined! Using a default name to save reduced data...')
    par.FILENAME='EBdyna_reduced_output_default.mat';
end

if isfield(output,'CX_flag')
    output.delta_CX_flag=output.CX_flag(:,end)-output.CX_flag(:,1);
end


output=remove_fields(output,{'v','x','Edep_th','x_ej_next','x_ej_prev','x_ej_prev2','x_ej_prev3','mm','vpll',...
    'pphi_an','vpll_ej','time_step_loss','cum_tor_ang_mom',...
    'x_CXn','x_CXi','CX_rate','CX_flag',...
    'Pdep_e','Pdep_i','Etot_e','Etot_i','Ekin','delta_Ekin','loss','loss_wall','x_gc','CX_NEUTRALS','ndd_val','ndd_cum',...
    'Delta_pphi','dpphi_dt','pphi_kin','psi_value_avg','birth_matrix'})

if ~exist('summary')
    error('Please run a file that has already been summarized')
end
save(FILENAME, 'ejected', 'input', 'output', 'par', 'pfc_losses', 'process_time', 'summary')

disp([' reduced data saved to file  : ' par.FILENAME])
