function [par,output]=reduce_time_stamps(par,output,new_time)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
narginchk(2,3)
if nargin==2 && isfield(par,'NB_STAMPS_saved')
    new_time=par.NB_STAMPS_saved;
end

if mod(par.NB_TIME_STAMPS,new_time)~=0
    warning('REDUCTION OF TIME STEP FAILED')
    return
end

%% New indexes
step_size=par.NB_TIME_STAMPS/new_time;
new_ind=step_size:step_size:par.NB_TIME_STAMPS;

%% Adjust output
fnames=fieldnames(output);
for i=1:length(fnames)
    t_dim=ndims(output.(fnames{i}));
    if size(output.(fnames{i}),t_dim)~=par.NB_TIME_STAMPS
        continue
    end 
    switch t_dim
        case 2
            output.(fnames{i})=output.(fnames{i})(:,new_ind);
        case 3  
            output.(fnames{i})=output.(fnames{i})(:,:,new_ind);
    end
end
%% Adjuse par
par.TIME_STAMP_PRECISION=step_size*par.TIME_STAMP_PRECISION;
par.time_scale=par.time_scale(new_ind);
par.NB_TIME_STAMPS=new_time;

end
