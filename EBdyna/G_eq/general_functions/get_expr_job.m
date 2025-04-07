function expr=get_expr_job(size_total,NB_PROCESS,PROCESS_NUMBER,expr_or_ind)
%%get_expr_job Returns expression which of the array this core should take 
%care of.
%   Distributes # values in an array (N_array) evenly accross cores.
%   Gives 1 extra to the first couple of nodes if a remainder after
%   devision is present.

if nargin==3
    expr_or_ind='expr';
end

N_array=prod(size_total) %% Total number of cells in the array
if N_array<PROCESS_NUMBER
	error('More processes indicated than number of particles / array values')
end

%% Find indexes for this job
PARTICLES_SPLIT=floor(N_array/NB_PROCESS) % Particles per processor
REM=N_array-NB_PROCESS*PARTICLES_SPLIT   % Remainder of particles

% Take in extra particle in last process if this number of particles are still left
if PROCESS_NUMBER<NB_PROCESS
    N_start=(PROCESS_NUMBER-1)*(PARTICLES_SPLIT)+1 % First particle
    N_end=(PROCESS_NUMBER)*(PARTICLES_SPLIT) % Last particle
else
    N_start=(PROCESS_NUMBER-1)*(PARTICLES_SPLIT)+1 % First particle
    N_end=(PROCESS_NUMBER)*(PARTICLES_SPLIT)+REM % Last particle
end

%% Set indexes to calculate to true
if strcmp(expr_or_ind,'expr')
    expr=false(size_total);
    expr(N_start:N_end)=true;
else
    expr=N_start:N_end;
end

end