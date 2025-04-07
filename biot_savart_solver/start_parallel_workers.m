function start_parallel_workers(c)

if exist('c','var')==0
    try
        c = str2double(getenv('NUMBER_OF_PROCESSORS')); %For windows only (hypertreading)
    catch
        c = feature('numcores');% For all other systems (only physical cores)
    end
end

p = gcp('nocreate'); % If pool with same amount of workers, do not create new one.
if isempty(p) || p.NumWorkers~=c;
    delete(p);
    feature('numcores'); % Show information of the core
    try
        d = parcluster('local');
        d.NumWorkers=c;
        parpool(d,d.NumWorkers); %See if profile allows more (logical) cores!
    catch
        delete(p);
        parpool(feature('numcores'));    %Use only physical cores
    end
end
