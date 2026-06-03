function log_msg(level, fmt, varargin)
    global par

    levels = struct('error',0, 'warn',0, 'info',1, 'verbose',2, 'debug',3);

    if ~isfield(par, 'verbose')
        par.verbose = 1;
    end

    if levels.(level) <= par.verbose
        fprintf([fmt '\n'], varargin{:});
    end
end