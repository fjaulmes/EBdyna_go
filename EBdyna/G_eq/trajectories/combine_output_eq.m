function combine_output_eq(ID, type_sim, SKIP_TS, OUTDIR)
%COMBINE_OUTPUT_EQ Combine distributed EBdyna equilibrium outputs.
%
%   combine_output_eq(ID)
%   combine_output_eq(ID, TYPE_SIM)
%   combine_output_eq(ID, TYPE_SIM, SKIP_TS)
%   combine_output_eq(ID, TYPE_SIM, SKIP_TS, OUTDIR)
%
% Inputs
%   ID        Simulation identifier (char or string scalar).
%   type_sim  'full', 'prec', or 'raw'. Default: 'full'.
%   SKIP_TS   Number of initial recorded time stamps to remove. Default: 0.
%   OUTDIR    Directory containing per-process output files. If omitted,
%             the historical geq_data/output_<ID> workflow is used.
%
% Notes
%   - Per-process data are assembled using input.particle_nr.
%   - Fields named loss, loss_wall, nr_profile, and profile_2D are summed.
%   - MATPREC files containing only ejected particles are recorded as
%     intentionally missing process outputs.
%   - The 'raw' branch requires global maps and dim to have been initialized.

    global par maps dim const

    DEFAULT_OUTDIR = '/scratch/project/open-27-30/geq_data';
    OUTPUT_ROOT = '/scratch/project/open-27-30';
    SINGLEFILE_DIR = '../output_singlefile';

    %% Inputs and paths
    if nargin == 0
        ID = '1915611';
        type_sim = 'prec';  % Preserve historical no-argument behavior.
    elseif nargin < 2 || isempty(type_sim)
        type_sim = 'full';
    end
    if nargin < 3 || isempty(SKIP_TS)
        SKIP_TS = 0;
    end
    if nargin < 4 || isempty(OUTDIR)
        OUTDIR = DEFAULT_OUTDIR;
    end

    ID = char(string(ID));
    type_sim = lower(char(string(type_sim)));
    OUTDIR = char(string(OUTDIR));

    validateattributes(SKIP_TS, {'numeric'}, {'scalar', 'integer', 'nonnegative'});
    if ~ismember(type_sim, {'full', 'prec', 'raw'})
        error('combine_output_eq:InvalidType', ...
            'type_sim must be ''full'', ''prec'', or ''raw'' (got ''%s'').', type_sim);
    end

    save_folder = resolve_output_folder(ID, OUTDIR, OUTPUT_ROOT);
    [fnames, fnames_matprec] = find_input_files(save_folder, ID, type_sim);

    if isempty(fnames)
        error('combine_output_eq:NoInputFiles', ...
            'No process output files were found in %s for type ''%s''.', ...
            save_folder, type_sim);
    end

    fprintf('Found %d process output file(s) in %s\n', numel(fnames), save_folder);
    if ~isempty(fnames_matprec)
        fprintf('Found %d MATPREC file(s) to inspect.\n', numel(fnames_matprec));
    end

    %% Load reference process and establish output names
    base = load(fnames{1});
    require_fields(base, {'par', 'input'}, fnames{1});

    par = base.par;
    base_fields = sort(fieldnames(base));
    N_total = base.input.N_total;

    [save_name, save_name_full, output_filename] = ...
        determine_output_names(base.par, ID, type_sim, SINGLEFILE_DIR);

    fprintf('Combined output will be saved in:\n  %s\n', save_name);

    %% Identify MATPREC processes that legitimately have no normal output
    ignored_process_numbers = inspect_matprec_files(fnames_matprec);

    %% Combine all process files
    total = struct();
    total.process_time = [];

    for index = 1:numel(fnames)
        new = load(fnames{index});
        require_fields(new, {'par', 'input'}, fnames{index});

        if SKIP_TS > 0
            new = remove_initial_timestamps(new, SKIP_TS);
        end

        process_number = new.par.PROCESS_NUMBER;
        if ismember(process_number, ignored_process_numbers)
            error('combine_output_eq:ProcessConflict', ...
                ['Process %d is present both as a normal output and as an ', ...
                 'all-ejected MATPREC process.'], process_number);
        end

        if ~isequal(sort(fieldnames(new)), base_fields)
            error('combine_output_eq:InconsistentFields', ...
                'Stored top-level fields differ in %s.', fnames{index});
        end

        new_fields = fieldnames(new);
        for j = 1:numel(new_fields)
            fname = new_fields{j};
            total = combine_fields(total, new, fname, ...
                new.input.particle_nr, new.input.N_job, N_total);
        end

        fprintf('Combined process %d (%d/%d, %.1f%%).\n', ...
            process_number, index, numel(fnames), 100 * index / numel(fnames));
    end

    %% Type-specific post-processing and saving
    switch type_sim
        case 'full'
            total = finalize_full_output(total, base, output_filename);
            save_struct_v73(save_name, total);

        case 'raw'
            append_raw_energy_profile(total, save_name_full, maps, dim, par, const);

        case 'prec'
            save_precession_outputs(total, ID, SINGLEFILE_DIR);
    end

    fprintf('Combined output file created.\n');
end


function save_folder = resolve_output_folder(ID, requested_folder, output_root)
% Preserve the historical geq_data -> output_<ID> workflow safely.

    target_folder = fullfile(output_root, ['output_', ID]);
    requested_folder = strip_trailing_filesep(requested_folder);
    target_folder = strip_trailing_filesep(target_folder);

    if strcmp(requested_folder, target_folder)
        save_folder = target_folder;
    elseif isfolder(requested_folder)
        if ~isfolder(target_folder)
            [ok, msg] = movefile(requested_folder, target_folder);
            if ~ok
                error('combine_output_eq:MoveFailed', ...
                    'Could not move %s to %s: %s', requested_folder, target_folder, msg);
            end
            fprintf('Moved %s -> %s\n', requested_folder, target_folder);
        else
            warning('combine_output_eq:TargetExists', ...
                ['Both %s and %s exist. Using %s and leaving the requested ', ...
                 'folder untouched to avoid accidental nesting.'], ...
                requested_folder, target_folder, target_folder);
        end
        save_folder = target_folder;
    elseif isfolder(target_folder)
        warning('combine_output_eq:FallbackFolder', ...
            'Requested folder %s was not found; using %s.', requested_folder, target_folder);
        save_folder = target_folder;
    else
        error('combine_output_eq:MissingFolder', ...
            'Neither %s nor %s exists.', requested_folder, target_folder);
    end
end


function [fnames, fnames_matprec] = find_input_files(save_folder, ID, type_sim)
    entries = dir(save_folder);
    entries = entries(~[entries.isdir]);

    fnames = {};
    fnames_matprec = {};

    for i = 1:numel(entries)
        filename = entries(i).name;
        lower_name = lower(filename);

        if strcmpi(filename, 'TEMP.mat') || contains(lower_name, 'filepart')
            continue
        end

        is_matprec = contains(lower_name, 'matprec');
        is_raw = contains(filename, [ID, '_raw']);

        if is_matprec
            if contains(filename, ID)
                fnames_matprec{end+1} = fullfile(save_folder, filename); %#ok<AGROW>
            end
            continue
        end

        switch type_sim
            case 'full'
                if ~is_raw
                    fnames{end+1} = fullfile(save_folder, filename); %#ok<AGROW>
                end
            case 'prec'
                fnames{end+1} = fullfile(save_folder, filename); %#ok<AGROW>
            case 'raw'
                if is_raw
                    fnames{end+1} = fullfile(save_folder, filename); %#ok<AGROW>
                end
        end
    end

    fnames = sort(fnames);
    fnames_matprec = sort(fnames_matprec);
end


function [save_name, save_name_full, output_filename] = ...
        determine_output_names(par, ID, type_sim, singlefile_dir)

    require_struct_fields(par, {'SAVENAME'}, 'par');
    save_name_full = truncate_savename_at_token(par.SAVENAME, 'full');
    output_filename = '';

    switch type_sim
        case 'full'
            ensure_directory(singlefile_dir);
            output_filename = ['EBdyna_', ID, '_full.mat'];
            save_name = fullfile(singlefile_dir, output_filename);

        case 'prec'
            ensure_directory(singlefile_dir);
            save_name = fullfile(singlefile_dir, [ID, '_prec.mat']);

        case 'raw'
            % Raw processing appends E_mean to the already-combined full file.
            save_name = save_name_full;
    end
end


function filename = truncate_savename_at_token(savename, token)
    idx = strfind(savename, token);
    if isempty(idx)
        error('combine_output_eq:BadSaveName', ...
            'Could not find token ''%s'' in save name: %s', token, savename);
    end
    idx = idx(1);
    filename = [savename(1:idx + numel(token) - 1), '.mat'];
end


function ignored = inspect_matprec_files(fnames_matprec)
    ignored = [];

    for i = 1:numel(fnames_matprec)
        data = load(fnames_matprec{i});
        require_fields(data, {'ejected', 'par'}, fnames_matprec{i});

        if all(data.ejected(:))
            process_number = data.par.PROCESS_NUMBER;
            warning('combine_output_eq:AllParticlesEjected', ...
                ['No normal output for process %d because every particle ', ...
                 'was ejected in equilibrium.'], process_number);
            ignored(end+1) = process_number; %#ok<AGROW>
        else
            warning('combine_output_eq:UnexpectedMatprec', ...
                'MATPREC file %s contains non-ejected particles.', fnames_matprec{i});
        end
    end

    ignored = unique(ignored);
end


function data = remove_initial_timestamps(data, skip_ts)
    require_fields(data, {'output', 'par'}, 'loaded process data');

    n_ts = data.par.NB_TIME_STAMPS;
    if skip_ts >= n_ts
        error('combine_output_eq:TooManySkippedTimestamps', ...
            'SKIP_TS=%d must be smaller than NB_TIME_STAMPS=%d.', skip_ts, n_ts);
    end

    output_fields = fieldnames(data.output);
    for j = 1:numel(output_fields)
        fname = output_fields{j};
        if strcmp(fname, 'x_ej_next')
            continue
        end

        value = data.output.(fname);
        if ndims(value) == 2 && size(value, 2) == n_ts
            value(:, 1:skip_ts) = [];
        elseif ndims(value) == 3 && size(value, 3) == n_ts
            value(:, :, 1:skip_ts) = [];
        end
        data.output.(fname) = value;
    end

    data.par.time_scale = data.par.time_scale(skip_ts + 1:end);
    data.par.NB_TS_RECORDED = numel(data.par.time_scale);
    data.par.SKIP_TS = skip_ts;
end


function total = finalize_full_output(total, base, output_filename)
    global par const

    if isfield(total, 'output') && ...
            isfield(total.output, 'time_step_loss') && ...
            ~isfield(total.output, 'loss')
        total.output.loss = zeros(base.par.NB_TIME_STEPS, 1);
        for i = 1:base.par.NB_TIME_STEPS
            total.output.loss(i) = sum(total.output.time_step_loss == i);
        end
    end

    total.par.FILENAME = output_filename;

    if isfield(total.par, 'COULOMB_COLL') && total.par.COULOMB_COLL
        if isempty(const)
            require_struct_fields(par, {'paths'}, 'par');
            require_struct_fields(par.paths, {'DATA_FOLDER'}, 'par.paths');
            const = load(fullfile(par.paths.DATA_FOLDER, 'physics_constants.mat'));
        end
        require_struct_fields(const, {'eV'}, 'physics constants');

        if ~isfield(total.output, 'v')
            error('combine_output_eq:MissingVelocity', ...
                'Coulomb-collision post-processing requires total.output.v.');
        end

        v2 = squeeze(sum(total.output.v.^2, 2));
        total.output.Ekin = 0.5 * (total.input.m / const.eV) .* v2;

        vR = squeeze(total.output.v(:, 1, :));
        vR0 = squeeze(total.input.v_ini(:, 1, :));
        birth_matrix = false(size(vR));
        for ts = 1:size(vR, 2)
            if isvector(vR0)
                birth_matrix(:, ts) = (vR(:, ts) ~= vR0(:));
            else
                birth_matrix(:, ts) = (vR(:, ts) ~= vR0(:, ts));
            end
        end
        total.output.birth_matrix = birth_matrix.';
    end

    if isfield(total.output, 'time_step_loss')
        total.output.time_stamp_loss = ceil( ...
            total.output.time_step_loss .* total.par.NB_TIME_STAMPS ./ total.par.NB_TIME_STEPS);
    end
end


function append_raw_energy_profile(total, save_name_full, maps, dim, par, const)
    if isempty(maps) || isempty(dim)
        error('combine_output_eq:RawNeedsMaps', ...
            ['The raw branch requires global maps and dim to be initialized ', ...
             'before calling combine_output_eq.']);
    end
    if isempty(const) || ~isfield(const, 'eV')
        error('combine_output_eq:RawNeedsConstants', ...
            'The raw branch requires global const.eV to be available.');
    end

    speed2 = sum(total.output.v.^2, 2);
    E = 0.5 * (total.input.m / const.eV) .* speed2;
    E_mean = NaN(dim.NB_PSI, par.NB_TIME_STAMPS);

    for t = 1:par.NB_TIME_STAMPS
        X_ind = ((total.output.x(:, 1, t) - dim.R0) * dim.DX_inv) + dim.mid_Xzero;
        Z_ind = ( total.output.x(:, 2, t)            * dim.DZ_inv) + dim.mid_Z;

        psi_norm = ba_interp2(maps(1).psi_norm_XZ, Z_ind, X_ind, 'linear');
        radial_bin = floor(psi_norm);

        for j = 1:dim.NB_PSI
            mask = radial_bin == (j - 1);
            values = E(mask, 1, t);
            values = values(isfinite(values));
            if ~isempty(values)
                E_mean(j, t) = mean(values);
            end
        end

        fprintf('Built raw E profile %d/%d (%.1f%%).\n', ...
            t, par.NB_TIME_STAMPS, 100 * t / par.NB_TIME_STAMPS);
    end

    saved = load(save_name_full, 'output');
    output = saved.output;
    output.E_mean = E_mean;
    save(save_name_full, '-append', 'output');
end


function save_precession_outputs(total, ID, singlefile_dir)
    ensure_directory(singlefile_dir);

    prec = total.prec;
    prec.input = total.input;
    output = total.output;

    filename_prec = fullfile(singlefile_dir, [ID, '_prec']);
    output_filename = fullfile(singlefile_dir, ['output_prec_', ID]);

    save(output_filename, 'output', '-v7.3');
    save(filename_prec, '-struct', 'prec', '-v7.3');
    clean_prec_file(filename_prec, output_filename);

    fprintf('Saved precession files.\n');
end


function save_struct_v73(filename, data)
    parent = fileparts(filename);
    if ~isempty(parent)
        ensure_directory(parent);
    end
    save(filename, '-v7.3', '-struct', 'data');
end


function total = combine_fields(total, new, fname, part_nr, N_job, N_total)
    value = new.(fname);

    if strcmp(fname, 'ind')
        if ~isfield(total, 'ind') || ~isstruct(total.ind)
            total.ind = struct();
        end
        names = fieldnames(new.ind);
        for j = 1:numel(names)
            subname = names{j};
            if ~isfield(total.ind, subname)
                total.ind.(subname) = [];
            end
            total.ind.(subname) = cat(1, total.ind.(subname), new.ind.(subname));
        end
        return
    end

    if isstruct(value)
        if ~isfield(total, fname) || ~isstruct(total.(fname))
            total.(fname) = struct();
        end
        names = fieldnames(value);
        for j = 1:numel(names)
            total.(fname) = combine_fields(total.(fname), value, names{j}, ...
                part_nr, N_job, N_total);
        end
        return
    end

    if strcmp(fname, 'process_time')
        total.process_time(new.par.PROCESS_NUMBER) = value;
        return
    end

    if ismember(fname, {'loss', 'loss_wall', 'nr_profile', 'profile_2D'})
        if ~isfield(total, fname)
            total.(fname) = value;
        else
            total.(fname) = add_accumulated_field(total.(fname), value, fname);
        end
        return
    end

    particle_dim = find_particle_dimension(value, N_job);
    if isempty(particle_dim)
        if ~isfield(total, fname)
            total.(fname) = value;
        end
        return
    end

    value = move_particle_dimension_first(value, particle_dim);

    if ~isfield(total, fname)
        out_size = size(value);
        out_size(1) = N_total;
        if islogical(value)
            total.(fname) = true(out_size);
        else
            total.(fname) = NaN(out_size);
        end
    end

    subs = repmat({':'}, 1, ndims(total.(fname)));
    subs{1} = part_nr;
    total.(fname)(subs{:}) = value;
end


function result = add_accumulated_field(current, incoming, fname)
    if isequal(size(current), size(incoming))
        result = current + incoming;
        return
    end

    if isvector(current) && isvector(incoming) && numel(incoming) <= numel(current)
        result = current;
        result(1:numel(incoming)) = result(1:numel(incoming)) + incoming;
        return
    end

    error('combine_output_eq:AccumulationSizeMismatch', ...
        'Cannot add field ''%s'' with sizes %s and %s.', ...
        fname, mat2str(size(current)), mat2str(size(incoming)));
end


function dim_index = find_particle_dimension(value, N_job)
    if isscalar(value)
        dim_index = [];
        return
    end

    matches = find(size(value) == N_job);

    if isempty(matches)
        dim_index = [];
    elseif matches(1) == 1
        dim_index = 1;
    elseif ismatrix(value) && numel(matches) == 1
        dim_index = matches(1);
    else
        error('combine_output_eq:AmbiguousParticleDimension', ...
            ['Could not uniquely identify the particle dimension for an array ', ...
             'of size %s with N_job=%d.'], mat2str(size(value)), N_job);
    end
end


function value = move_particle_dimension_first(value, particle_dim)
    if particle_dim == 1
        return
    end

    order = 1:ndims(value);
    order([1, particle_dim]) = order([particle_dim, 1]);
    value = permute(value, order);
end


function require_fields(data, fields, source_name)
    for i = 1:numel(fields)
        if ~isfield(data, fields{i})
            error('combine_output_eq:MissingField', ...
                'Missing field ''%s'' in %s.', fields{i}, source_name);
        end
    end
end


function require_struct_fields(data, fields, source_name)
    if ~isstruct(data)
        error('combine_output_eq:ExpectedStruct', '%s must be a struct.', source_name);
    end
    require_fields(data, fields, source_name);
end


function ensure_directory(folder)
    if isempty(folder) || isfolder(folder)
        return
    end
    [ok, msg] = mkdir(folder);
    if ~ok
        error('combine_output_eq:CreateDirectoryFailed', ...
            'Could not create directory %s: %s', folder, msg);
    end
end


function path_out = strip_trailing_filesep(path_in)
    path_out = path_in;
    while numel(path_out) > 1 && (path_out(end) == '/' || path_out(end) == '\\')
        path_out(end) = [];
    end
end
