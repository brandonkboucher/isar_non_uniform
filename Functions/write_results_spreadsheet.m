function results = write_results_spreadsheet(...
    results, ...        % accumulated results so far; pass [] on the first call
    scenario_rows, ...  % row(s) of the scenarios table this run covers
    x_hat, ...          % the reconstruction struct (.omp / .nomp / .bp / ...)
    filename, ...       % output .xlsx path
    sim_config ...      % OPTIONAL struct of derived config to record as columns
    )

% WRITE_RESULTS_SPREADSHEET  Append one row per scenario and rewrite the output
% spreadsheet.
%
% Each row is one simulated configuration. Scenarios that differ only by image
% formation algorithm are a single row, with one set of result columns per
% algorithm (OMP_Error, NOMP_Error, BP_EQ_Error, ...), so the techniques sit
% side by side for comparison.
%
% Records the scenario configuration exactly as given in the scenarios table,
% any derived simulation settings passed in SIM_CONFIG, and the outputs of
% calculate_reconstruction_error for every algorithm present in X_HAT.
%
% The file is rewritten on every call, so partial results survive an
% interrupted run. Typical use inside the scenario loop:
%
%   results = [];                                   % before the loop
%   ...
%   results = write_results_spreadsheet( ...        % end of each iteration
%       results, ...
%       scenarios(scenario_group == scenario_group(isc), :), ...
%       x_hat, ...
%       'results/simulation_results.xlsx', ...
%       sim_config);
%
% SIM_CONFIG is any struct whose fields are scalars or short strings; each
% field becomes a column. Pass struct() or omit it if you do not want extras.

    if nargin < 5 || isempty(sim_config)
        sim_config = struct();
    end
    if isempty(results)
        results = {};
    end

    % algorithms understood by this reporter, in the order their columns
    % appear. Column 1 is the x_hat field, column 2 is the column-name prefix
    % (must be a valid MATLAB identifier, hence BP_EQ not BP-EQ), column 3 is
    % the ImageFormationAlgorithm spelling used in the scenarios table.
    known_algorithms = { ...
        'omp',   'OMP',   'OMP'; ...
        'nomp',  'NOMP',  'NOMP'; ...
        'bp',    'BP_EQ', 'BP-EQ'; ...
        'lasso', 'LASSO', 'LASSO'; ...
        'sbl',   'SBL',   'SBL'};

    % configuration columns are everything in the scenarios table except the
    % per-row descriptors, so added columns are picked up automatically
    cfg_vars = setdiff(scenario_rows.Properties.VariableNames, ...
        {'Case_', 'ImageFormationAlgorithm', 'Prediction'}, 'stable');

    row = struct();

    % ---- identification --------------------------------------------------
    % all case numbers this simulation covers (the algorithm-only duplicates
    % collapse into this single row)
    row.Cases = char(strjoin(string(scenario_rows.Case_(:)).', ', '));

    % ---- scenario configuration, verbatim from the spreadsheet -----------
    for iv = 1:numel(cfg_vars)
        row.(cfg_vars{iv}) = scalar_value(scenario_rows.(cfg_vars{iv}));
    end

    % ---- derived simulation configuration --------------------------------
    extra = fieldnames(sim_config);
    for ie = 1:numel(extra)
        row.(extra{ie}) = scalar_value(sim_config.(extra{ie}));
    end

    % ---- one block of columns per image formation algorithm --------------
    for ia = 1:size(known_algorithms, 1)

        field  = known_algorithms{ia,1};
        prefix = known_algorithms{ia,2};
        matrix_label = known_algorithms{ia,3};

        if ~isfield(x_hat, field) || ~isfield(x_hat.(field), 'error')
            continue
        end

        res = x_hat.(field);

        % prediction belonging to this algorithm, when the matrix has one
        % prediction = '';
        % if ismember('ImageFormationAlgorithm', scenario_rows.Properties.VariableNames) ...
        %         && ismember('Prediction', scenario_rows.Properties.VariableNames)
        %     for irow = 1:height(scenario_rows)
        %         alg_here = scalar_value(scenario_rows.ImageFormationAlgorithm(irow));
        %         if strcmpi(alg_here, matrix_label)
        %             prediction = scalar_value(scenario_rows.Prediction(irow));
        %             break
        %         end
        %     end
        % end
        % row.([prefix '_Prediction']) = prediction;

        row.([prefix '_Error']) = res.error;
    end

    results{end+1} = row;

    % create the output folder on first write
    outdir = fileparts(filename);
    if ~isempty(outdir) && ~isfolder(outdir)
        mkdir(outdir);
    end

    % rewrite in full so the file is always consistent with `results`
    if isfile(filename)
        delete(filename);
    end
    writetable(rows_to_table(results), filename);

end


function T = rows_to_table(rows)
    % Assemble the accumulated row structs into one table. Rows are normally
    % identical in structure, but an algorithm enabled partway through a run
    % would introduce new fields, so take the union and pad the rest.
    all_fields = {};
    for i = 1:numel(rows)
        f = fieldnames(rows{i});
        all_fields = [all_fields; f(~ismember(f, all_fields))]; %#ok<AGROW>
    end

    % a prototype per field decides how a missing entry is padded
    proto = struct();
    for j = 1:numel(all_fields)
        for i = 1:numel(rows)
            if isfield(rows{i}, all_fields{j})
                proto.(all_fields{j}) = rows{i}.(all_fields{j});
                break
            end
        end
    end

    S = struct();
    for i = 1:numel(rows)
        for j = 1:numel(all_fields)
            fn = all_fields{j};
            if isfield(rows{i}, fn)
                S(i).(fn) = rows{i}.(fn);
            elseif ischar(proto.(fn))
                S(i).(fn) = '';
            else
                S(i).(fn) = NaN;
            end
        end
    end

    T = struct2table(S, 'AsArray', true);
end


function v = scalar_value(x)
    % Reduce one table column / struct field to a single scalar or char row
    % vector that struct2table and writetable can both handle.
    if istable(x)
        x = x{:,1};
    end
    if iscell(x)
        x = x{1};
    elseif ~ischar(x) && numel(x) > 1
        x = x(1);
    end
    if isstring(x) || iscategorical(x)
        v = char(x);
    elseif ischar(x)
        v = x;
    elseif islogical(x)
        v = double(x);
    else
        v = x;
    end
end
