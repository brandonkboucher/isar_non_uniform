function test_scenario_regression(mode)
% TEST_SCENARIO_REGRESSION  Lock the reconstruction outputs against a baseline.
%
%   Run it as  tests/test_scenario_regression  from the project root, or add
%   tests/ to the path. It resolves everything relative to its own location.
%
%   TEST_SCENARIO_REGRESSION runs every unique scenario in
%   documentation/scenarios.xlsx through ISAR_RUN_SCENARIO with fixed options
%   and a fixed seed, fingerprints the results, and compares them against a
%   stored baseline. It prints one line per scenario and raises an error if
%   anything moved, so it can be dropped into a CI step or run by hand after
%   touching the solvers.
%
%   TEST_SCENARIO_REGRESSION('regenerate') runs the same sweep and writes the
%   baseline instead of comparing. Do this deliberately -- after a change you
%   have decided is correct -- and say so in the commit, because it is the
%   only thing standing between a real regression and a silent one.
%
%   The fingerprint covers, per scenario: the grid dimensions, the norm of the
%   synthesized measurement, the true scatterer positions, and for every
%   algorithm that ran, its RMS error, its recovered positions, the L2 norm
%   and absolute sum of its image, its reflectivity norm, and the size and
%   norm of its refinement history. That is enough to catch a changed
%   derivative, a changed dictionary, or a changed stopping rule, while
%   staying small enough to keep in the repository.
%
%   The sweep runs at a deliberately reduced problem size (see OPTS below) so
%   it finishes in a minute or two. The point is to detect change, not to
%   reproduce publication figures -- isar_testing_v2.m is for that.
%
%   See also ISAR_RUN_SCENARIO, ISAR_DEFAULT_OPTIONS.

    if nargin < 1 || isempty(mode)
        mode = 'compare';
    end
    mode = lower(char(string(mode)));
    if ~ismember(mode, {'compare', 'regenerate'})
        error('test_scenario_regression:mode', ...
            'mode must be ''compare'' or ''regenerate'', got ''%s''', mode);
    end

    % this file lives in tests/, so the project root is one level up
    root = fileparts(fileparts(mfilename('fullpath')));
    addpath(root);
    addpath(fullfile(root, 'tests'));
    addpath(fullfile(root, 'gui'));      % isar_default_options lives here

    % functions/ has subfolders now (omp_and_variants), so add the subtree
    % rather than the one folder. archive/ holds superseded copies of live
    % functions and must never shadow them.
    fn_dirs = strsplit(genpath(fullfile(root, 'functions')), pathsep);
    fn_dirs = fn_dirs(~cellfun(@isempty, fn_dirs));
    fn_dirs = fn_dirs(~contains(lower(fn_dirs), [filesep 'archive']));
    addpath(strjoin(fn_dirs, pathsep));

    scenario_file = fullfile(root, 'documentation', 'scenarios.xlsx');

    % the baseline sits next to this test rather than in results/, so the
    % test and the thing it asserts against move together
    baseline_file = fullfile(root, 'tests', 'regression_baseline.mat');

    % relative tolerance. the sweep is deterministic, so a correct run
    % reproduces bit for bit; the tolerance only absorbs the last bit or two
    % of drift between MATLAB releases
    tol = 1e-9;

    % ---- fixed simulation settings ------------------------------------
    % Held in tests/isar_test_options.m rather than here, and complete rather
    % than a set of overrides, so the sweep is insulated from
    % ISAR_DEFAULT_OPTIONS. A default that moves under the test turns a
    % settings change into an unexplained numeric difference in every case.
    opts = isar_test_options();
    opts.log_fcn      = @(~) [];
    opts.progress_fcn = @(~,~) [];

    % advisory only: a knob that exists but is not pinned would be inherited,
    % which is the failure mode isar_test_options exists to prevent
    unpinned = setdiff(fieldnames(isar_default_options()), fieldnames(opts));
    if ~isempty(unpinned)
        fprintf(['\n  note: %d option(s) not pinned by isar_test_options and ' ...
            'inherited from defaults:\n        %s\n'], ...
            numel(unpinned), strjoin(unpinned(:).', ', '));
    end

    if ~isfile(scenario_file)
        error('test_scenario_regression:noScenarios', ...
            'scenario file not found: %s', scenario_file);
    end

    scenarios = readtable(scenario_file);
    rows      = unique_scenario_rows(scenarios);

    fprintf('\n%s\n', repmat('-', 1, 78));
    fprintf('scenario regression  (%s)\n', mode);
    fprintf('  scenarios : %s (%d rows, %d unique)\n', ...
        'documentation/scenarios.xlsx', height(scenarios), numel(rows));
    fprintf('  baseline  : %s\n', strrep(baseline_file, [root filesep], ''));
    fprintf('%s\n\n', repmat('-', 1, 78));

    % ---- run the sweep -------------------------------------------------
    current = struct('case_id', {}, 'fp', {});
    for ii = 1:numel(rows)
        isc = rows(ii);
        cfg = scenario_config(scenarios, isc);

        t0  = tic;
        out = isar_run_scenario(cfg, opts);
        dt  = toc(t0);

        current(ii).case_id = scenarios.Case_(isc);
        current(ii).fp      = fingerprint(out);    

        fprintf('  case %-3s ran in %5.1f s   %s\n', ...
            string(scenarios.Case_(isc)), dt, describe(cfg));
    end

    % ---- regenerate ----------------------------------------------------
    if strcmp(mode, 'regenerate')
        baseline = struct('created', char(datetime('now')), ...
            'matlab', version, 'opts', strip_handles(opts), ...
            'records', current);

        if ~isfolder(fileparts(baseline_file))
            mkdir(fileparts(baseline_file));
        end
        save(baseline_file, 'baseline');

        fprintf('\nwrote baseline for %d scenarios to %s\n\n', ...
            numel(current), baseline_file);
        return
    end

    % ---- compare -------------------------------------------------------
    if ~isfile(baseline_file)
        error('test_scenario_regression:noBaseline', ...
            ['no baseline at %s\n' ...
             'run  test_scenario_regression(''regenerate'')  to create one'], ...
            baseline_file);
    end

    loaded   = load(baseline_file, 'baseline');
    baseline = loaded.baseline;

    fprintf('\n  baseline written %s (MATLAB %s)\n\n', ...
        baseline.created, baseline.matlab);

    failures = {};

    % a changed option changes every number, so report that instead of
    % drowning the user in per-scenario diffs
    % compare_struct returns [ok, worst, detail]; taking two outputs bound the
    % worst relative difference to opts_detail, so the message below printed
    % blank instead of naming the setting that moved
    [opts_ok, ~, opts_detail] = compare_struct(strip_handles(opts), baseline.opts, 0);
    if ~opts_ok
        failures{end+1} = sprintf('simulation options changed: %s', opts_detail);
    end

    if numel(current) ~= numel(baseline.records)
        failures{end+1} = sprintf('scenario count changed: %d now, %d in baseline', ...
            numel(current), numel(baseline.records));
    else
        fprintf('  %-8s %-8s %-14s %s\n', 'case', 'result', 'worst rel', 'first difference');
        fprintf('  %s\n', repmat('-', 1, 74));

        for ii = 1:numel(current)
            [ok, worst, detail] = compare_struct( ...
                current(ii).fp, baseline.records(ii).fp, tol);

            if ok
                fprintf('  %-8s %-8s %-14.3g\n', ...
                    string(current(ii).case_id), 'PASS', worst);
            else
                fprintf('  %-8s %-8s %-14.3g %s\n', ...
                    string(current(ii).case_id), 'FAIL', worst, detail);
                failures{end+1} = sprintf('case %s: %s', ...
                    string(current(ii).case_id), detail);
            end
        end
    end

    fprintf('\n');
    if isempty(failures)
        fprintf('  all %d scenarios match the baseline (tol %g)\n\n', ...
            numel(current), tol);
        return
    end

    fprintf('  %d failure(s):\n', numel(failures));
    for ii = 1:numel(failures)
        fprintf('    - %s\n', failures{ii});
    end
    fprintf('\n');
    error('test_scenario_regression:mismatch', ...
        '%d scenario(s) differ from the baseline', numel(failures));
end

%% ------------------------------------------------------------------ local

function rows = unique_scenario_rows(scenarios)
% UNIQUE_SCENARIO_ROWS  Row indices with distinct simulation parameters.
%
%   The spreadsheet repeats a configuration once per image formation
%   algorithm, but this sweep runs every algorithm on every scenario, so those
%   rows would be identical work. Case number, algorithm and prediction are
%   descriptive, not parameters, so they are excluded from the comparison.

    descriptive = {'Case_', 'ImageFormationAlgorithm', 'Prediction'};
    param_vars  = setdiff(scenarios.Properties.VariableNames, descriptive, 'stable');

    keys = strings(height(scenarios), 1);
    for i = 1:height(scenarios)
        vals    = table2cell(scenarios(i, param_vars));
        keys(i) = lower(strjoin(cellfun(@(v) char(string(v)), vals, ...
            'UniformOutput', false), '|'));
    end

    [~, rows] = unique(keys, 'stable');
    rows      = rows(:).';
end

function cfg = scenario_config(scenarios, isc)
% SCENARIO_CONFIG  One row of the table as an ISAR_RUN_SCENARIO cfg struct.
%
%   scenarios.xlsx predates the TargetSpacing column, so it is defaulted here
%   rather than failing. If the column is ever added, its value is used.

    cfg = struct( ...
        'NumberOfAmbiguitiesHavingScatterers', ...
            scenarios.NumberOfAmbiguitiesHavingScatterers(isc), ...
        'NumberOfAmbiguitiesInImageFormer', ...
            scenarios.NumberOfAmbiguitiesInImageFormer(isc), ...
        'ScattererLocations',     char(string(scenarios.ScattererLocations(isc))), ...
        'ImageFormerGridDensity', char(string(scenarios.ImageFormerGridDensity(isc))), ...
        'AngleRate',              char(string(scenarios.AngleRate(isc))), ...
        'Noise',                  char(string(scenarios.Noise(isc))), ...
        'TargetSpacing',          'normal');

    if ismember('TargetSpacing', scenarios.Properties.VariableNames)
        cfg.TargetSpacing = char(string(scenarios.TargetSpacing(isc)));
    end
end

function s = describe(cfg)
    s = sprintf('%s / %s / %s / noise %s', ...
        cfg.ScattererLocations, cfg.ImageFormerGridDensity, ...
        cfg.AngleRate, cfg.Noise);
end

function fp = fingerprint(out)
% FINGERPRINT  The numbers this test locks down, as a flat struct.

    fp = struct();
    fp.Nx               = out.sim_config.Nx;
    fp.Ny               = out.sim_config.Ny;
    fp.cross_range_res  = out.sim_config.CrossRangePixelRes_m;
    fp.range_res        = out.sim_config.RangePixelRes_m;
    fp.measurement_l2   = norm(out.y);
    fp.target_locations = out.target_locations;

    for alg = ["omp" "mod_omp" "nomp" "promp" "bp"]
        a = char(alg);
        if ~isfield(out.x_hat, a)
            continue
        end
        r = out.x_hat.(a);

        if isfield(r, 'error'),     fp.([a '_error'])     = r.error;     end
        if isfield(r, 'positions'), fp.([a '_positions']) = r.positions; end

        if isfield(r, 'image')
            fp.([a '_image_l2'])    = norm(r.image(:));
            fp.([a '_image_absum']) = sum(abs(r.image(:)));
        end
        if isfield(r, 'alpha')
            fp.([a '_alpha_l2']) = norm(r.alpha(:));
        end
        if isfield(r, 'p_hat_hist')
            fp.([a '_hist_size']) = size(r.p_hat_hist);
            fp.([a '_hist_l2'])   = norm(r.p_hat_hist(:));
        end
    end
end

function o = strip_handles(opts)
% STRIP_HANDLES  The frozen settings, minus what cannot be compared or saved.
%
%   The options stored in the baseline are isar_test_options' own fields, not
%   the set resolved through ISAR_DEFAULT_OPTIONS -- resolving would couple
%   the baseline back to the file this test is deliberately insulated from,
%   and a default change would fail the options guard even though it cannot
%   reach the sweep. Function handles compare and save meaninglessly, so they
%   come back out.

    o = rmfield(opts, intersect(fieldnames(opts), {'log_fcn', 'progress_fcn'}));
end

function [ok, worst, detail] = compare_struct(now_s, base_s, tol)
% COMPARE_STRUCT  Field-by-field numeric comparison with a relative tolerance.
%
%   Returns the worst relative difference seen and a description of the first
%   field that exceeded TOL. A changed field set or a changed array size is
%   reported directly, since a tolerance cannot speak to either.

    ok    = true;
    worst = 0;
    detail = '';

    fn_now  = fieldnames(now_s);
    fn_base = fieldnames(base_s);

    added   = setdiff(fn_now,  fn_base);
    removed = setdiff(fn_base, fn_now);
    if ~isempty(added) || ~isempty(removed)
        ok = false;
        detail = sprintf('fields added {%s} removed {%s}', ...
            strjoin(added(:).', ', '), strjoin(removed(:).', ', '));
        return
    end

    for i = 1:numel(fn_now)
        f = fn_now{i};
        A = now_s.(f);
        B = base_s.(f);

        if ischar(A) || isstring(A)
            if ~isequal(char(string(A)), char(string(B)))
                ok = false;
                if isempty(detail)
                    detail = sprintf('%s: ''%s'' -> ''%s''', f, ...
                        char(string(B)), char(string(A)));
                end
            end
            continue
        end

        if ~isequal(size(A), size(B))
            ok = false;
            if isempty(detail)
                detail = sprintf('%s: size %s -> %s', f, ...
                    mat2str(size(B)), mat2str(size(A)));
            end
            continue
        end

        A = double(A(:));
        B = double(B(:));

        % scale each difference by the baseline magnitude, with a floor so a
        % value that should be zero is judged on its absolute size
        d = abs(A - B) ./ max(abs(B), 1e-12);
        m = max([d; 0]);

        if m > worst
            worst = m;
        end
        if m > tol
            ok = false;
            if isempty(detail)
                detail = sprintf('%s differs by %.3g (rel)', f, m);
            end
        end
    end
end
