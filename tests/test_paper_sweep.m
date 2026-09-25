function test_paper_sweep(mode, n_iter)
% TEST_PAPER_SWEEP  Lock the paper's Monte Carlo numbers against a baseline.
%
%   Run it as  tests/test_paper_sweep  from the project root, or add tests/ to
%   the path. It resolves everything relative to its own location.
%
%   TEST_PAPER_SWEEP runs RUN_PAPER_SWEEP with the frozen settings in
%   PAPER_SWEEP_SETTINGS and compares every draw's RMS position error against
%   tests/paper_sweep_baseline.mat, the numbers the paper quotes. It prints
%   the summary table (median, p90, gross failures, timings) and raises an
%   error if any draw moved, so a changed derivative, dictionary, stopping
%   rule or Newton step cannot alter a published figure unnoticed.
%
%   TEST_PAPER_SWEEP('compare', n_iter) compares only the first n_iter draws,
%   for a quicker check; the baseline holds all of them. The full sweep takes
%   a few minutes, so 20 draws is a reasonable smoke test.
%
%   TEST_PAPER_SWEEP('regenerate') runs the full sweep and writes the baseline
%   instead of comparing. Do this deliberately -- after a change you have
%   decided is correct -- and say so in the commit, because it is the only
%   thing standing between a real regression and a silent one.
%
%   Two invariants are asserted besides the per-draw errors, because both are
%   properties of the code rather than of the baseline:
%
%     - "NOMP orig, wrapper" must equal "NOMP (paper)" on every draw, which is
%       what makes nomp_newton a faithful copy of nomp_vec;
%     - no arm may report a missed scatterer, i.e. the miss penalty is never
%       charged (every arm returns all sc.num_of_scatterers positions).
%
%   Timings are printed but never asserted: they are wall clock and depend on
%   the machine. The per-measurement column adds the dictionary build to the
%   run, which is the cost when a dictionary cannot be reused across
%   measurements.
%
%   Every run also writes results/paper_sweep_summary.csv (one row per arm)
%   and results/paper_sweep_errors.csv (one row per draw), via
%   WRITE_SWEEP_CSV, whatever the mode.
%
%   See also RUN_PAPER_SWEEP, PAPER_SWEEP_SETTINGS, TEST_SCENARIO_REGRESSION.

    if nargin < 1 || isempty(mode), mode = 'compare'; end
    mode = lower(char(string(mode)));
    if ~ismember(mode, {'compare', 'regenerate'})
        error('test_paper_sweep:mode', ...
            'mode must be ''compare'' or ''regenerate'', got ''%s''', mode);
    end
    if nargin < 2, n_iter = []; end
    if strcmp(mode, 'regenerate') && ~isempty(n_iter)
        error('test_paper_sweep:partialRegenerate', ...
            'regenerate always runs the full sweep; drop the n_iter argument');
    end

    % this file lives in tests/, so the project root is one level up
    root = fileparts(fileparts(mfilename('fullpath')));
    addpath(root);
    addpath(fullfile(root, 'tests'));

    % the exact-Hessian Newton step and the mod-OMP/NOMP variants that take a
    % Newton routine still live in claude_scratch; promote them into
    % functions/omp_and_variants and drop this line when they settle
    addpath(fullfile(root, 'claude_scratch'));

    % functions/ has subfolders, so add the subtree. archive/ holds superseded
    % copies of live functions and must never shadow them.
    fn_dirs = strsplit(genpath(fullfile(root, 'functions')), pathsep);
    fn_dirs = fn_dirs(~cellfun(@isempty, fn_dirs));
    fn_dirs = fn_dirs(~contains(lower(fn_dirs), [filesep 'archive']));
    addpath(strjoin(fn_dirs, pathsep));

    baseline_file = fullfile(root, 'tests', 'paper_sweep_baseline.mat');

    % the sweep is deterministic, so a correct run reproduces bit for bit;
    % the tolerance only absorbs the last bit or two of drift between releases
    tol = 1e-9;

    % a draw is a gross failure when a scatterer is in the wrong place at all:
    % a wrong ambiguity costs ~5.1 m RMS, a range wrap ~4.6 m, refinement
    % differences stay below a centimetre
    gross = 1.0;

    if strcmp(mode, 'compare')
        if ~isfile(baseline_file)
            error('test_paper_sweep:noBaseline', ...
                ['baseline not found: %s\nrun ' ...
                 'test_paper_sweep(''regenerate'') first'], baseline_file);
        end
        B = load(baseline_file);
        if isempty(n_iter), n_iter = B.res.n_iter; end
        if n_iter > B.res.n_iter
            error('test_paper_sweep:tooManyDraws', ...
                'baseline holds %d draws, asked for %d', B.res.n_iter, n_iter);
        end
    end

    fprintf('\n%s\n', repmat('-', 1, 86));
    fprintf('paper sweep  (%s)\n', mode);
    fprintf('  baseline : %s\n', strrep(baseline_file, [root filesep], ''));
    fprintf('%s\n\n', repmat('-', 1, 86));

    t0  = tic;
    res = run_paper_sweep(n_iter, true);
    fprintf('\n  %d draws in %.1f s\n', res.n_iter, toc(t0));

    % run_paper_sweep already printed the table and wrote the CSVs
    check_invariants(res, gross);

    % ---- regenerate ----------------------------------------------------
    if strcmp(mode, 'regenerate')
        meta = struct('created', char(datetime('now')), 'matlab', version);
        save(baseline_file, 'res', 'meta');
        fprintf('\n  wrote %s\n\n', strrep(baseline_file, [root filesep], ''));
        return
    end

    % ---- compare -------------------------------------------------------
    if ~isequal(cellstr(res.labels), cellstr(B.res.labels))
        error('test_paper_sweep:armsMoved', ...
            ['the arms differ from the baseline''s:\n  now      : %s\n' ...
             '  baseline : %s'], strjoin(cellstr(res.labels), ', '), ...
            strjoin(cellstr(B.res.labels), ', '));
    end

    D = abs(res.E - B.res.E(1:n_iter, :));
    scale = max(abs(B.res.E(1:n_iter, :)), 1);   % relative, floored at 1 m
    worst = max(D ./ scale, [], 1);
    bad   = find(worst > tol);

    fprintf('\n  worst relative difference from the baseline: %.2e\n', max(worst));
    if isempty(bad)
        fprintf('  %d draws x %d arms match the baseline\n\n', n_iter, numel(res.labels));
        return
    end

    for a = bad(:).'
        [~, it] = max(D(:,a));
        fprintf('    %-26s draw %3d (seed %d): %.6f -> %.6f\n', res.labels(a), it, ...
            B.res.base_seed + it - 1, B.res.E(it,a), res.E(it,a));
    end
    error('test_paper_sweep:regression', ...
        '%d of %d arms moved from the baseline', numel(bad), numel(res.labels));
end

% ------------------------------------------------------------------------
function check_invariants(res, gross)
% properties of the code, independent of the recorded numbers

    ia = find(res.labels == "NOMP (paper)");
    ib = find(res.labels == "NOMP orig, wrapper");
    d  = max(abs(res.E(:,ia) - res.E(:,ib)));
    fprintf('\n  nomp_newton reproduces nomp_vec: max |diff| = %.2e\n', d);
    if d > 0
        error('test_paper_sweep:wrapperDrift', ...
            ['"%s" no longer reproduces "%s" (max |diff| %.2e), so it is not ' ...
             'the same algorithm with the Newton step swapped'], ...
            res.labels(ib), res.labels(ia), d);
    end

    % the miss penalty is ~26.7 m, so a charged miss lands near 15 m for one
    % of three scatterers: far above anything a misplacement produces
    miss_floor = res.miss_penalty / sqrt(res.settings.sc.num_of_scatterers) - 1;
    n_missed = sum(res.E(:) >= miss_floor);
    fprintf('  no arm missed a scatterer (errors >= %.1f m): %d\n\n', ...
        miss_floor, n_missed);
    if n_missed > 0
        error('test_paper_sweep:missedScatterer', ...
            ['%d draw(s) scored at or above the miss floor %.1f m, so an arm ' ...
             'returned fewer scatterers than the scene holds'], n_missed, miss_floor);
    end

    if all(res.E(:) < gross)
        warning('test_paper_sweep:noGrossFailures', ...
            'no arm produced a gross failure; the sweep may not be exercising the ambiguity search');
    end
end
